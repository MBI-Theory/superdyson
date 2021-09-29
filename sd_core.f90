!
!  14 March 2009, SP     - Initial version, from fission of superdyson.f90
!  09 September 2009, SP - Try to exploit frozen cores and identical MO spaces
!  21 May 2010, SP       - Incorporate new GAMESS import routines, incl. 1-e 
!                          integrals. This should allow calculation of completely
!                          general matrix elements, with no restriction on basis
!                          sets or molecular geometries.
!  12 July 2010, SP      - Add an option not to calculate transition dipoles
!                          with the new-style input
!  18 June 2011, SP      - Add option to disable S^2 and normalization checks
!  05 July 2013, SP      - Added support for even newer-style namelist input 
!                          and access to more one-electron matrix elements
!
!  Core routines for dealing with 1-electron operators between
!  general CI wavefunctions.
!
!  WARNING: The format of the input determinant list is perhaps misleading. 
!           The signs of the amplitudes for the determinants are expected
!           for the alpha/beta ordering of the molecular orbitals (ie, the
!           determinants should have all orbitals associated with spin alpha
!           at the beginning, followed by all spin-beta orbitals). This is 
!           the sign choice used by GAMESS CI routines and internally by my 
!           routines here. The commonly used "textbook" sign choice (with
!           orbitals appearing in the fixed canonical order, and alpha/beta
!           spins intemixed) can be obtained by multiplying amplitudes by
!           (-1)**parity, where parity is the number of permutations leading
!           from the alpha/beta to the canonical ordering.
!
module sd_core
  use accuracy
  use timer
  use import_gamess
  use lapack
  use math
  use block_diag
  use block_matrices
  use block_determinant
  use sort_tools
  implicit none
  public               ! The sd_core module is not supposed to lead an
                       ! independent existence - it is always a part of
                       ! a specific integral package. Hence, everything 
                       ! must be public.
  private sd_v2
  !
  character(len=clen), save :: rcsid_sd_core = "$Id: sd_core.f90,v 1.30 2021/09/29 13:43:22 ps Exp $"
  !
  integer, parameter      :: lk = kind(.true.)
  integer, parameter      :: lchar = 1024
  !
  !  Local data for dyson
  !
  integer(ik)             :: verbose = 1        ! Debug level
  integer(ik), parameter  :: ioscr = 10         ! Local use in subroutines
  real(rk), parameter     :: eps_s = 1e-4_rk    ! Tolerance for self-overlap matrix elements
  real(rk)                :: eps_cdet = 1e-5_rk ! Cut-off for determinant weights product
  real(rk)                :: eps_integral = 1e-10_rk
                                                ! Tolerance for blocking overlap and property matrices
  !
  integer(ik)             :: naos_bra = -1      ! Total number of AOs in the "reference"
  integer(ik)             :: naos_ket = -1      ! Total number of AOs in the "ion"
                                                ! -1 means that the value will be obtained from 
                                                ! the MOs input block
  integer(ik)             :: nmobra   = 0       ! Total number of orbitals for the
  integer(ik)             :: nmoket   = 0       ! reference and the ion
  integer(ik)             :: ndetbra  = 0       ! Total number of determinants for the
  integer(ik)             :: ndetket  = 0       ! reference and the ion
  real(rk)                :: sign_bra = 1.0_rk  ! Overall scale factor to apply to the bra wavefunction
  real(rk)                :: sign_ket = 1.0_rk  ! Overall scale factor to apply to the ket wavefunction
  !
  integer(ik)             :: nelbra             ! Number of electrons in the reference
  integer(ik)             :: nelket             ! Number of electrons in the ion
  logical(lk)             :: have_overlaps      ! Overlap integrals are available; MO sets may differ
  logical(lk)             :: have_operator = .true.
                                                ! One-electron integrals are available, calculate
                                                ! corresponding matrix elements
  logical(lk)             :: mode_overlap       ! Same number of electrons in both wfs
  logical(lk)             :: same_basis         ! Is the basis set identical for the ref and ion?
  logical(lk)             :: same_mos           ! MOs used to reprent the two WFs are identical
  integer(ik)             :: frozen_orbitals    ! Number of identical doubly-occupied orbitals in
                                                ! all determinants. 
  logical(lk)             :: dont_check_s2 = .false.
                                                ! True if normalization and S^2 checks are disabled
  character(len=20)       :: task = 'braket'    ! Task to perform; currently one of:
                                                !  'braket' or '1-rdm'
  character(len=20)       :: braket_operator &
                           = 'dipole'           ! Operator to apply for task='bracket'. Can be one of:
                                                !  'none'     - 
                                                !  'kinetic'  - kinetic energy (-0.5 \Nabla)
                                                !  'dipole'   - coordinate expectation
                                                !  'velocity' - velocity \Delta (ie momentum expectation 
                                                !               sans -I*hbar)
                                                !  'r-d/dr'   - r_i (d/d r_j) operator. The i index runs fast
                                                ! Operators below require coordinates of the special centre:
                                                !  '1/r'      - Coulomb nuclear attraction
                                                !  'r'        - Linear attraction
                                                !  'r**2'     - Harmonic attraction
                                                ! Operators below require both coordinates of the centre, and
                                                ! Gaussian exponent:
                                                !  'gauss'    - Gaussian attraction
                                                !  'gauss r'  - Gaussian-damped linear attraction
                                                ! See import_gamess.f90 for the complete list of 
                                                ! operators potentially available in the integrals 
                                                ! package.
                                                ! Regardless of the choice of operator, overlap
                                                ! matrix element will be evaluated as well.
  real(rk)                :: op_centre(3) = 0.  ! Coordinates of the special centre (if any)
  real(rk)                :: op_gauss     = 1.  ! Gaussian exponent (if any)
  integer(ik)             :: op_components = 3  ! Number of components in the operator; must
                                                ! match the choice of operator in braket_operator above
  logical                 :: use_symmetry = .true.
                                                ! Set to true to try to use structure of the overlap
                                                ! matrix for eliminating unnecessary work. There is
                                                ! a slight overhead associated with this, so there
                                                ! may be cases where symmetry use is not beneficial.
  logical                 :: detail_timers = .false.
                                                ! Record timing of low-level routines, potentially with timing
                                                ! overhead exceeding the actual work done
  !
  !  Composition of the determinants is encoded by:
  !   2 = alpha+beta
  !  +1 = alpha
  !  -1 = beta
  !   0 = unoccupied
  !
  integer(ik), allocatable :: occbra(:,:)       ! Occupation numbers (0/1) for reference and
  integer(ik), allocatable :: occket(:,:)       ! ion. The last index is determinant.
  real(rk), allocatable    :: cdetref(:)        ! Weights of the determinants for the reference
  real(rk), allocatable    :: cdetion(:)        ! and the ion
  !
  !  For some of our routines (analyze_spinorbit_integral_blocks, which is in a critical path)
  !  it is advantageous to record sorted lists of orbitals separately. We also need to record
  !  the sorting sequence, so that we can flip back and forth.
  !
  !  The spinbra, spinket, ord_bra, and ord_ket arrays are only needed if symmetry is on.
  !
  integer(ik), allocatable :: scnt_bra  (:,:)   ! Counts of spin-orbitals for the bra determinants
                                                ! First index: 1=alpha, 2=beta
                                                ! Second index: the determinant
  integer(ik), allocatable :: scnt_ket  (:,:)   ! ditto for the ket
  integer(ik), allocatable :: spinbra (:,:,:)   ! First index: up to nelbra spin-orbitals, sorted, grouped by spin
                                                ! Second index: 1=alpha, 2=beta
                                                ! Third index: determinant
  integer(ik), allocatable :: spinket (:,:,:)   ! ditto, for the ket
  integer(ik), allocatable :: ord_bra (:,:,:)   ! Original spin-MO index for each sorted spin-MO in spinbra
  integer(ik), allocatable :: ord_ket (:,:,:)   ! ditto, ket
  !
  real(rk), allocatable    :: cbra (:,:)        ! Reference spatial orbitals
  real(rk), allocatable    :: cket (:,:)        ! Ion spatial orbitals
  real(rk), allocatable    :: sao_ref(:,:)      ! Reference overlap matrix in AO basis
  real(rk), allocatable    :: sao_ion(:,:)      ! Ion overlap matrix in AO basis
  real(rk), allocatable    :: sao  (:,:)        ! Cross-system overlap matrix in AO basis
  real(rk), allocatable    :: dipao(:,:,:)      ! Dipole matrix in AO basis, X/Y/Z along last ind.
  !
  character(len=lchar)     :: comment     = ' ' ! Comment line
  character(len=14), allocatable &
                           :: aolabels(:)       ! Atomic orbital labels, for human readability
  !
  !  File names
  !
  character(len=lchar)     :: file_sao    = ' ' ! File name for the overlaps and (optionally) dipoles
  character(len=lchar)     :: file_cbra   = ' ' ! File name for the reference MOs
  character(len=lchar)     :: file_cket   = ' ' ! File name for the ion MOs
  character(len=lchar)     :: file_detbra = ' ' ! File name for the reference CI coefficients
  character(len=lchar)     :: file_detket = ' ' ! File name for the ion CI coeffficients
  !
  !  The parent MO index is on the left; the ion MO index is on the right.
  !
  real(rk), allocatable    :: sx    (:,:)       ! 1-particle overlaps between MOs of reference and target
  real(rk), allocatable    :: dipx  (:,:,:)     ! Dipole matrix elements between the MOs. Last index is X/Y/Z
                                                ! Could be some other operator, too (see braket_operator),
                                                ! but since the logic does not change, who cares?
  type(bm_structure)       :: sx_ms             ! Block structure of the overlap matrix
  type(bm_structure)       :: dipx_ms           ! Block structure of the property matrices
  !
  !  Stuff for the new-style input:
  !
  logical(lk)              :: own_integrals = .true.
                                                ! We'll calculate our own integrals using GAMESS ckpt.
  type(gam_structure)      :: g_bra             ! Structure and basis description from GAMESS checkpoint
  type(gam_structure)      :: g_ket             ! ditto, the ion
  !
  namelist /sd_v2/ verbose, &
                   task, braket_operator, op_centre, op_gauss, &
                   comment, &
                   nmobra, nmoket, ndetbra, ndetket, eps_cdet, &
                   dont_check_s2, &
                   file_cbra, file_cket, file_detbra, file_detket, &
                   eps_integral, use_symmetry, detail_timers, &
                   sign_bra, sign_ket
  
contains
  subroutine read_input_new
    integer(ik) :: info
    !
    call math_init
    call TimerStart('read_input_new')
    !
    read (input,nml=sd_v2,iostat=info)
    if (info/=0) then
      write (out,"(/'**** ENCOUNTERED ERROR ',i4,' READING INPUT NAMELIST ****'/)") info
    end if
    write (out,"(/'==== Options ====')")
    write (out,nml=sd_v2)
    write (out,"()")
    !
    if (nmobra <=0 .or. nmoket<=0 .or. ndetbra<=0 .or. ndetket<=0 ) then
      write (out,"('ERROR: All of nmobra, nmoket, ndetbra, and ndetket must be given.')")
      stop 'sd_core%read_input_new - missing values with no defaults'
    end if
    !
    call TimerStop('read_input_new')
  end subroutine read_input_new
  !
  subroutine read_core_input
    call math_init
    call TimerStart('read_core_input')
    !
    !  LINE 0: comment, to be added to the output
    !
    read(input,"(a)") comment
    !
    !  LINE 1: naos nmobra nmoket ndetbra ndetket cutoff
    !          If naos is set to -1, this will trigger the new-style
    !          input format. naos of -2 will be internally converted
    !          to -1 together with have_operator = .false.
    !
    read(input,*) naos_bra, nmobra, nmoket, ndetbra, ndetket, eps_cdet
    !
    naos_ket      = naos_bra
    own_integrals = .true.
    dont_check_s2 = .true.
    select case (naos_bra)
      case ( -1) ; have_overlaps = .true. ;  have_operator  = .true. ;  dont_check_s2 = .false.
      case ( -2) ; have_overlaps = .true. ;  have_operator  = .false. ; dont_check_s2 = .false.
      case (-11) ; have_overlaps = .true. ;  have_operator  = .true. ;  
      case (-12) ; have_overlaps = .true. ;  have_operator  = .false. ; 
      case default
        own_integrals = .false.
    end select
    !
    !  LINE 2: name of the file (in GAMESS NPRINT=3 format) containing
    !          the overlap integrals. Not used in the new format.
    !
    if (.not.own_integrals) then
      read(input,*) file_sao
    end if
    !
    !  LINE 3: name of the file (in GAMESS PUNCH format) containing reference MOs
    !          In the old format, everything up to and including the $VEC line
    !          and everyting starting at the $END line should be removed.
    !          In the new format, the file should contain a C1-symmetric geometry
    !          with an inline basis specification (see import_gamess.f90).
    !
    read(input,*) file_cbra
    !
    !  LINE 4: name of the file (in GAMESS PUNCH format) containing ion MOs
    !
    read(input,*) file_cket
    !
    !  LINE 5: name of the file containing determinants for the parent
    !           Composition of the determinants for the parent wavefunction.
    !           Each determinant is defined by:
    !             weight (real)
    !             orbital occupations (nmobra values; 2, 0, +1, or -1)
    !
    read(input,*) file_detbra
    !
    !  LINE 6: Name of the file, containing ion wavefunction. This is not
    !          produced directly by GAMESS; See det_conv.awk for an example
    !          of how to post-process GAMESS output for some interesting
    !          cases.
    !
    read(input,*) file_detket
    call TimerStop('read_core_input')
  end subroutine read_core_input
  !
  subroutine math_init
    real(rk) :: dummy
    !
    if (log10(huge(1._rk))>=4930._rk) then
      ! Quad precision; don't forget "factorial_slack"
      dummy = MathFactorial(1750_ik-5_ik)
      dummy = MathLogFactorial(20000_ik)
    else if (log10(huge(1._rk))>=308._rk) then
      ! Double precision
      dummy = MathFactorial(170_ik-5_ik)
      dummy = MathLogFactorial(10000_ik)
    else
      ! Single precision
      dummy = MathFactorial(34_ik-5_ik)
      dummy = MathLogFactorial(500_ik)
    end if
    dummy = MathDoubleFactorial(80)
  end subroutine math_init
  !
  subroutine initialize_core_data
    integer(ik) :: detref, detion, status
    integer(ik) :: ic, jc
    !
    call TimerStart('initialize_core_data')
    !
    !  Echo core input parameters
    !
    if (verbose>=0) then
      write(out,"(a)") trim(comment)
      write (out,"('          task = ',a)") trim(task)
      write (out,"('      operator = ',a)") trim(braket_operator)
      write (out,"('          naos = ',i6)") naos_bra
      write (out,"('        nmobra = ',i6)") nmobra
      write (out,"('        nmoket = ',i6)") nmoket
      write (out,"('       ndetbra = ',i6)") ndetbra
      write (out,"('       ndetket = ',i6)") ndetket
      write (out,"('      eps_cdet = ',g12.6)") eps_cdet
      write (out,"(' dont_check_s2 = ',l8)") dont_check_s2
      write (out,"(' own_integrals = ',l8)") own_integrals
      !
      if (.not.own_integrals) then
        write (out,"('         AO file = ',a)") trim(file_sao)
      end if
      write (out,"('       Bra MO file = ',a)") trim(file_cbra)
      write (out,"('       Ket MO file = ',a)") trim(file_cket)
      write (out,"('  Bra determinants = ',a)") trim(file_detbra)
      write (out,"('  Ket determinants = ',a)") trim(file_detket)
    end if
    !
    !  If we are using the new format, we'll have to parse checkpoint 
    !  files before doing anything else.
    !
    if (own_integrals) then
      call gamess_load_orbitals(file=trim(file_cbra),structure=g_bra)
      call gamess_load_orbitals(file=trim(file_cket),structure=g_ket)
      naos_bra = g_bra%nbasis
      naos_ket = g_ket%nbasis
      write (out,"('Number of AOs in the reference: ',i6)") naos_bra
      write (out,"('      Number of AOs in the ion: ',i6)") naos_ket
    end if
    !
    !  How many components does our operator have?
    !
    select case (braket_operator)
      case default
        op_components = 1
      case ('none')
        op_components = 0
        have_operator = .false.
      case ('dipole','velocity')
        op_components = 3
      case ('r-d/dr')
        op_components = 9
    end select
    !
    !  Allocate storage for orbital quantities now.
    !
    allocate (occbra(nmobra,ndetbra),cbra(naos_bra,nmobra),cdetref(ndetbra))
    allocate (occket(nmoket,ndetket),cket(naos_ket,nmoket),cdetion(ndetket))
    allocate (sao(naos_bra,naos_ket),dipao(naos_bra,naos_ket,op_components))
    allocate (sao_ref(naos_bra,naos_bra),sao_ion(naos_ket,naos_ket))
    allocate (sx(nmobra,nmoket),dipx(nmobra,nmoket,op_components))
    allocate (aolabels(naos_bra))
    !
    !  1e integrals in the AO basis
    !
    if (own_integrals) then
      have_overlaps = .true.
      call gamess_1e_integrals('AO OVERLAP', sao,         g_bra,g_ket)
      if (have_operator) then 
        select case (braket_operator)
          case default
            write (out,"('initialize_core_data: operator ',a,' is not recognized')") trim(braket_operator)
            stop 'sd_core%initialize_core_data - bad operator'
          case ('kinetic')
            call gamess_1e_integrals('AO KINETIC',   dipao(:,:,1),g_bra,g_ket)
          case ('dipole')
            call gamess_1e_integrals('AO DIPOLE X',  dipao(:,:,1),g_bra,g_ket)
            call gamess_1e_integrals('AO DIPOLE Y',  dipao(:,:,2),g_bra,g_ket)
            call gamess_1e_integrals('AO DIPOLE Z',  dipao(:,:,3),g_bra,g_ket)
          case ('velocity')
            call gamess_1e_integrals('AO D/DX',      dipao(:,:,1),g_bra,g_ket)
            call gamess_1e_integrals('AO D/DY',      dipao(:,:,2),g_bra,g_ket)
            call gamess_1e_integrals('AO D/DZ',      dipao(:,:,3),g_bra,g_ket)
          case ('r-d/dr')
            fill_rddr: do jc=1,3
              do ic=1,3
                call gamess_1e_integrals('AO R-D/DR',dipao(:,:,ic+3*(jc-1)),g_bra,g_ket,op_index=(/ic,jc/))
              end do
            end do fill_rddr
          case ('1/r')
            call gamess_1e_integrals('AO 3C 1/R',    dipao(:,:,1),g_bra,g_ket,op_xyz=op_centre)
          case ('r')
            call gamess_1e_integrals('AO 3C R',      dipao(:,:,1),g_bra,g_ket,op_xyz=op_centre)
          case ('r**2')
            call gamess_1e_integrals('AO 3C R**2',   dipao(:,:,1),g_bra,g_ket,op_xyz=op_centre)
          case ('gauss')
            call gamess_1e_integrals('AO 3C GAUSS',  dipao(:,:,1),g_bra,g_ket,op_xyz=op_centre,op_param=(/op_gauss/))
          case ('gauss r')
            call gamess_1e_integrals('AO 3C GAUSS R',dipao(:,:,1),g_bra,g_ket,op_xyz=op_centre,op_param=(/op_gauss/))
        end select
      end if
      call gamess_1e_integrals('AO OVERLAP', sao_ref,     g_bra,g_bra)
      call gamess_1e_integrals('AO OVERLAP', sao_ion,     g_ket,g_ket)
      aolabels(:) = g_bra%bas_labels
    else  ! .not.own_integrals; can only handle overlaps
      have_overlaps = .true.
      call read_AO_integrals(file_sao,' OVERLAP MATRIX',sao,aolabels,found=have_overlaps)
      sao_ref = sao
      sao_ion = sao
      !
      !  Dipole integrals may or may not be present on input, depending on the
      !  type of caluclation used to generate the output. If no dipoles are
      !  provided, we won't be able to calculate transition matrix elements
      !  or exchange corrections.
      !
      braket_operator = 'dipole'
      have_operator   = .true.
      call read_AO_integrals(file_sao,' X DIPOLE INTEGRALS',dipao(:,:,1),found=have_operator)
      call read_AO_integrals(file_sao,' Y DIPOLE INTEGRALS',dipao(:,:,2),found=have_operator)
      call read_AO_integrals(file_sao,' Z DIPOLE INTEGRALS',dipao(:,:,3),found=have_operator)
    end if ! own_integrals
    !
    !  MO coeffcients. If using the new format, these have been loaded before.
    !  Otherwise, read them now.
    !
    if (own_integrals) then
      if (nmobra>g_bra%nvectors) then
        write (out,"('Need ',i6,' MOs in the reference, but only ',i6,' were found in ',a)") &
               nmobra, g_bra%nvectors, trim(file_cbra)
        stop 'sd_core%read_core_input - not enough ref MOs'
      end if
      if (nmoket>g_ket%nvectors) then
        write (out,"('Need ',i6,' MOs in the ion, but only ',i6,' were found in ',a)") &
               nmoket, g_ket%nvectors, trim(file_cket)
        stop 'sd_core%read_core_input - not enough ion MOs'
      end if
      cbra = g_bra%vectors(:,1:nmobra)
      cket = g_ket%vectors(:,1:nmoket)
    else ! .not.own_integrals
      call read_orbitals(file_cbra,nmobra,cbra)
      call read_orbitals(file_cket,nmoket,cket)
    end if
    !
    !  If orbitals are identical, I'd like to know this.
    !
    call compare_orbitals(same_mos)
    if (.not.same_mos .and. .not.have_overlaps) then
      write (out,"('Reference and ion wavefunctions do not share the MO space,"// &
                 " and no overlap integrals are available')") 
      stop 'insufficient data'
    end if
    !
    !  Ready to read in the reterminants
    !
    open(ioscr,file=file_detbra,action='read',status='old',recl=lchar,iostat=status)
    if (status/=0) then
      write (out,"('Cannot open reference (bra) wavefunction ',a)") trim(file_detbra)
      stop 'file_detbra'
    end if
    read_ref_dets: do detref=1,ndetbra
      read(ioscr,*) cdetref(detref), occbra(:,detref)
    end do read_ref_dets
    close (ioscr)
    open(ioscr,file=file_detket,action='read',status='old',recl=lchar,iostat=status)
    if (status/=0) then
      write (out,"('Cannot open ion (ket) wavefunction ',a)") trim(file_detket)
      stop 'file_detket'
    end if
    read_ion_dets: do detion=1,ndetket
      read(ioscr,*) cdetion(detion), occket(:,detion)
    end do read_ion_dets
    close (ioscr)
    call TimerStop('initialize_core_data')
    !
    !  We can now finally know the number of electrons ...
    !
    nelbra = sum(abs(occbra(:,1)))
    nelket = sum(abs(occket(:,1)))
    !
    !  Precompute sorted spin-orbital lists for symmetry handling
    !
    if (use_symmetry) then
      call TimerStart('initialize_symmetry')
      allocate (scnt_bra(2,ndetbra),spinbra(nelbra,2,ndetbra),ord_bra(nelbra,2,ndetbra))
      allocate (scnt_ket(2,ndetket),spinket(nelket,2,ndetket),ord_ket(nelket,2,ndetket))
      call symmetry_sort_spinorbitals(nelbra,ndetbra,occbra,scnt_bra,spinbra,ord_bra)
      call symmetry_sort_spinorbitals(nelket,ndetket,occket,scnt_ket,spinket,ord_ket)
      call TimerStop('initialize_symmetry')
    end if
  end subroutine initialize_core_data
  !
  subroutine check_input_consistency
    integer(ik) :: detref, detion
    !
    call TimerStart('check_input_consistency')
    !
    !  Check that determinants have matching number of electrons, and 
    !  that ion has one less electron than the reference
    !
    check_ref_dets: do detref=1,ndetbra
      if ( any(occbra(:,detref)/=0 .and. occbra(:,detref)/=2 .and. abs(occbra(:,detref))/=1) ) then
        write (out,"('Reference determinant ',i6,' has invalid occupations')") detref
        stop 'bad ref det'
      end if
      if (sum(abs(occbra(:,detref)))/=nelbra) then
        write (out,"('Reference determinant ',i6,' has invalid electron count')") detref
        stop 'bad ref nel'
      end if
    end do check_ref_dets
    !
    check_ion_dets: do detion=1,ndetket
      if ( any(occket(:,detion)/=0 .and. occket(:,detion)/=2 .and. abs(occket(:,detion))/=1) ) then
        write (out,"('Reference determinant ',i6,' has invalid occupations')") detion
        stop 'bad ref det'
      end if
      if (sum(abs(occket(:,detion)))/=nelket) then
        write (out,"('Reference determinant ',i6,' has invalid electron count')") detion
        stop 'bad ref nel'
      end if
    end do check_ion_dets
    !
    write (out,"('Number of electrons in reference and ion: ',2i6)") nelbra, nelket
    if (nelket  ==nelbra) then
      write (out,"('Overlap mode')")
      mode_overlap = .true.
    else if (nelket+1==nelbra) then
      write (out,"('Dyson mode')")
      mode_overlap = .false.
    else
      stop 'bad electron counts'
    end if
    !
    !  Check normalization of the reference and ions' MOs - this should
    !  guard against any stupid errors in reading the data.
    !
    call check_orthonormality('parent',nmobra,cbra,sao_ref)
    call check_orthonormality('ion   ',nmoket,cket,sao_ion)
    !
    !  Check total wavefunction norm and s^2 expectation value
    !
    if (dont_check_s2) then
      write (out,"(/'WARNING: CHECK FOR S^2 AND NORMALIZATION IS OFF'/)")
    else
      call check_s2('parent',cdetref,occbra)
      call check_s2('ion   ',cdetion,occket)
    end if
    !
    call TimerStop('check_input_consistency')
  end subroutine check_input_consistency
  !
  !  Find the extent of the frozen core, taking into account the MO overlap 
  !  integrals and the determinant populations.
  !
  subroutine locate_frozen_core
    integer(ik) :: ref_docc, ion_docc ! Number of doubly occupied MOs
    integer(ik) :: ndocc, imo
    real(rk)    :: err
    !
    call TimerStart('locate_frozen_core')
    ref_docc = count_doubly_occupied(occbra)
    ion_docc = count_doubly_occupied(occket)
    ndocc    = min(ref_docc,ion_docc)
    scan_overlap: do imo=1,ndocc
      err = abs(sx(imo,imo)-1.0_rk)
      if (err>spacing(1e3_rk)) then
        ndocc = imo - 1
        exit scan_overlap
      end if
    end do scan_overlap
    frozen_orbitals = ndocc
    write (out,"(' Doubly occupied reference orbitals: ',i6)") ref_docc
    write (out,"('       Doubly occupied ion orbitals: ',i6)") ion_docc
    write (out,"('               Frozen core orbitals: ',i6)") frozen_orbitals
    call TimerStop('locate_frozen_core')
  end subroutine locate_frozen_core
  !
  !  Count the number of doubly-occupied MOs in all determinants in the list
  !
  function count_doubly_occupied(dets) result(ndocc)
    integer(ik), intent(in) :: dets(:,:)
    integer(ik)             :: ndocc
    !
    integer(ik) :: nmos, ndets, idet, imo
    !
    nmos  = size(dets,dim=1)
    ndets = size(dets,dim=2)
    ndocc = nmos
    scan_dets: do idet=1, ndets
      scan_mos: do imo=1, nmos
        if (dets(imo,idet)/=2) then
          ndocc = min(ndocc,imo-1)
          exit scan_mos
        end if
      end do scan_mos
    end do scan_dets
  end function count_doubly_occupied
  !
  !  Read one-electron integrals from GAMESS NPRINT=3 output
  !
  subroutine read_AO_integrals(name,section,xao,aolabels,found)
    character(len=*), intent(in)            :: name        ! File containing integrals
    character(len=*), intent(in)            :: section     ! Section header
    real(rk), intent(out)                   :: xao(:,:)    ! Matrix containing the same
    character(len=*), intent(out), optional :: aolabels(:) ! AO labels, for human-readable output
    logical(lk), intent(inout), optional    :: found       ! Set to .FALSE. if section was not found
    !
    character(len=lchar) :: buf
    character(len=lchar) :: errbuf
    integer(ik)          :: status, file_line
    logical(lk)          :: go
    integer(ik)          :: row, column, base_column
    !
    if (verbose>=1) write (out,"('Reading integral section ',a,' from ',a)") trim(section), trim(name)
    file_line = 0
    error_block: do
      errbuf = 'Error opening '//trim(name)//' for reading'
      open(ioscr,file=name,action='read',status='old',recl=lchar,iostat=status)
      if (status/=0) exit error_block
      !
      !  Skip until the header in front of the overlap integrals matrix
      !
      errbuf = 'Error searching for "'//trim(section)//'" in '//trim(name)
      go     = .true.
      skip_header: do while(go)
        buf = ' ' 
        file_line = file_line + 1
        read(ioscr,"(a)",iostat=status) buf
        if (status/=0) exit error_block
        go = (buf/=section)
      end do skip_header
      !
      !  Read the integrals. The whole thing will end with a line
      !  containing a word 'INTEGRALS'
      !
      errbuf = 'Error processing '//trim(section)//' from '//trim(name)
      go     = .true.
      read_integrals: do while(go)
        buf = ' '
        file_line = file_line + 1
        read(ioscr,"(a)",iostat=status) buf
        if (status/=0) exit error_block
        if (buf==' ') cycle read_integrals
        !
        if (index(buf,'INTEGRALS')/=0) then
          go = .false.
          cycle read_integrals
        end if
        !
        if (buf(1:10)==' ') then
          !
          !  This is the header - we'll read the first integer, which
          !  is our base column for the next few lines.
          !
          read(buf,*,iostat=status) base_column
          if (status/=0) exit error_block
          if (base_column<=0 .or. base_column>naos_bra) then  ! Old-style; naos_bra == naos_ket
            errbuf = 'Invalid base column value'
            exit error_block
          end if
        else
          !
          !  Integral lines. The first integer is the orbital index; then
          !  we have labels (to be ignored) and the overlap integrals.
          !
          read(buf,"(i5)",iostat=status) row
          if (status/=0) exit error_block
          if (row<=0 .or. row>naos_bra) then ! Old-style; naos_bra == naos_ket
            errbuf = 'Invalid row value'
            exit error_block
          end if
          if (base_column==1 .and. present(aolabels)) then
            read(buf,"(5x,a11)",iostat=status) aolabels(row)
          end if
          read(buf,"(15x,5f11.6)",iostat=status) &
            (xao(row,column),column=base_column,min(row,base_column+4))
          if (status/=0) exit error_block
          !
          !  Add the transpose
          !
          xao(base_column:min(row,base_column+4),row) = xao(row,base_column:min(row,base_column+4))
        end if
      end do read_integrals
      !
      !  Close the file, and we are done
      !
      errbuf = 'Error closing '//trim(name)
      close(ioscr,iostat=status)
      if (status/=0) exit error_block
      return
    end do error_block
    write (out,"(' Line ',i10,': ',a)") file_line, trim(errbuf)
    if (present(found)) then
      found = .false.
    else
      stop 'read_AO_integrals'
    end if
  end subroutine read_AO_integrals
  !
  !  Load a set of molecular orbitals.
  !
  subroutine read_orbitals(name,nmos,c)
    character(len=*), intent(in) :: name   ! Name of the file
    integer(ik), intent(in)      :: nmos   ! Number of orbitals to read
    real(rk), intent(out)        :: c(:,:) ! Space for the orbitals
    !
    integer(ik)          :: status, file_line
    character(len=lchar) :: errbuf
    integer(ik)          :: imo, minao, maxao
    integer(ik)          :: chk_mo, chk_line, line
    !
    if (verbose>=1) write (out,"('Reading orbitals from ',a)") trim(name)
    file_line = 0
    error_block: do
      errbuf = 'Error opening '//trim(name)//' for reading'
      open(ioscr,file=name,action='read',status='old',recl=lchar,iostat=status)
      if (status/=0) exit error_block
      !
      errbuf = 'Error reading orbitals from '//trim(name)
      mos_loop: do imo=1,nmos
        !
        ! The orbitals are split over a number of lines, in a rather 
        ! ugly format - so we'll have to jump through a few hoops.
        !
        minao = 1
        line  = 1
        line_loop: do while(minao<=naos_bra)  ! Note that with the old input, naos_bra and naos_ket are identical
          maxao = min(minao+4,naos_bra)
          file_line = file_line + 1
          read(ioscr,"(i2,i3,5e15.0)",iostat=status) chk_mo, chk_line, c(minao:maxao,imo)
          if (status/=0) exit error_block
          if (mod(chk_mo,100)/=mod(imo,100) .or. chk_line/=line) then
            write (errbuf,"('In file ',a,' chk_mo = ',i5,' imo = ',i5,' chk_line = ',i5,' line = ',i5)") &
                   trim(name), chk_mo, imo, chk_line, line
            exit error_block
          end if
          minao = maxao + 1
          line  = line  + 1
        end do line_loop
      end do mos_loop
      !
      errbuf = 'Error closing '//trim(name)
      close(ioscr,iostat=status)
      return
    end do error_block
    write (out,"(' Line ',i10,': ',a)") file_line, trim(errbuf)
    stop 'read_orbitals'
  end subroutine read_orbitals
  !
  !  Compare orbital sets for being identical - this permits many short-cuts
  !
  subroutine compare_orbitals(same)
    logical(lk), intent(out) :: same ! Orbital sets are identical
    !
    integer(ik) :: nmo, iat, ish, is, p1, p2, p1i, p2i
    real(rk)    :: diff, maxc
    !
    same = .false.
    same_basis = .false.
    check_orbitals: do
      !
      !  If using the new input, we have to make sure that geometries and the
      !  number of AOs are the same.
      !
      if (own_integrals) then
        write (out,"('           Reference atoms = ',i8)") g_bra%natoms
        write (out,"('                 Ion atoms = ',i8)") g_ket%natoms
        if (g_bra%natoms /= g_ket%natoms) exit check_orbitals
        !
        scan_atoms: do iat=1,g_bra%natoms
          diff = maxval(abs(g_bra%atoms(iat)%xyz-g_ket%atoms(iat)%xyz))
          if (diff>spacing(100._rk)) then
            write (out,"('Reference atom ',i4,' is at ',3f16.8)") iat, g_bra%atoms(iat)%xyz
            write (out,"('      Ion atom ',i4,' is at ',3f16.8)") iat, g_ket%atoms(iat)%xyz
            exit check_orbitals
          end if
          if (g_bra%atoms(iat)%nshell/=g_ket%atoms(iat)%nshell) then
            write (out,"('Reference atom ',i4,' has ',i4,' basis shells')") iat, g_bra%atoms(iat)%nshell
            write (out,"('      Ion atom ',i4,' has ',i4,' basis shells')") iat, g_ket%atoms(iat)%nshell
            exit check_orbitals
          end if
          scan_shells: do is=1,g_bra%atoms(iat)%nshell
            if ( g_bra%atoms(iat)%sh_l(is) /= g_ket%atoms(iat)%sh_l(is) ) then
              write (out,"('Reference atom ',i4,' shell ',i3,' has angular momentum ',i2)") iat, is, g_bra%atoms(iat)%sh_l(is)
              write (out,"('      Ion atom ',i4,' shell ',i3,' has angular momentum ',i2)") iat, is, g_ket%atoms(iat)%sh_l(is)
              exit check_orbitals
            end if
            p1  = g_bra%atoms(iat)%sh_p(is)
            p1i = g_ket%atoms(iat)%sh_p(is)
            p2  = g_bra%atoms(iat)%sh_p(is+1) - 1
            p2i = g_ket%atoms(iat)%sh_p(is+1) - 1
            if ( p1/=p1i .or. p2/=p2i ) then
              write (out,"('Reference atom ',i4,' shell ',i4,' spans primitives ',i4,' - ',i4)") iat, is, p1,  p2
              write (out,"('      Ion atom ',i4,' shell ',i4,' spans primitives ',i4,' - ',i4)") iat, is, p1i, p2i
              exit check_orbitals
            end if
            maxc = max(maxval(abs(g_bra%atoms(iat)%p_zet(p1:p2))),maxval(abs(g_ket%atoms(iat)%p_zet(p1:p2))))
            diff = maxval(abs(g_bra%atoms(iat)%p_zet(p1:p2)-g_ket%atoms(iat)%p_zet(p1:p2)))
            if (diff>100._rk*maxc) then
              write (out,"('Reference and ion atoms ',i4,' shell ',i4,' have different primitive exponents')") iat, ish
              exit check_orbitals
            end if
            maxc = max(maxval(abs(g_bra%atoms(iat)%p_c(p1:p2))),maxval(abs(g_ket%atoms(iat)%p_c(p1:p2))))
            diff = maxval(abs(g_bra%atoms(iat)%p_c(p1:p2)-g_ket%atoms(iat)%p_c(p1:p2)))
            if (diff>100._rk*maxc) then
              write (out,"('Reference and ion atoms ',i4,' shell ',i4,' have different contraction weights')") iat, ish
              exit check_orbitals
            end if
          end do scan_shells
        end do scan_atoms
      end if
      !
      same_basis = .true.
      nmo  = min(nmoket,nmobra)
      maxc = max(maxval(abs(cbra)),maxval(abs(cket)))
      diff = maxval(abs(cbra(:,:nmo)-cket(:,:nmo)))
      same = (diff <= 100._rk * spacing(maxc)) 
      write (out,"('Comparing ',i6,' lowest MOs out of ',i6,'/',i6)") nmo, nmobra, nmoket
      write (out,"('  Max absolute coefficient = ',g12.6)") maxc
      write (out,"('Max coefficient difference = ',g12.6)") diff
      exit check_orbitals
    end do check_orbitals
    write (out,"('        The two sets match = ',l8)") same
  end subroutine compare_orbitals
  !
  !  Calculate and check overlaps between MOs - we must get a unit matrix
  !  AO overlaps are in sao(:,:).
  !
  subroutine check_orthonormality(title,nmos,mos,sao)
    character(len=*), intent(in) :: title    ! For error messages 
    integer(ik), intent(in)      :: nmos     ! Number of MOs
    real(rk), intent(in)         :: mos(:,:) ! MO coefficients
    real(rk), intent(in)         :: sao(:,:) ! Overlap matrix (possibly ion or ref-specific)
    !
    integer(ik) :: lmo, rmo ! MO indices
    real(rk)    :: ov       ! Calculated MO overlap
    real(rk)    :: refval   ! Expected MO overlap
    !
    if (.not.have_overlaps) then
      write (out,"('Atomic overlap integrals not present; unable to check '" // &
                 ",a,' orbitals for orthonormality')") trim(title)
      return
    end if
    !
    left_scan: do lmo=1,nmos
      right_scan: do rmo=1,nmos
        ov = dot_product(mos(:,lmo),matmul(sao,mos(:,rmo)))
        !
        refval = 0._rk
        if (lmo==rmo) refval = 1._rk
        !
        if (abs(ov-refval)>eps_s) then
          write (out,"(5x,'warning: ',a,' <',i5,'|',i5,'> = ',g16.9,' (expected ',g16.9,')')") &
                 trim(title), lmo, rmo, ov, refval
        end if
      end do right_scan
    end do left_scan
    !
  end subroutine check_orthonormality
  !
  !  Check norm and S^2
  !
  subroutine check_s2(title,c,d)
    character(len=*), intent(in) :: title  ! Name of the w.f.
    real(rk), intent(in)         :: c(  :) ! Amplitudes
    integer(ik), intent(in)      :: d(:,:) ! Determinants
    !
    real(rk)    :: ov, s2, s2exact
    integer(ik) :: is
    !
    call calculate_s2(c,d,ov,s2)
    !
    write (out,"(/1x,'         ',a,' norm is: ',f12.8 )") trim(title), ov
    write (out,"( 1x,a,' S^2 expectation is: ',f12.8, '(',f12.8,')'/)") trim(title), s2, s2/ov
    !
    if (abs(ov-1._rk)>0.03_rk) then
      write (out,"('sd_core%check_s2 - unnormalized wavefunction')")
      stop 'sd_core%check_s2 - unnormalized wavefunction'
    end if
    !
    scan_spins: do is=0,10
      s2exact = 0.25_rk * is * (is+2)
      if (abs(s2exact-s2/ov)<0.03_rk) return
    end do scan_spins
    stop 'sd_core%check_s2 - wavefunction is not an S^2 eigenfunction'
  end subroutine check_s2
  !
  !  Give user some idea as to how similar the MOs are between parent
  !  and the ion.
  !
  subroutine report_mo_correlations(smo)
    real(rk), intent(in) :: smo(:,:) ! Cross-system overlap integrals, MO basis
    !
    integer(ik) :: refmo, simmo
    !
    write (out,"(/t5,'Reference - ion MO correlations')")
    write (out,"(2x,a7,2x,a7,2x,a12)") 'Ref. MO', 'Ion MO', 'Overlap', &
                                       '-------', '------', '-------'
    mo_ref: do refmo=1,nmobra
      simmo = maxloc(abs(smo(refmo,:)),dim=1)
      write (out,"(2x,i7,2x,i7,2x,f12.5)") refmo, simmo, smo(refmo,simmo)
    end do mo_ref
    write (out,"()")
  end subroutine report_mo_correlations
  !
  !  Construct a one-particle operator matrix elements in cross-state 
  !  molecular orbital basis. Reference state is on the left; target 
  !  is on the right; all alpha-spin appears before all beta-spin
  !  spin-orbitals.
  !
  subroutine mat_ao_to_mo(title,xmo,xao)
    character(len=*), intent(in) :: title    ! Title; identification only
    real(rk), intent(out)        :: xmo(:,:) ! Matrix elements in the MO basis
    real(rk), intent(in)         :: xao(:,:) ! Matrix elements in the AO basis
    !
    integer(ik) :: refmo, ionmo ! Orbitals of the reference and ion
    real(rk)    :: val
    !
    call TimerStart('mat_ao_to_mo')
    !
    !  This is a very inefficient way of making the transformation.
    !  Much better perfomance can be obtained with MATMUL and a buffer or DGEMM
    !
    if (verbose>=2) write (out,"(/'Cross-state ',a,' elements')") trim(title)
    mo_ref: do refmo=1,nmobra
      mo_ion: do ionmo=1,nmoket
        val = dot_product(cbra(:,refmo),matmul(xao,cket(:,ionmo)))
        if (verbose>=2) then
          write (out,"(5x,a,' <',i5,'|',i5,'> = ',f25.12)") trim(title), refmo, ionmo, val
        end if
        xmo(refmo,ionmo) = val
      end do mo_ion
    end do mo_ref
    if (verbose>=2) write (out,"()")
    !
    call TimerStop('mat_ao_to_mo')
  end subroutine mat_ao_to_mo
  !
  !  Transform 1-electron quantity back from MO to the AO basis.
  !  Should reverse the action mat_ao_to_mo
  !
  subroutine mat_mo_to_ao(title,xmo,xao)
    character(len=*), intent(in) :: title    ! Title; identification only
    real(rk), intent(in)         :: xmo(:,:) ! Matrix elements in the MO basis
    real(rk), intent(out)        :: xao(:,:) ! Matrix elements in the AO basis
    !
    integer(ik) :: refao, ionao ! Atomic orbitals of the reference and the ion
    real(rk)    :: val
    !
    call TimerStart('mat_mo_to_ao')
    !
    !  This is a very inefficient way of making the transformation.
    !  Much better perfomance can be obtained with MATMUL and a buffer or DGEMM
    !
    if (verbose>=2) write (out,"(/'Cross-state ',a,' elements')") trim(title)
    ao_ref: do refao=1,naos_bra
      ao_ion: do ionao=1,naos_ket
        val = dot_product(cbra(refao,:),matmul(xmo,cket(ionao,:)))
        if (verbose>=2) then
          write (out,"(5x,a,' <',i5,'|',i5,'> = ',f25.12)") trim(title), refao, ionao, val
        end if
        xao(refao,ionao) = val
      end do ao_ion
    end do ao_ref
    if (verbose>=2) write (out,"()")
    !
    call TimerStop('mat_mo_to_ao')
  end subroutine mat_mo_to_ao
  !
  !  Initialize unit matrix
  !
  subroutine mat_unit(m)
    real(rk), intent(out) :: m(:,:) 
    integer(ik)           :: i
    m = 0._rk
    do i=1,min(size(m,dim=1),size(m,dim=2))
      m(i,i) = 1._rk
    end do
  end subroutine mat_unit
  !
  !  Calculate the distance between two determinants
  !
  function determinant_distance(ref,ion) result(n)
    integer(ik), intent(in) :: ref(:), ion(:)
    integer(ik)             :: n
    !
    integer(ik) :: len_r, len_i, len_c, idet
    integer(ik) :: nr, ni
    !
    len_r = size(ref) 
    len_i = size(ion)
    len_c = min(len_r,len_i)
    n     = 0
    sum_common: do idet=1,len_c
      nr = ref(idet)
      ni = ion(idet)
      if (nr==ni) cycle sum_common
      !
      if (abs(nr)/=1) then
        n = n + abs(nr - abs(ni))
      else if (abs(ni)/=1) then
        n = n + abs(ni - abs(nr))
      else
        n = n + 2
      end if
    end do sum_common
    if (len_r>len_c) n = n + sum(abs(ref(len_c+1:)))
    if (len_i>len_c) n = n + sum(abs(ion(len_c+1:)))
  end function determinant_distance
  !
  !  Prepare translation table for the orbital indices in a determinant 
  !  and the full space. Should be used as a companion to fill_spinorbit_integrals
  !
  subroutine fill_spinorbit_indices(occ,ind)
    integer(ik), intent(in)  :: occ(:)    ! Occupation table
    integer(ik), intent(out) :: ind(:)    ! Mapping between spin-orbital indices in the 
                                          ! determinant, and spatial MOs in the full
                                          ! orbital space
    integer(ik) :: mo   ! Spatial MO indices
    integer(ik) :: spin ! Spin+space MO indices
    !
    !  Spin-alpha orbitals are first
    !
    spin = 0
    alpha: do mo=1,size(occ)
      if (occ(mo)/=1 .and. occ(mo)/=2) cycle alpha
      !
      spin = spin + 1
      ind(spin) = mo
    end do alpha
    !
    !  Spin-beta orbitals - map to the same set of spatial orbitals
    !
    beta: do mo=1,size(occ)
      if (occ(mo)/=-1 .and. occ(mo)/=2) cycle beta
      !
      spin = spin + 1
      ind(spin) = -mo
    end do beta
    if (spin/=size(ind)) then
      write (out,"('fill_spinorbit_indices: Counted ',i8,' electrons, but expected ',i8)") spin, size(ind)
      stop 'sd_core%fill_spinorbit_indices - bad index size'
    end if
  end subroutine fill_spinorbit_indices
  !
  !  Calculate 1-particle overlap for all spin-orbitals in a given
  !  pair of determinants. If dipoles are present, simultaneously 
  !  arrange dipole integrals in a matching sequence.
  !
  subroutine fill_spinorbit_integrals(occbra,occket,adet,amo)
    integer(ik), intent(in) :: occbra(:)            ! Parent determinant
    integer(ik), intent(in) :: occket(:)            ! Ion determinant
    real(rk), intent(out)   :: adet(:,:)            ! Integral matrix for a given determinant
    real(rk), intent(in)    :: amo (:,:)            ! Integral matrix for all MOs
    !
    integer(ik) :: moref, moion     ! Spatial MO indices
    integer(ik) :: spinref, spinion ! Spin+space MO indices
    integer(ik) :: nalpharef        ! Number of spin-alpha electrons in reference
    !
    if (detail_timers) call TimerStart('fill_spinorbit_integrals')
    adet = 0._rk
    !
    !  Alpha-alpha overlaps - upper left corner of sdet
    !
    spinion = 0
    spinref = 0
    ion_alpha: do moion=1,nmoket
      if (occket(moion)/=1 .and. occket(moion)/=2) cycle ion_alpha
      !
      spinion = spinion + 1
      spinref = 0
      ref_alpha: do moref=1,nmobra
        if (occbra(moref)/=1 .and. occbra(moref)/=2) cycle ref_alpha
        !
        spinref = spinref + 1
        adet(spinref,spinion) = amo(moref,moion)
      end do ref_alpha
    end do ion_alpha
    nalpharef = spinref
    !
    !  Beta-beta overlaps - bottom right corner of sdet
    !
    ion_beta: do moion=1,nmoket
      if (occket(moion)/=-1 .and. occket(moion)/=2) cycle ion_beta
      !
      spinion = spinion + 1
      spinref = nalpharef
      ref_beta: do moref=1,nmobra
        if (occbra(moref)/=-1 .and. occbra(moref)/=2) cycle ref_beta
        !
        spinref = spinref + 1
        adet(spinref,spinion) = amo(moref,moion)
      end do ref_beta
    end do ion_beta
    if (detail_timers) call TimerStop('fill_spinorbit_integrals')
  end subroutine fill_spinorbit_integrals
  !
  subroutine symmetry_sort_spinorbitals(nel,ndet,occ,scnt,spin,ord)
    integer(ik), intent(in)  :: nel         ! Number of occupied spin-orbitals
    integer(ik), intent(in)  :: ndet        ! Number of determinants
    integer(ik), intent(in)  :: occ (:,:)   ! Determinant list, "textbook" (2/-1/+1) orbilal order
    integer(ik), intent(out) :: scnt  (:,:) ! Counts of alpha/beta spin-orbitals per determinant
    integer(ik), intent(out) :: spin(:,:,:) ! Sorted orbital lists, separately for alpha and beta
    integer(ik), intent(out) :: ord (:,:,:) ! Sorting sequence 
    !
    integer(ik) :: idet
    integer(ik) :: spintmp(nel), ordtmp(nel)! Temporary arrays, for combined alpha and beta orbital lists
    !
    !$omp parallel do default(none) private(idet,spintmp,ordtmp) shared(ndet,occ,scnt,spin,ord)
    convert_determinants: do idet=1,ndet
      call fill_spinorbit_indices(occ(:,idet),spintmp)
      call sort(spintmp,ordtmp)
      call copy_plus (spintmp,ordtmp,scnt(1,idet),spin(:,1,idet),ord(:,1,idet))
      call copy_minus(spintmp,ordtmp,scnt(2,idet),spin(:,2,idet),ord(:,2,idet))
    end do convert_determinants
    !$omp end parallel do
    contains
    !
    subroutine copy_plus(spin_all,ord_all,cnt,spin,ord)
      integer(ik), intent(in)  :: spin_all(:) ! Alpha (positive) and beta (negative) spin-orbitals, sorted
      integer(ik), intent(in)  :: ord_all (:) ! Sorting order, matching spin_all
      integer(ik), intent(out) :: cnt         ! Number of positive spin-orbitals
      integer(ik), intent(out) :: spin(:)     ! Positive spin-orbitals, sorted
      integer(ik), intent(out) :: ord (:)     ! Matching sorting order
      !
      integer(ik) :: i
      !
      cnt = 0
      scan_orbitals: do i=1,size(spin_all)
        if (spin_all(i)<0) cycle scan_orbitals
        cnt       = cnt + 1
        spin(cnt) = spin_all(i)
        ord (cnt) = ord_all (i)
      end do scan_orbitals
    end subroutine copy_plus
    !
    subroutine copy_minus(spin_all,ord_all,cnt,spin,ord)
      integer(ik), intent(in)  :: spin_all(:) ! Alpha (positive) and beta (negative) spin-orbitals, sorted
      integer(ik), intent(in)  :: ord_all (:) ! Sorting order, matching spin_all
      integer(ik), intent(out) :: cnt         ! Number of negative spin-orbitals
      integer(ik), intent(out) :: spin(:)     ! Negative spin-orbitals, with sign flipped, sorted
      integer(ik), intent(out) :: ord (:)     ! Matching sorting order
      !
      integer(ik) :: i
      !
      cnt = 0
      scan_orbitals: do i=size(spin_all),1,-1
        if (spin_all(i)>0) cycle scan_orbitals
        cnt       = cnt + 1
        spin(cnt) =-spin_all(i)
        ord (cnt) = ord_all (i)
      end do scan_orbitals
    end subroutine copy_minus
  end subroutine symmetry_sort_spinorbitals
  !
  !  WARNING: analyze_spinorbit_integral_blocks is:
  !           a) rather non-trivial
  !           b) on the critical path for many calculations
  !  It absolutely must be efficient - so everything what could be usefully
  !  precomputed and cached, is precomputed and cached. Effectively, half
  !  of this routine is found above, in symmetry_sort_spinorbitals()
  !
  subroutine analyze_spinorbit_integral_blocks(detbra,detket,bms,bmo)
    integer(ik), intent(in)           :: detbra    ! Parent determinant index
    integer(ik), intent(in)           :: detket    ! Ion determinant index
    type(bm_structure), intent(in)    :: bms       ! Block descriptor, input matrix
    type(bm_structure), intent(inout) :: bmo       ! Block descriptor, output matrix
    !
    integer(ik) :: na_bra, nb_bra     ! Number of alpha- and beta- orbitals in bra and ket
    integer(ik) :: na_ket, nb_ket     ! 
    integer(ik) :: row_a(bms%n_rows)  ! List of rows and columns in an alpha-spin block
    integer(ik) :: col_a(bms%n_cols)
    integer(ik) :: row_b(bms%n_rows)  ! List of rows and columns in an beta-spin block
    integer(ik) :: col_b(bms%n_cols)
    integer(ik) :: n_row_a, n_col_a, n_row_b, n_col_b
    integer(ik) :: ib
    integer(ik) :: r_pos, r_end       ! First and last posion in bms%rows
    integer(ik) :: c_pos, c_end       ! First and last posion in bms%cols
    !
    ! This routine is so much on the critical path, timing buggers it up completely
    ! if (detail_timers) call TimerStart('analyze_spinorbit_integral_blocks')
    !
    !  We could potentially double the number of blocks present in pure-space matrix
    !
    if (bmo%n_blocks/=-1) stop 'sd_core%analyze_spinorbit_integral_blocks - descriptor already in use?!'
    !
    bmo%n_blocks = 0
    bmo%n_rows   = 0
    bmo%n_cols   = 0
    bmo%max_row  = nelbra
    bmo%max_col  = nelket
    na_bra       = scnt_bra(1,detbra)
    nb_bra       = scnt_bra(2,detbra)
    na_ket       = scnt_ket(1,detket)
    nb_ket       = scnt_ket(2,detket)
    scan_blocks: do ib=1,bms%n_blocks
      r_pos = bms%block_data(3,ib)
      c_pos = bms%block_data(4,ib)
      r_end = bms%block_data(1,ib) + r_pos - 1
      c_end = bms%block_data(2,ib) + c_pos - 1
      !
      !  Find which rows and columns from this block have been taken, 
      !  separately for alpha and beta spins. This is basically an intersection of 
      !  row/column indices for this block, and corresponding spinbra/spinket
      !
      call bm_intersection(spinbra(:na_bra,1,detbra),ord_bra(:na_bra,1,detbra),bms%rows(r_pos:r_end),n_row_a,row_a)
      call bm_intersection(spinbra(:nb_bra,2,detbra),ord_bra(:nb_bra,2,detbra),bms%rows(r_pos:r_end),n_row_b,row_b)
      call bm_intersection(spinket(:na_ket,1,detket),ord_ket(:na_ket,1,detket),bms%cols(c_pos:c_end),n_col_a,col_a)
      call bm_intersection(spinket(:nb_ket,2,detket),ord_ket(:nb_ket,2,detket),bms%cols(c_pos:c_end),n_col_b,col_b)
      !
      !  See whether we actually have any non-empty blocks here
      !  stuff_block will alson take care of adjusting n_rows/n_cols counters
      !
      call bm_add_block(bmo,row_a(:n_row_a),col_a(:n_col_a))
      call bm_add_block(bmo,row_b(:n_row_b),col_b(:n_col_b))
    end do scan_blocks
    !
    ! if (detail_timers) call TimerStop('analyze_spinorbit_integral_blocks')
    !
  end subroutine analyze_spinorbit_integral_blocks
  !
  function hopeless_blocks(drop_row,rep_row,drop_col,rep_col,sms) result (hopeless)
    integer(ik), intent(in)        :: drop_row, rep_row ! Number of rows which may be dropped or replaced
    integer(ik), intent(in)        :: drop_col, rep_col ! Number of columns which may be dropped or replaced
    type(bm_structure), intent(in) :: sms               ! Overlap matrix structure
    logical                        :: hopeless
    !
    integer(ik) :: ib, nbr, nbc
    integer(ik) :: defect, defect_row, defect_col
    integer(ik) :: rep_any
    !
    !  Replacing a row may fix up a column, and vice versa
    !
    rep_any = rep_row + rep_col
    !
    !  If matrix has zero rows, we need to be able to drop or replace enough rows to compensate
    !
    defect_row = sms%max_row-sms%n_rows
    if (drop_row+rep_any<defect_row) then
      ! write (out,"('defect_row= ',i0,' drop_row= ',i0,' rep_any= ',i0)") defect_row, drop_row, rep_any
      hopeless = .true. ; return
    end if
    !
    !  ditto for columns
    !
    defect_col = sms%max_col-sms%n_cols
    if (drop_col+rep_any<defect_col) then
      ! write (out,"('defect_col= ',i0,' drop_col= ',i0,' rep_any= ',i0)") defect_col, drop_col, rep_any
      hopeless = .true. ; return
    end if
    !
    !  If matrix is non-square, we must be able to drop enough rows/columns to make it square
    !
    defect = sms%max_row - sms%max_col
    if (defect>drop_row .or. -defect>drop_col) then
      ! write (out,"('defect= ',i0,' drop_row= ',i0,' drop_col= ',i0)") defect, drop_row, drop_col
      hopeless = .true. ; return
    end if
    !
    !  Non-rectangular blocks will also yield zeros; however we now must differentiate
    !  between dropped rows/columns (which can only fix one sub-block) and replaced 
    !  rows/columns (which may fix the whole matrix, depending on where they land.
    !
    if (rep_any>0) then
      hopeless = .false. ; return
    end if
    !
    defect_row = 0
    defect_col = 0
    scan_blocks: do ib=1,sms%n_blocks
      nbr = sms%block_data(1,ib)
      nbc = sms%block_data(2,ib)
      defect_row = defect_row + max(nbr-nbc,0)
      defect_col = defect_col + max(nbc-nbr,0)
    end do scan_blocks
    hopeless = (defect_row>drop_row) .or. (defect_col>drop_col)
    ! if (hopeless) then
    !   write (out,"('defect_row= ',i0,' defect_col= ',i0,' drop_row= ',i0,' drop_col= ',i0)") defect_row,defect_col,drop_row,drop_col
    ! end if
  end function hopeless_blocks
  !
  !  Calculate a overlap integral of two determinants
  !
  function determinant_overlap(sdet,sdet_ms) result(overdet)
    real(rk), intent(in)           :: sdet(:,:)
    type(bm_structure), intent(in) :: sdet_ms    ! Structure of the overlap matrix
    real(rk)                       :: overdet
    !
    if (detail_timers) call TimerStart('determinant_overlap')
    if (use_symmetry) then
      overdet = block_det(sdet,sdet_ms)
    else 
      overdet = linpack_determinant(sdet)
    end if
    if (detail_timers) call TimerStop('determinant_overlap')
  end function determinant_overlap
  !
  !  Calculate 1-electron matrix element between two determinants
  !
  function determinant_transition(smo,vmo,smo_ms) result(vdet)
    real(rk), intent(in)           :: smo(:,:) ! Overlap matrix elements in the MO basis
    real(rk), intent(in)           :: vmo(:,:) ! 1-electron matrix elements in the MO basis
    type(bm_structure), intent(in) :: smo_ms   ! Structure of the overlap matrix
    real(rk)                       :: vdet     ! n-electron matrix element
    !
    real(rk)           :: scr(nelket,nelket)   ! Buffer
    integer(ik)        :: iscan 
    !
    if (detail_timers) call TimerStart('determinant_transition')
    vdet = 0
    if (use_symmetry) then
      !
      !  We have two important special cases, which are advantageous to handle separately
      !
      select case (nelket-smo_ms%n_rows)
        case default
        case (0)
          scr = smo
          accumulate_sym: do iscan=1,nelket
            vdet = vdet + substitute_determinant(iscan)
          end do accumulate_sym
        case (1) 
          iscan = bm_find_zero_row(smo_ms)
          if (iscan<=0) stop 'sd_core%determinant_transition - no zero row in a defective matrix?!'
          scr = smo
          vdet = substitute_determinant(iscan)
      end select
    else
      ! Original brute-force version
      scr = smo
      accumulate_nosym: do iscan=1,nelket
        scr(iscan,:) = vmo(iscan,:)
        vdet = vdet + linpack_determinant(scr)
        scr(iscan,:) = smo(iscan,:)
      end do accumulate_nosym
    end if
    if (detail_timers) call TimerStop('determinant_transition')
    !
    contains 
    function substitute_determinant(irow) result(det)
      integer(ik), intent(in) :: irow ! Row to substitute
      real(rk)                :: det
      !
      integer(ik)        :: ic
      integer(ik)        :: cols(nelket)  ! Connected columns
      integer(ik)        :: ncols         ! Number of connected columns
      type(bm_structure) :: scr_ms        ! Matrix structure for smat with one row substitution
      ! real(rk)           :: det_chk
      !
      ncols = 0
      count_columns: do ic=1,nelket
        if (abs(vmo(irow,ic))<eps_integral) cycle count_columns
        ncols = ncols + 1
        cols(ncols) = ic
      end do count_columns
      !
      if (ncols<=0) then  ! This is a zero row, no point in continuing
        det = 0
        return
      end if
      !
      call bm_substitute_row(smo_ms,irow,cols(:ncols),scr_ms)
      scr(iscan,:) = vmo(iscan,:)
      det = block_det(scr,scr_ms)
      ! det_chk = linpack_determinant(scr)
      ! write (out,"('Block = ',g24.16,' brute = ',g24.16,' diff = ',g24.16)") det, det_chk, det-det_chk
      scr(iscan,:) = smo(iscan,:)
      call bm_destroy(scr_ms)
    end function substitute_determinant
  end function determinant_transition
  !
  !  Drop row of a rectangular matrix
  !
  subroutine drop_row(row,src,dst)
    integer(ik), intent(in) :: row      ! Row to be dropped
    real(rk), intent(in)    :: src(:,:) ! Source matrix
    real(rk), intent(out)   :: dst(:,:) ! Destination matrix
    !
    dst(  1:row-1,:) = src(    1:row-1,:)
    dst(row:,     :) = src(row+1:,     :)
  end subroutine drop_row
  !
  !  Get a minor of the specified element of a square matrix
  !
  subroutine get_minor(row,col,src,dst)
    integer(ik), intent(in) :: row, col ! Element
    real(rk), intent(in)    :: src(:,:) ! Source matrix
    real(rk), intent(out)   :: dst(:,:) ! Destination matrix
    !
    if (size(src,dim=1)/=size(src,dim=2)) then
      stop 'sd_core%get_minor - src is not square'
    end if
    if (size(dst,dim=1)/=(size(src,dim=1)-1) .or. &
        size(dst,dim=2)/=(size(src,dim=2)-1) ) then
      stop 'sd_core%get_minor - src and dst are not compatible'
    end if
    !
    dst(  1:row-1,  1:col-1) = src(    1:row-1,    1:col-1)
    dst(  1:row-1,col:     ) = src(    1:row-1,col+1:     )
    dst(row:     ,  1:col-1) = src(row+1:     ,    1:col-1)
    dst(row:     ,col:     ) = src(row+1:     ,col+1:     )
  end subroutine get_minor
  !
  !  Calculate a determinant of a square matrix
  !
  real(rk) function determinant(mat)
    real(rk), intent(inout) :: mat(:,:) ! Matrix destroyed on exit
    !
    ! GAMESS built-in routine
    ! call poldet(determinant,mat,size(mat,dim=1),size(mat,dim=1))
    !
    ! LINPACK routines
    if (verbose>=4) then
      write (out,"(/'Requested determinant for matrix: ')")
      call print_matrix(mat)
    end if
    !
    determinant = linpack_determinant_trash_input(mat)
  end function determinant
  !
  !  Consistency test - normalization and S^2 for the total wavefunction.
  !  We use ladder operator representation of S^2:
  !
  !    S^2 = 0.5 * S_+ S_- + 0.5 * S_ S_+ + S_z^2 
  !
  !  This code is not terribly efficient, but this should not really matter (?)
  !
  subroutine calculate_s2(c,d,ov,s2)
    real(rk), intent(in)    :: c(  :)  ! Amplitudes of each determinant in the w.f.
    integer(ik), intent(in) :: d(:,:)  ! Orbital list for each determinant
    real(rk), intent(out)   :: ov      ! Overlap integral
    real(rk), intent(out)   :: s2      ! S^2 expectation value
    !
    integer(ik) :: ndet, id, is, norb, io, jo
    integer(ik) :: td (size(d,dim=1))   ! Space for a temporary determinant
    integer(ik) :: par(size(d,dim=2))   ! "textbook" parity factor for each determinant
    !
    call TimerStart('S^2 evaluation')
    norb = size(d,dim=1)
    ndet = size(d,dim=2)
    ov   = 0._rk
    s2   = 0._rk
    !$omp parallel default(none) &
    !$omp&         shared(c,d,ndet,norb,par) private(id,is,td,io,jo) reduction(+:ov,s2)
    !$omp do
    get_parity: do id=1,ndet
      par(id) = det_parity(d(:,id))
    end do get_parity
    !$omp end do
    ! print_parity: do id=1,ndet
    !   write (out,*) 'id = ',id,' par = ',par(id)
    ! end do print_parity
    ov   = 0._rk
    s2   = 0._rk
    !$omp do
    scan_determinants: do id=1,ndet
      ov  = ov + c(id)**2
      orbital_i: do io=1,norb
        if (abs(d(io,id))/=1) cycle orbital_i
        !
        !  Only partially open subshells could contribute to s^2; closed 
        !  subshells give identical zero upon action on the total spin operator
        !
        !  Diagonal part of S+ S- + S- S+ ; determinant does not change
        !
        s2 = s2 + 0.5_rk * c(id)**2
        !
        orbital_j: do jo=1,norb
          if (abs(d(jo,id))/=1) cycle orbital_j
          !
          !  S_z^2; determinant does not change
          !
          s2 = s2 + 0.25_rk * d(io,id) * d(jo,id) * c(id)**2 
          !
          !  Off-diagonal part of S+ S- + S- S+; determinant changes
          !
          if (d(io,id)/=d(jo,id)) then
            td = d(:,id) ; td(io) = -td(io) ; td(jo) = -td(jo)
            find_match: do is=1,ndet
              if (any(d(:,is)/=td)) cycle find_match
              !  S_+ S_- have interchanged the location of alpha and beta spins; to bring 
              !  us back into the canonical alpha/beta order, we have to swap the spin-
              !  and spatial orbitals into the expected order.
              !
              !  The simplest way of accomplishing this operation is to multiply the
              !  alpha/beta to "textbook" parities for the initial and final states.
              !
              s2 = s2 + 0.5_rk * c(id) * c(is) * par(id) * par(is)
            end do find_match
          end if
        end do orbital_j
      end do orbital_i
    end do scan_determinants
    !$omp end do
    !$omp end parallel
    call TimerStop('S^2 evaluation')
  end subroutine calculate_s2
  !
  !  Calculate parity factor connecting amplitudes of the same determinant,
  !  with orbitals in the "textbook" (alpha and beta factors intermixed, spatial
  !  orbitals in the increasing order) and the alpha/beta (all alpha before all
  !  beta) orders.
  !
  !  Since this routine is not on a critical path, we'll be very stupid: we'll
  !  assign an "alpha/beta" index to each "textbook" orbital, then buble-sort
  !  the resulting array while counting the number of interchanges.
  !
  integer(ik) function det_parity(det)
    integer(ik), intent(in) :: det(:) ! Determinant, what else?
    !
    integer(ik) :: ord(2*size(det)) ! Ordering array - can't get any bigger than this
    integer(ik) :: io               ! Orbital number
    integer(ik) :: it               ! Current orbital number, "textbook" order
    integer(ik) :: ia, ib           ! Current alpha/beta orbital number, "alpha/beta" order
    integer(ik) :: nel, nalp, nbet  ! Number of electrons - total, alpha, and beta
    logical     :: sorted
    !
    nel  = sum(abs(det))
    nalp = count(det==1) + count(det==2)
    nbet = nel - nalp
    !
    !  For each "textbook" orbital, figure out the desired "alpha/beta" index
    !
    ia = 0 ; ib = nalp ; it = 0
    order_initialize: do io=1,size(det)
      if (any(det(io)==(/ 1,2/))) then
        it = it + 1 ; ia = ia + 1 ; ord(it) = ia
      end if
      if (any(det(io)==(/-1,2/))) then
        it = it + 1 ; ib = ib + 1 ; ord(it) = ib
      end if
    end do order_initialize
    if (ia/=nalp .or. ib/=nel .or. it/=nel) then
      write (out,"('     Got: ia = ',i4,' ib = ',i4,' it = ',i4)") ia, ib, it
      write (out,"('Expected: ia = ',i4,' ib = ',i4,' it = ',i4)") nalp, nel, nel
      stop 'sd_core%det_parity - count error'
    end if
    !
    !  Count ordering permutations
    !
    det_parity = 1
    sorted = .false.
    ordering_sweep: do while(.not.sorted)
      sorted = .true.
      scan_order: do it=2,nel
        if (ord(it-1)<ord(it)) cycle scan_order
        ia        = ord(it-1)
        ord(it-1) = ord(it)
        ord(it)   = ia
        sorted    = .false. 
        det_parity= -det_parity
      end do scan_order
    end do ordering_sweep
  end function det_parity
  !
  !  Punch one vector in restart format
  !
  subroutine punch_vector(imo,vec)
    integer(ik), intent(in) :: imo    ! Identifying index
    real(rk), intent(in)    :: vec(:) ! Vector to punch
    !
    integer(ik) :: ibas, tbas, line
    real(rk)    :: buf(5)
    !
    line = 1
    punch_line: do ibas=1,size(vec),5
      tbas = min(size(vec),ibas+4)
      buf = 0._rk
      buf(1:tbas-ibas+1) = vec(ibas:tbas)
      write (out,"(i2,i3,5e15.8)") mod(imo,100), line, buf
      line = line + 1
    end do punch_line
  end subroutine punch_vector
  !
  subroutine print_matrix(m)
    real(rk), intent(in) :: m(:,:) ! Matrix to be printed
    !
    integer(ik), parameter :: cb = 16 ! Batch size
    integer(ik)            :: nr, nc  ! Number of rows and columns
    integer(ik)            :: r,  c   ! Row and column
    integer(ik)            :: c1, c2  ! Current batch of columns
    !
    nr = size(m,dim=1)
    nc = size(m,dim=2)
    col_batch: do c1=1,nc,cb
      c2 = min(nc,c1+cb-1)
      write (out,"()")
      write (out,"(1x,5x,16(1x,i10))") (c,c=c1,c2)
      !
      rows: do r=1,nr
        write (out,"(1x,i5,16(1x,f10.6))") r, (m(r,c),c=c1,c2)
      end do rows
      !
      write (out,"()")
    end do col_batch
  end subroutine print_matrix
  !
  subroutine report_progress(detref)
    !$ use OMP_LIB
    integer(ik), intent(in) :: detref ! Current progress indicator
    logical, save           :: first = .true.
    integer(ik), save       :: p_interval, p_max, p_last
    integer(ik)             :: p_now
    !
    !$ if (omp_get_thread_num()/=0) return
    if (first) then
      first = .false.
      call system_clock(p_last,p_interval,p_max)
      p_interval = 10*p_interval  ! 10 seconds
    else
      !
      !  Progress indicator - try not to be too verbose
      !
      call system_clock(p_now)
      if ( (p_last+p_interval< p_max .and. p_now>=p_last+p_interval) .or. &
           (p_last+p_interval>=p_max .and. p_now< p_last) ) then
        write (out,"('Processing reference determinant ',i6,' of ',i6)") detref, ndetbra
        p_last = p_now
        call flush(out)
      end if
    end if
  end subroutine report_progress
end module sd_core
