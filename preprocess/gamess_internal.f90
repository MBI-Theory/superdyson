!
!  2016 Sept 14 - Extended to support H and I orbitals (L=5,6)
!
!  Internal definitions needed for working directly with GAMESS orbitals
!  and data files. This module is tightly coupled to import_gamess and
!  ecp_gamess. It should not be directly relied upon by anything else.
!
!  We will keep everything here is extended precision; routines requiring
!  lower precision must typecast all constants at the point of use.
!  The only exception are the molecular orbitals - these are kept in
!  "high" precision - our input files do not provide more in any event
!
!  In order to maintain computational costs under control for routines
!  which do not need high precision, we'll also keep a copy of the
!  frequently-used terms in standard precision.
!
 module gamess_internal
   use accuracy
   use constants
   implicit none
   public
   private a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q
   !
   character(len=clen), save :: rcsid_gamess_internal = "$Id: gamess_internal.f90,v 1.4 2021/09/29 13:43:22 ps Exp $"
   !
    integer(ik), parameter :: gam_max_oversample     = 5     ! Maximum oversampling order (per Cartesian direction)
    integer(ik), parameter :: gam_file               = 57    ! A semi-randomly chosen unit
    integer(ik), parameter :: gam_max_line           = 192
    integer(ik), parameter :: gam_max_shells         = 600
    integer(ik), parameter :: gam_max_ecp_terms      = 100
    integer(ik), parameter :: gam_max_primitive      = 2000
    integer(ik), parameter :: gam_max_atoms          = 1000
    integer(ik), parameter :: gam_max_contraction    = 50
    integer(ik), parameter :: gam_max_batch_size     = 60 ! Split up very long batches; otherwise, memory requirements 
                                                          ! become -very- large
    !
    integer(ik), parameter :: gam_lmax               = 6  ! Highest angular momentum
    integer(ik), parameter :: gam_nxyz               = 84 ! Number of possible XYZ products
    ! orbcnt(l) = (l+1)*(l+2)/2
    integer(ik), parameter :: gam_orbcnt(0:gam_lmax) = (/ 1, 3, 6, 10, 15, 21, 28 /)
    ! normc(l)  = Gamma[l+1/2] / sqrt(pi) = (2l-1)!!/2**l
    real(xrk), parameter   :: gam_normc (0:gam_lmax) = (/ 1._xrk, 1._xrk/2._xrk, 3._xrk/4._xrk, 15._xrk/8._xrk, 105._xrk/16._xrk, &
                                                          945._xrk/32._xrk, 10395._xrk/64._xrk /)
    !
    !  Definition of the angular parts of GAMESS basis functions
    ! 
    real(xrk), parameter   :: a = 1.73205080756887729352744634150587236694281_xrk ! sqrt( 3.0)
    real(xrk), parameter   :: b = 2.23606797749978969640917366873127623544062_xrk ! sqrt( 5.0)
    real(xrk), parameter   :: c = 3.87298334620741688517926539978239961083292_xrk ! sqrt(15.0)
    real(xrk), parameter   :: d = 2.64575131106459059050161575363926042571026_xrk ! sqrt( 7.0)
    real(xrk), parameter   :: e = 3.41565025531986612774034622684035066357834_xrk ! sqrt(35.0/3.0)
    real(xrk), parameter   :: f = 5.91607978309961604256732829156161704841550_xrk ! sqrt(35.0)    
    real(xrk), parameter   :: g = 3.00000000000000000000000000000000000000000_xrk !
    real(xrk), parameter   :: h = 4.58257569495584000658804719372800848898446_xrk !
    real(xrk), parameter   :: i = 7.93725393319377177150484726091778127713078_xrk !
    real(xrk), parameter   :: j = 10.2469507659595983832210386805210519907350_xrk !
    real(xrk), parameter   :: k = 3.31662479035539984911493273667068668392709_xrk !
    real(xrk), parameter   :: l = 5.74456264653802865985061146821892931822026_xrk !
    real(xrk), parameter   :: m = 9.94987437106619954734479821001206005178127_xrk !
    real(xrk), parameter   :: n = 6.79705818718657137846900052738478763345227_xrk !
    real(xrk), parameter   :: o = 1.00000000000000000000000000000000000000000_xrk !
    real(xrk), parameter   :: p = 15.1986841535706636316687442060396450375823_xrk !
    real(xrk), parameter   :: q = 19.6214168703485834685260037891817743263128_xrk !
    !
    ! ang_loc(l) = sum(gam_orbcnt(0:l-1))
    integer(ik), parameter :: ang_loc(0:gam_lmax+1) = (/ 0, 1, 4, 10, 20, 35, 56, 84 /)
    ! ang_n? defines function ordering. Ad hoc, so be careful with the table!
                                                     !                      1                   2                   3
                                                     !  0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4
                                                     !                                          X Y Z X X Y Y Z Z X X Y X Y Z
                                                     !                      X Y Z X X Y Y Z Z X X Y Z X X Y Y Z Z X X Y X Y Z
                                                     !          X Y Z X X Y X Y Z X X Y Y Z Z Y X Y Z X X Y Y Z Z Y Z Z Y X X
                                                     !  S X Y Z X Y Z Y Z Z X Y Z Y Z X Z X Y Z X Y Z Y Z X Z X Y Y Z Z Z Z Y
    integer(ik), parameter :: ang_nx(0:gam_nxyz-1) = (/ 0,1,0,0,2,0,0,1,1,0,3,0,0,2,2,1,0,1,0,1,4,0,0,3,3,1,0,1,0,2,2,0,2,1,1, &
                                                        5,0,0,4,4,1,0,1,0,3,3,2,0,2,0,3,1,1,2,2,1, &
                                                        6,0,0,5,5,1,0,1,0,4,4,2,0,2,0,4,1,1,3,3,0,3,3,2,1,2,1,2/)
    integer(ik), parameter :: ang_ny(0:gam_nxyz-1) = (/ 0,0,1,0,0,2,0,1,0,1,0,3,0,1,0,2,2,0,1,1,0,4,0,1,0,3,3,0,1,2,0,2,1,2,1, &
                                                        0,5,0,1,0,4,4,0,1,2,0,3,3,0,2,1,3,1,2,1,2, &
                                                        0,6,0,1,0,5,5,0,1,2,0,4,4,0,2,1,4,1,3,0,3,2,1,3,3,1,2,2 /)
    integer(ik), parameter :: ang_nz(0:gam_nxyz-1) = (/ 0,0,0,1,0,0,2,0,1,1,0,0,3,0,1,0,1,2,2,1,0,0,4,0,1,0,1,3,3,0,2,2,1,1,2, &
                                                        0,0,5,0,1,0,1,4,4,0,2,0,2,3,3,1,1,3,1,2,2, &
                                                        0,0,6,0,1,0,1,5,5,0,2,0,2,4,4,1,1,4,0,3,3,1,2,1,2,3,3,2 /)
    real(xrk), parameter   :: ang_c (0:gam_nxyz-1) = (/ o,o,o,o,o,o,o,a,a,a,o,o,o,b,b,b,b,b,b,c,o,o,o,d,d,d,d,d,d,e,e,e,f,f,f, &
                                                        o,o,o,g,g,g,g,g,g,h,h,h,h,h,h,i,i,i,j,j,j, &
                                                        o,o,o,k,k,k,k,k,k,l,l,l,l,l,l,m,m,m,n,n,n,p,p,p,p,p,p,q /)
    real(rk), parameter    :: ang_c_rk(0:gam_nxyz-1) = real(ang_c,kind=rk)
    integer(ik), parameter :: ang_nxyz(0:gam_nxyz-1,3) = reshape((/ang_nx,ang_ny,ang_nz/),(/gam_nxyz,3/))
    !
    character(len=1), parameter :: shell_name(0:gam_lmax) = (/ 'S', 'P', 'D', 'F', 'G', 'H', 'I' /)
    character(len=6), parameter :: ang_name(0:gam_nxyz-1) = (/ &
        ' S    ', &
        ' X    ', ' Y    ', ' Z    ', &
        'XX    ', 'YY    ', 'ZZ    ', 'XY    ', 'XZ    ', 'YZ    ', &
        'XXX   ', 'YYY   ', 'ZZZ   ', 'XXY   ', 'XXZ   ', 'YYX   ', 'YYZ   ', 'ZZX   ', 'ZZY   ', 'XYZ   ', &
        'XXXX  ', 'YYYY  ', 'ZZZZ  ', 'XXXY  ', 'XXXZ  ', 'YYYX  ', 'YYYZ  ', 'ZZZX  ', 'ZZZY  ', 'XXYY  ', &
        'XXZZ  ', 'YYZZ  ', 'XXYZ  ', 'YYXZ  ', 'ZZXY  ', &
        'X5    ', 'Y5    ', 'Z5    ', 'X4Y   ', 'X4Z   ', 'Y4X   ', 'Y4Z   ', 'Z4X   ', 'Z4Y   ', 'X3Y2  ', &
        'X3Z2  ', 'Y3X2  ', 'Y3Z2  ', 'Z3X2  ', 'Z3Y2  ', 'X3YZ  ', 'Y3XZ  ', 'Z3XY  ', 'X2Y2Z ', 'X2Z2Y ', &
        'Y2Z2X ',  &
        'X6    ', 'Y6    ', 'Z6    ', 'X5Y   ', 'X5Z   ', 'Y5X   ', 'Y5Z   ', 'Z5X   ', 'Z5Y   ', 'X4Y2  ', &
        'X4Z2  ', 'Y4X2  ', 'Y4Z2  ', 'Z4X2  ', 'Z4Y2  ', 'X4YZ  ', 'Y4XZ  ', 'Z4XY  ', 'X3Y3  ', 'X3Z3  ', &
        'Y3Z3  ', 'X3Y2Z ', 'X3Z2Y ', 'Y3X2Z ', 'Y3Z2X ', 'Z3X2Y ', 'Z3Y2X ', 'X2Y2Z2' /)
    !
    !  The "drop_?" arrays make it easier to construct Obara-Saika recurrences 
    !  for the integrals. One could also compute drop tables on the fly from the
    !  ang_?? tables - however, having them as a parameter lets the compiler to
    !  unroll recurrence loops and convert them to optimal data flow graphs.
          !                                 1                             2                             3              
          !   0  1  2  3  4  5  6  7  8  9  0  1  2  3  4  5  6  7  8  9  0  1  2  3  4  5  6  7  8  9  0  1  2  3  4
          !                                                               X  Y  Z  X  X  Y  Y  Z  Z  X  X  Y  X  Y  Z
          !                                 X  Y  Z  X  X  Y  Y  Z  Z  X  X  Y  Z  X  X  Y  Y  Z  Z  X  X  Y  X  Y  Z
          !               X  Y  Z  X  X  Y  X  Y  Z  X  X  Y  Y  Z  Z  Y  X  Y  Z  X  X  Y  Y  Z  Z  Y  Z  Z  Y  X  X
          !   S  X  Y  Z  X  Y  Z  Y  Z  Z  X  Y  Z  Y  Z  X  Z  X  Y  Z  X  Y  Z  Y  Z  X  Z  X  Y  Y  Z  Z  Z  Z  Y
    integer(ik), parameter :: drop_x(0:gam_nxyz-1) = &
          (/ -1, 0,-1,-1, 1,-1,-1, 2, 3,-1, 4,-1,-1, 7, 8, 5,-1, 6,-1, 9,10,-1,-1,13,14,11,-1,12,-1,15,17,-1,19,16,18, &
             20,-1,-1,23,24,21,-1,22,-1,29,30,25,-1,27,-1,32,26,28,33,34,31,35,-1,-1,38,39,36,-1,37,-1,44,45,40,-1,42, &
             -1,50,41,43,46,48,-1,53,54,51,47,52,49,55 /)
    integer(ik), parameter :: drop_y(0:gam_nxyz-1) = &
          (/ -1,-1, 0,-1,-1, 2,-1, 1,-1, 3,-1, 5,-1, 4,-1, 7, 9,-1, 6, 8,-1,11,-1,10,-1,15,16,-1,12,13,-1,18,14,19,17, &
             -1,21,-1,20,-1,25,26,-1,22,23,-1,29,31,-1,28,24,33,27,32,30,34,-1,36,-1,35,-1,40,41,-1,37,38,-1,46,47,-1, &
             43,39,51,42,44,-1,49,50,45,53,55,48,52,54 /)
    integer(ik), parameter :: drop_z(0:gam_nxyz-1) = &
          (/ -1,-1,-1, 0,-1,-1, 3,-1, 1, 2,-1,-1, 6,-1, 4,-1, 5, 8, 9, 7,-1,-1,12,-1,10,-1,11,17,18,-1,14,16,13,15,19, &
             -1,-1,22,-1,20,-1,21,27,28,-1,24,-1,26,30,31,23,25,34,29,32,33,-1,-1,37,-1,35,-1,36,42,43,-1,39,-1,41,48, &
             49,38,40,52,-1,45,47,44,50,46,51,54,55,53 /)
    integer(ik), parameter :: drop_xyz(0:gam_nxyz-1,3) = reshape((/drop_x,drop_y,drop_z/),(/gam_nxyz,3/))
    !
    type gam_atom
      character(len=20)    :: name                         ! Name of the atom
      real(xrk)            :: znuc                         ! Effective nuclear charge. In the case of an ECP, this charge
                                                           ! gives the asymptotic long-range potential of this nucleus
      real(xrk)            :: xyz(3)                       ! Coordinates, in Angstrom
      real(rk)             :: xyz_rk(3)                    ! ditto, standard real
      integer(ik)          :: nshell                       ! Number of contracted shells
      integer(ik)          :: sh_l    (gam_max_shells)     ! Per-shell angular momentum
      integer(ik)          :: sh_p    (gam_max_shells+1)   ! Index of the first primitive
      real(xrk)            :: p_zet   (gam_max_primitive)  ! Primitive exponents
      real(rk)             :: p_zet_rk(gam_max_primitive)  ! ditto, standard real
      real(xrk)            :: p_c     (gam_max_primitive)  ! Contraction coefficients
      real(rk)             :: p_c_rk  (gam_max_primitive)  ! ditto, standard real
      !
      character(len=20)    :: ecp_name                     ! Name of the ECP
      real(xrk)            :: ecp_zcore                    ! Number of electrons in the ECP
      integer(ik)          :: ecp_lmaxp1                   ! LMAX+1 for this ECP
      integer(ik)          :: ecp_nterms                   ! Number of distinct Gaussians in the ECP
      integer(ik)          :: ecp_l(gam_max_ecp_terms)     ! Angular momentum of each projector. Negative
                                                           ! value means a universal contribution, applicable
                                                           ! to all channels.
      integer(ik)          :: ecp_n(gam_max_ecp_terms)     ! Power of r in an ECP term, plus 2
      real(xrk)            :: ecp_c(gam_max_ecp_terms)     ! Weight of the Gaussian term
      real(xrk)            :: ecp_d(gam_max_ecp_terms)     ! Exponent of the Gaussian term
    end type gam_atom
    !
    type gam_structure
      integer(ik)                    :: natoms   = 0          ! Number of atoms in the structure
      integer(ik)                    :: nbasis   = 0          ! Number of Cartesian basis functions
      integer(ik)                    :: nvectors = 0          ! Number of vectors in this set
      type(gam_atom), allocatable    :: atoms(:)              ! Table of atoms
      character(len=13), allocatable :: bas_labels(:)         ! Basis labels
      real(rk), allocatable          :: vectors(:,:)          ! Molecular orbitals
      real(xrk)                      :: rotmat(3,3)           ! Rotation matrix for this structure
      real(rk)                       :: rotmat_rk(3,3)        ! ditto, standard precision
      real(rk)                       :: dx(3)                 ! Grid spacing for oversampling.
                                                              ! All-negative input disables oversampling.
      !
      !  The following two fields are used only for evaluation of 2-electron integrals
      !
      integer(ik)                    :: nbatches   = 0        ! Number of batches in this structure
                                                              ! Consecuitive shells with the same L will be grouped
                                                              ! into the same batch - this dramatically improves 
                                                              ! utilization of primitive integrals for contracted basis sets.
      integer(ik)                    :: max_shells = 0        ! Number of shells in the largest batch
      integer(ik)                    :: max_batch  = 0        ! Number of orbitals in the largest batch
      integer(ik)                    :: max_contractions = 0  ! Longest contraction across the entire structure
      integer(ik), allocatable       :: batches(:,:)          ! First index:
                                                              ! 1   = atom of the batch, index into atoms()
                                                              ! 2   = number of contracted shells in this batch
                                                              ! 3   = index of the first orbital of the batch in the list of
                                                              !       basis functions
                                                              ! 4   = number of orbitals in this batch
                                                              ! 4+i = shell of the batch, index into atoms%sh_l(:), atoms%sh_p(:)
                                                              !       i = 1 ... [2]
    end type gam_structure
    !
    !  Data structure for passing around extra information for 3-centre integrals.
    !  We can't make this data static: this would prevent parallel execution of the
    !  integral package.
    !
    type gam_operator_data
      character(len=20)              :: op_name               ! Name of the operator
      real(xrk)                      :: op_xyz(3)             ! Operator-carrying centre, in Bohr
      real(xrk)                      :: omega                 ! exponent of the exp(-omega*r**2) factor
      real(xrk)                      :: imag_rc               ! Imaginary part of op_xyz (see os_basic_integral)
      real(xrk)                      :: kvec(3)               ! Used by planewave operator
      integer(ik)                    :: op_i, op_j            ! Additional indices, used by 'r-d/dr' operator
    end type gam_operator_data
    !
    character(len=GAM_MAX_LINE)      :: gam_line_buf          ! Line buffer for reading GAMESS files
    !
    contains
    !
    !  The normalization of the primitives is for the function:
    !
    !    z**l exp(-zeta*r**2)
    !
    real(xrk) function primitive_norm(l,z)
      integer(ik), intent(in) :: l ! Angular momentum
      real(xrk), intent(in)   :: z ! Orbital exponent
      !
      primitive_norm =  sqrt( (2*z)**(1.5_xrk+l) / (gam_normc(l)*pi_xrk**1.5_xrk) ) ;
    end function primitive_norm
    !
    !  Expectation of r for a normalized primitive
    !
    real(xrk) function primitive_r_expectation(l,z)
      integer(ik), intent(in) :: l ! Angular momentum
      real(xrk), intent(in)   :: z ! Orbital exponent
      !
      primitive_r_expectation = (1.5_xrk+l)/z
    end function primitive_r_expectation
    !
    !  Generate basis function names
    !
    subroutine label_basis_set(gam)
      type (gam_structure), intent(inout) :: gam
      !
      integer(ik) :: ip, iat, ish, ish_l, im
      !
      if (allocated(gam%bas_labels)) deallocate(gam%bas_labels)
      allocate (gam%bas_labels(gam%nbasis))
      ip = 1
      label_basis: do iat=1,gam%natoms
        label_shells: do ish=1,gam%atoms(iat)%nshell
          ish_l = gam%atoms(iat)%sh_l(ish)
          label_functions: do im=ang_loc(ish_l),ang_loc(ish_l+1)-1
            write (gam%bas_labels(ip),"(1x,a2,i3,1x,a6)") gam%atoms(iat)%name(1:2), iat, ang_name(im)
            ip = ip + 1
          end do label_functions
        end do label_shells
      end do label_basis
      if (ip-1/=gam%nbasis) stop 'gamess_internal%label_basis_set - count error'
    end subroutine label_basis_set
    !
    subroutine report_gamess_data(iu,gam,label,xlabel,xv)
      integer(ik), intent(in)                 :: iu     ! Unit to write the report to
      type(gam_structure), intent(in)         :: gam    ! Structure to report
      character(len=*), intent(in)            :: label  ! Comment to put in the $DATA section
      character(len=*), intent(in), optional  :: xlabel ! Section title for the extra data
      real(rk), intent(in), optional          :: xv(:)  ! Extra data to write out
      !
      integer(ik) :: iat, ivec
      !
      !  Data section
      !
      write (iu,"(' $DATA')")
      write (iu,"(a)") trim(label)
      write (iu,"('C1')")
      print_atoms: do iat=1,gam%natoms
        call print_one_atom(iu,gam%atoms(iat))
      end do print_atoms
      write (iu,"(' $END')")
      !
      !  Extra data section (if specified)
      !
      if (present(xlabel).and.present(xv)) then
        write (iu,"(1x,a)") trim(xlabel)
        write (iu,"(5(1x,g15.9))") xv
        write (iu,"(' $END')")
      end if
      !
      !  Orbitals section (if allocated)
      !
      if (allocated(gam%vectors)) then
        write (iu,"(' $VEC')") 
        print_vectors: do ivec=1,gam%nvectors
          call punch_vector(iu,ivec,gam%vectors(:,ivec))
        end do print_vectors
        write (iu,"(' $END')")
      end if
    end subroutine report_gamess_data
    !
    !  Print one atom, including the basis set
    !
    subroutine print_one_atom(iu,at)
      integer(ik), intent(in)    :: iu  ! Unit to use
      type(gam_atom), intent(in) :: at  ! Atom to report
      !
      integer(ik) :: is, sh_l, p1, pn, ip
      real(xrk)   :: norm
      !
      write (iu,"(1x,a11,1x,f10.6,1x,3(1x,f14.8))") trim(at%name), at%znuc, at%xyz
      print_shells: do is=1,at%nshell
        sh_l = at%sh_l(is)
        p1   = at%sh_p(is)
        pn   = at%sh_p(is+1)-1
        write (iu,"(a1,1x,i4)") shell_name(sh_l), pn-p1+1
        print_primitives: do ip=p1,pn
          norm = primitive_norm(sh_l,at%p_zet(ip))
          write (iu,"(i3,1x,g44.34,2x,g44.34)") ip-p1+1, at%p_zet(ip), at%p_c(ip)/norm
        end do print_primitives
      end do print_shells
      write (iu,"()")
    end subroutine print_one_atom
    !
    !  Punch one vector in GAMESS-US restart format
    !
    subroutine punch_vector(iu,imo,vec)
      integer(ik), intent(in) :: iu     ! I/O unit
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
        write (iu,"(i2,i3,5e15.8)") mod(imo,100), line, buf
        line = line + 1
      end do punch_line
    end subroutine punch_vector
    !
    !  Punch one vector in our spin-orbit restart format
    !  Since we are breaking GAMESS-US format anyways, we'll also
    !  store more significant digits.
    !
    subroutine punch_vector_complex(iu,imo,vec)
      integer(ik), intent(in) :: iu     ! I/O unit
      integer(ik), intent(in) :: imo    ! Identifying index
      complex(rk), intent(in) :: vec(:) ! Vector to punch
      !
      integer(ik) :: ibas, tbas, line
      complex(rk) :: buf(4)
      !
      line = 1
      punch_line: do ibas=1,size(vec),4
        tbas = min(size(vec),ibas+3)
        buf = 0._rk
        buf(1:tbas-ibas+1) = vec(ibas:tbas)
        write (iu,"(i2,i3,8e20.13)") mod(imo,100), line, buf
        line = line + 1
      end do punch_line
    end subroutine punch_vector_complex

    subroutine gam_readline(eof)
      logical, intent(out), optional :: eof
      integer(ik)                    :: ios
      !
      if (present(eof)) eof = .false.
      read (gam_file,"(a)",iostat=ios) gam_line_buf
      if (ios<0 .and. present(eof)) then
        eof = .true.
      else
        if (ios/=0) then
          write (out,"('gamess_internal%gam_readline: error ',i8,' reading .dat file')") ios
          stop 'gamess_internal%gam_readline - Read error'
        end if
      end if
    end subroutine gam_readline

 end module gamess_internal
