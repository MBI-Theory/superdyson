!
!  29 January 2006, Serguei Patchkovskii - Initial version
!  23 November 2006, SP - Extended to handle overlap integrals, transition
!                         dipoles and non-orthogonal exchange correction
!  05 December 2008, SP - Added OpenMP parallelization
!  14 March 2009, SP    - Separated out "core" routines, which can be
!                         used in other related applications
!
!  Calculation of Dyson orbitals for small CI expansions.
!
!  This code is based on the Dyson orbital routines, which were included in
!  GAMESS. Because the calculation we do here would not be easy to set up
!  inside GAMESS, I've decided to separate the code entirely. The final
!  result of the calculation, Dyson orbital expanded in AOs is still intended
!  to be added to the GAMESS punch file, for later processing and 
!  visualization.
!
!  The calculation is a naive evaluation of contributions arising from
!  overlaps between individual determinants contributing to two wavefunctions:
!  the parent molecule's wavefunction (with N electrons) and the target ion's
!  wavefunction (with N-1 electrons).
!
!  In many cases, using expansion based on natural orbitals of the target state
!  would give a very short sequence of determinants (ideally, just one on each
!  side).
!
!  There are five pieces of input needed to perform this calculations:
!
!  1. Orbitals of the parent molecule
!  2. Orbitals of the ion
!  3. Overlap and (optionally) dipole integrals in the AO basis
!  4. List of determinants for the parent
!  5. List of determinants for the ion
!
!  See read_core_input code in sd_core.f90 for the ad-hoc description of the input
!
module superdyson
  use sd_core
  implicit none
  private
  public call_superdyson
  public rcsid_superdyson
  !
  character(len=clen), save :: rcsid_superdyson = "$Id: superdyson.f90,v 1.30 2021/09/29 13:43:22 ps Exp $"
  !
  !  Local data for dyson
  !
  real(rk)                :: det_included     ! Number of determinant pairs included
  real(rk)                :: det_cutoff       ! Number of determinant pairs cut off
  real(rk)                :: det_distant      ! Number of too distant determinant pairs
  real(rk)                :: det_earlyzero    ! Number of determinants found zero by block structure
  !
  real(rk)                :: amplitude        ! 2-norm of the Dyson orbital
  !
  !  Our results. The last index is spin
  !
  real(rk)                :: overterm           ! Overlap between two given determinats (mode_overlap)
  real(rk)                :: overall            ! Running sum for the overlap (mode_overlap)
  real(rk), allocatable   :: dipall(:)          ! ... same for the transition dipole
  real(rk), allocatable   :: wdyson    (:,:,:)  ! Coefficients for the Dyson/exchange orbital in MO basis
  real(rk), allocatable   :: aodyson   (:,:,:)  ! Coefficients for the Dyson/exchnage orbital in AO basis
!
contains
!
  subroutine call_superdyson
    integer(ik)            :: detref, detion, ic
    integer(ik)            :: distance, tolerance
    integer(ik)            :: drop_row, rep_row    ! Number of rows we may be dropping and replacing
    integer(ik)            :: drop_col, rep_col    ! ditto for columns
    character(len=40)      :: tag
    real(rk), allocatable  :: sdet  (:,:)          ! 1-particle overlaps of spin-orbitals
    real(rk), allocatable  :: dipdet(:,:,:)        ! Spin-orbital dipole integrals
    real(rk), allocatable  :: wdysondet (:,:)      ! One determinant's contribution to the Dyson/exchange orbital
    real(rk), allocatable  :: wdyson_thread(:,:,:) ! Thread-private copy of wdyson
    real(rk), allocatable  :: dipterm(:)           ! Transition dipole ...
    real(rk), allocatable  :: dipall_thread(:)     ! Transition dipole ...
    type(bm_structure)     :: sdet_ms              ! Block-descriptor of sdet
    !
    call TimerStart('Superdyson')
    call TimerStart('Initialization')
    !
    !  Done reading input, begin calculations
    !
    call check_input_consistency
    !
    !
    !  In situations where left and right basis sets might differ, we'll 
    !  project the result on the ion orbitals.
    !
    allocate (wdyson(nmobra,2,1+op_components),aodyson(naos_bra,2,1+op_components),dipall(op_components))
    !
    overall = 0._rk
    dipall  = 0._rk
    !
    !  Transform integrals to the MO basis
    !
    call TimerStart('AO->MO')
    if (have_overlaps) then
      call mat_ao_to_mo('<parent|ion>',sx,sao)
      call report_mo_correlations(sx)
    end if
    if (same_mos) then
      write (out,"('Setting MO overlaps to a unit matrix')")
      call mat_unit(sx)
    endif
    if (use_symmetry) then
      call bm_scan_matrix(sx,sx_ms,eps=eps_integral)
      call bm_print(sx_ms,'MO overlap')
      ! Be paranoid - check whether scan produced sensible result!
      call bm_verify(sx,sx_ms,'MO overlap block structure check',eps=eps_integral)
    end if
    if (have_operator) then
      convert_operator: do ic=1,op_components
        write (tag,"('<bra|',a,'-',i0,'|ket>')") trim(braket_operator), ic
        call mat_ao_to_mo(trim(tag),dipx(:,:,ic),dipao(:,:,ic))
      end do convert_operator
      if (use_symmetry) then
        call bm_scan_matrix(sum(abs(dipx),dim=3),dipx_ms,eps=eps_integral)
        call bm_print(dipx_ms,'MO property')
        ! Be paranoid - check whether scan produced sensible result!
        call bm_verify(sum(abs(dipx),dim=3),dipx_ms,'MO property block structure check',eps=eps_integral)
      end if
    end if
    call TimerStop('AO->MO')
    !
    call locate_frozen_core
    !
    !  Figure out how much distance between determinants still gives non-zero
    !  contribution. We use the same tolerances to figure out how many zero
    !  rows/columns we can survive without the result being zero. The latter
    !  test is independent of Slater rules, which provide a convenient shortcut
    !  for the spin symmetry when orbitals are the same
    !
    if (same_mos) then
      tolerance = 0
      if (.not.mode_overlap) then
        tolerance = tolerance + 2
      end if
      if (have_operator) tolerance = tolerance + 2
      write (out,"('This calculation has determinant distance tolerance of ',i2)") tolerance
    end if
    drop_row = 0 ; rep_row = 0 ; drop_col = 0 ; rep_col = 0 ;
    if (have_operator) rep_row = 1
    if (.not.mode_overlap) drop_row = 1
    !
    !  For each pair of determinants in the parent and the ion, calculate
    !  the overlap - these will form contributions to the final Dyson
    !  orbital.
    !
    call TimerStop('Initialization')
    write (out,"(/'Done with preliminaries, beginning the calculation'/)")
    call TimerStart('Calculation')
    wdyson = 0._rk
    dipall = 0._rk
    det_included = 0
    det_cutoff   = 0
    det_distant  = 0
    det_earlyzero= 0
    !
    !$omp parallel default(none) &
    !$omp& private(detref,detion,sdet,dipdet,wdysondet,overterm,dipterm) &
    !$omp& private(distance,wdyson_thread,dipall_thread) &
    !$omp& private(sdet_ms) &
    !$omp& shared(tolerance,nelbra,nelket,ndetket,cdetion,cdetref) &
    !$omp& shared(eps_cdet,same_mos,occbra,occket,mode_overlap,have_operator) &
    !$omp& shared(nmobra,ndetbra,wdyson,verbose,op_components,dipall) &
    !$omp& shared(drop_row,drop_col,rep_row,rep_col,use_symmetry) &
    !$omp& reduction(+:det_cutoff,det_included,det_distant,det_earlyzero,overall)
    !
    !  Allocate thread-private storage now
    !
    allocate (sdet(nelbra,nelket), &
              dipdet(nelbra,nelket,op_components), &
              wdysondet(nelbra,1+op_components), &
              wdyson_thread(nmobra,2,1+op_components), &
              dipterm(op_components),dipall_thread(op_components))
    wdyson_thread = 0._rk
    dipall_thread = 0._rk
    sdet_ms%n_blocks = -1  ! This should not be necessary, but gfortran OpenMP is weird ...
    !$omp do schedule(dynamic)
    dyson_ref: do detref=1,ndetbra
      call report_progress(detref)
      dyson_ion: do detion=1,ndetket
        if ( abs(cdetion(detion)*cdetref(detref)) <= eps_cdet ) then
          det_cutoff = det_cutoff + 1
          cycle dyson_ion
        end if
        !
        !  Try a short-cut for the same MOs
        !
        if (same_mos) then
          distance = determinant_distance(occbra(:,detref),occket(:,detion))
          if (distance>tolerance) then
            det_distant = det_distant + 1
            cycle dyson_ion
          end if
        end if
        !
        if (.not.pick_determinant_integrals(drop_row,rep_row,drop_col,rep_col,detref,detion,sdet,dipdet,sdet_ms)) then 
          !
          !  We know this determinant will contribute only zeros; drop it
          !  Note that pick_determinant_integrals() will not copy the integrals
          !  if it figures out the result is zero; therefore continuing past
          !  this point is a mistake.
          !
          det_earlyzero = det_earlyzero + 1
          cycle dyson_ion
        end if
        !
        det_included = det_included + 1
        !
        !
        !  Now calculation takes different path for the overlaps and Dyson orbitals
        !
        if (mode_overlap) then
          overterm = determinant_overlap(sdet,sdet_ms)
          if (verbose>=2) then
            write (out,"(' ref det. = ',i5,' ion det. = ',i5,' weight = ',g16.9,' overlap = ',g16.9)") &
                   detref, detion, cdetion(detion)*cdetref(detref), overterm
          end if
          overall = overall + cdetion(detion)*cdetref(detref)*overterm
          !
          if (have_operator) then
            get_dipterms: do ic=1,op_components
              dipterm(ic) = determinant_transition(sdet,dipdet(:,:,ic),sdet_ms)
            end do get_dipterms
            if (verbose>=2) then
              write (out,"(' ref det. = ',i5,' ion det. = ',i5,' weight = ',g16.9,' dipole = ',3g16.9)") &
                     detref, detion, cdetion(detion)*cdetref(detref), dipterm
            end if
            dipall_thread = dipall_thread + cdetion(detion)*cdetref(detref)*dipterm
          end if
        else ! .not.mode_overlap
          !
          !  This determinant's contribution to the Dyson orbital
          !
          call determinant_dyson_orbital(sdet,sdet_ms,wdysondet)
          if (have_operator) then
            call determinant_exchange_orbitals(sdet,sdet_ms,dipdet,wdysondet)
          end if
          !
          !  Report partial amplitude
          !
          ! write (out,"('Ref. det. ',i5,' (',f10.5,') and ion det. ',i5,' (',f10.5,') contribute ',f12.8,' to cross-section')") &
          !        detref, cdetref(detref), detion, cdetion(detion), (cdetion(detion)*cdetref(detref))**2*sum(wdysondet(:)**2)
          !
          !  Debug output
          !
          if (verbose>=2) then
            write (out,"(' ref det. = ',i5,' ion det. = ',i5,' weight = ',g16.9)") &
                   detref, detion, cdetion(detion)*cdetref(detref)
            write (out,"(' ref = ',80i2)") occbra(:,detref)
            write (out,"(' ion = ',80i2)") occket(:,detion)
            write (out,"(('dyson:',8(1x,f10.6)))") wdysondet(:,1)
            print_wd: do ic=1,op_components
              write (out,"('exc ',i0,':')",advance='no') ic
              write (out,"(t7,8(1x,f10.6))") wdysondet(:,1+ic)
            end do print_wd
          end if
          !
          !  Accumulate contributions to the total Dyson orbital.
          !
          call accumulate_dyson(cdetion(detion)*cdetref(detref),occbra(:,detref),wdysondet,wdyson_thread)
        end if ! mode_overlap
        !
        !  Clean up auxiliary data
        !
        if (use_symmetry) then 
          call bm_destroy(sdet_ms)
        end if
      end do dyson_ion
    end do dyson_ref
    !$omp end do
    !
    !  Accumulate contributions to the total Dyson orbital. Accumulating arrays
    !  is not thread-safe, and must be inside a critical section. I wish OpenMP
    !  allowed reduction operators for deferred-shape arrays ...
    !
    !$omp critical
    if (mode_overlap) then
      dipall = dipall + dipall_thread
    else ! .not.mode_overlap
      wdyson = wdyson + wdyson_thread
    end if
    !$omp end critical
    !
    deallocate (sdet,dipdet,wdysondet,wdyson_thread,dipall_thread,dipterm)
    !$omp end parallel
    call TimerStop('Calculation')
    !
    call TimerStart('Output')
    !
    !  Cut-off statistics
    !
    write (out,"('          Determinant pairs processed: ',f23.0)") det_included
    write (out,"('  Determinant pairs cut off by weight: ',f23.0)") det_cutoff
    write (out,"('Determinant pairs cut off by distance: ',f23.0)") det_distant
    write (out,"('Determinant pairs cut off by symmetry: ',f23.0)") det_earlyzero
    !
    if (mode_overlap) then
      overall = sign_bra * sign_ket * overall
      write (out,"('Total overlap is ',g16.9)") overall
      if (have_operator) then
        dipall = sign_bra * sign_ket * dipall
        select case (braket_operator)
          case ('dipole')
            write (out,"('Transition dipole is',3(1x,f20.12))") dipall
          case default
            write (out,"('Matrix element of ',a,' operator is:'/(3(1x,f20.12)))") trim(braket_operator), dipall
        end select
      end if
    else 
      wdyson = sign_bra * sign_ket * wdyson
      !
      !  Calculate AO representation of the Dyson orbital
      !
      transform_dyson: do ic=1,1+op_components
        aodyson(:,1,ic) = matmul(cbra,wdyson(:,1,ic))
        aodyson(:,2,ic) = matmul(cbra,wdyson(:,2,ic))
      end do transform_dyson
      !
      !  Print our results. We'll do it two ways:
      !  a) Expansion over the MOs, which is good for understanding
      !  b) Expansion over the AOs, which gives a compact expression for
      !     numerical use.
      !
      call print_mo_results
      call print_ao_results
      call punch_ao_results
    end if ! mode_overlap
    call TimerStop('Output')
    call TimerStop('Superdyson')
    call TimerReport
  end subroutine call_superdyson
  !
  !  Prepare 1-particle integrals for a given determinant pair
  !
  function pick_determinant_integrals(drop_row,rep_row,drop_col,rep_col,detref,detion,sdet,dipdet,sdet_ms) result(success)
    integer(ik), intent(in)           :: drop_row, rep_row ! How many rows may be dropped and/or replaced?
    integer(ik), intent(in)           :: drop_col, rep_col ! ditto for columns
    integer(ik), intent(in)           :: detref, detion
    real(rk), intent(out)             :: sdet (:,:)
    real(rk), intent(out)             :: dipdet(:,:,:)
    type(bm_structure), intent(inout) :: sdet_ms         ! Block descriptor of sdet
    logical                           :: success
    !
    integer(ik)             :: ic
    !
    !  Try to bail out early by using block structure of smat
    !
    if (use_symmetry) then
      call analyze_spinorbit_integral_blocks(detref,detion,sx_ms,sdet_ms)
      if (verbose>=2) call bm_print(sdet_ms,'sdet')
      if (hopeless_blocks(drop_row,rep_row,drop_col,rep_col,sdet_ms)) then
        !
        !  All our manipulations below will not lead to a non-zero determinant; bail out now.
        !
        call bm_destroy(sdet_ms)
        success = .false.
        return
      end if
    end if
    !
    !  Spin-orbit 1-particle overlaps and dipoles for this determinant
    !
    call fill_spinorbit_integrals(occbra(:,detref),occket(:,detion),sdet,sx)
    if (have_operator) then
      fill_components: do ic=1,op_components
        call fill_spinorbit_integrals(occbra(:,detref),occket(:,detion),dipdet(:,:,ic),dipx(:,:,ic))
      end do fill_components
    end if
    !
    if (verbose>=3) then
      write (out,"('1-e overlap matrix for dets ',2i5)") detref, detion
      call print_matrix(sdet)
      if (have_operator) then
        print_components: do ic=1,op_components
          write (out,"(1x,a,'-',i0,'matrix for dets ',2i5)") trim(braket_operator), ic, detref, detion
          call print_matrix(dipdet(:,:,ic))
        end do print_components
      end if
    end if
    success = .true.
  end function pick_determinant_integrals
  !
  !  Calculate coefficients of a Dyson orbital
  !
  subroutine determinant_dyson_orbital(sdet,sdet_ms,wdysondet)
    real(rk), intent(in)           :: sdet(:,:)      ! 1-particle overlaps of spin-orbitals
    type(bm_structure), intent(in) :: sdet_ms        ! Block-structure of sdet
    real(rk), intent(inout)        :: wdysondet(:,:) ! ne determinant's contribution to the Dyson/exchange orbital
    !
    integer(ik) :: imo                ! Spin-orbital index
    real(rk)    :: scr(nelket,nelket) ! Buffer for the 1-particle overlap matrix 
    !
    call TimerStart('determinant_dyson_orbital')
    if (use_symmetry) then
      !
      !  We have two important special cases, which are advantageous to handle separately
      !
      wdysondet(:,1) = 0
      select case (sdet_ms%max_row-sdet_ms%n_rows)
        case default ! More than one zero row means all determinants will be zero, exit right away
        case (0)
          accumulate_sym: do imo=1,nelbra
            call drop_row_determinant(imo)
          end do accumulate_sym
        case (1)
          imo = bm_find_zero_row(sdet_ms)
          if (imo<=0) stop 'superdyson%determinant_dyson_orbital - no zero row in a defective matrix?!'
          call drop_row_determinant(imo)
      end select
    else
      mo_weights: do imo=1,nelbra
        call drop_row(imo,sdet,scr)
        wdysondet(imo,1) = ((-1)**(nelbra+imo))*determinant(scr)
      end do mo_weights
    end if
    call TimerStop('determinant_dyson_orbital')
    !
    contains
    !
    subroutine drop_row_determinant(imo)
      integer(ik), intent(in) :: imo  ! Row to drop from the sdet() matrix
      !
      type(bm_structure) :: scr_ms    ! Structure of the matrix with deleted row
      !
      call bm_drop_row(sdet_ms,imo,scr_ms)
      if (hopeless_blocks(0_ik,0_ik,0_ik,0_ik,scr_ms)) then
        call bm_destroy(scr_ms)
        return
      end if
      call drop_row(imo,sdet,scr)
      ! call bm_verify(scr,scr_ms,eps=eps_integral,text="drow_row")
      wdysondet(imo,1) = ((-1)**(nelbra+imo))*block_det(scr,scr_ms)
      call bm_destroy(scr_ms)
    end subroutine drop_row_determinant
  end subroutine determinant_dyson_orbital
  !
  !  Calculate coefficients of the exchange correction orbitals
  !
  subroutine determinant_exchange_orbitals(sdet,sdet_ms,dipdet,wdysondet)
    real(rk), intent(in)           :: sdet(:,:)      ! 1-particle overlaps of spin-orbitals
    type(bm_structure), intent(in) :: sdet_ms        ! Block-structure of sdet
    real(rk), intent(in)           :: dipdet(:,:,:)  ! Spin-orbital dipole integrals
    real(rk), intent(inout)        :: wdysondet(:,:) ! ne determinant's contribution to the Dyson/exchange orbital
    !
    integer(ik)        :: imo                  ! Spin-orbital index
    integer(ik)        :: ic                   ! Component index
    real(rk)           :: scr(nelket,nelket,2) ! Buffer for the 1-particle overlap and dipole matrices
    type(bm_structure) :: bms
    !
    call TimerStart('determinant_exchange_orbitals')
    mo_weights: do imo=1,nelbra
      if (use_symmetry) then
        wdysondet(imo,2:) = 0
        call bm_drop_row(sdet_ms,imo,bms)
        if (hopeless_blocks(1_ik,1_ik,0_ik,0_ik,bms)) then
          call bm_destroy(bms)
          cycle mo_weights
        end if
      end if ! use_symmetry
      call drop_row(imo,sdet,scr(:,:,1))
      process_components: do ic=1,op_components
        call drop_row(imo,dipdet(:,:,ic),scr(:,:,2))
        wdysondet(imo,1+ic) = ((-1)**(nelbra+imo))*determinant_transition(scr(:,:,1),scr(:,:,2),bms)
      end do process_components
      if (use_symmetry) call bm_destroy(bms)
    end do mo_weights
    call TimerStop('determinant_exchange_orbitals')
    !
  end subroutine determinant_exchange_orbitals
  !
  !  Add contributions to the total Dyson and exchange orbital
  !
  subroutine accumulate_dyson(wgt,occbra,wdysondet,wdyson)
    real(rk), intent(in)    :: wgt            ! Weight of this contribution
    integer(ik), intent(in) :: occbra(:)      ! Parent determinant composition, bra
    real(rk), intent(in)    :: wdysondet(:,:) ! One determinant's contribution to the Dyson/exchange orbital
    real(rk), intent(inout) :: wdyson(:,:,:)  ! Thread-private copy of wdyson
    !
    integer(ik) :: moref   ! Spatial MO indices
    integer(ik) :: spinref ! Spin+space MO indices
    !
    call TimerStart('accumulate_dyson')
    !
    !  Spin-alpha part
    !
    spinref = 0
    ref_alpha: do moref=1,nmobra
      if (occbra(moref)/=1 .and. occbra(moref)/=2) cycle ref_alpha
      !
      spinref = spinref + 1
      wdyson(moref,1,:) = wdyson(moref,1,:) + wgt*wdysondet(spinref,:)
      if (verbose>=2 .and. any(abs(wdysondet(spinref,:))>=1e-6_rk)) then
         write (out,"(' alpha MO ',i5,' gets ',4(f16.10,1x),' from spin-orbital ',i5)") &
                moref, wgt*wdysondet(spinref,:), spinref
      end if
    end do ref_alpha
    !
    !  Spin-beta part
    !
    ref_beta: do moref=1,nmobra
      if (occbra(moref)/=-1 .and. occbra(moref)/=2) cycle ref_beta
      !
      spinref = spinref + 1
      wdyson(moref,2,:) = wdyson(moref,2,:) + wgt*wdysondet(spinref,:)
      if (verbose>=2 .and. any(abs(wdysondet(spinref,:))>=1e-6_rk)) then
         write (out,"('  beta MO ',i5,' gets ',4(f16.10,1x),' from spin-orbital ',i5)") &
                moref, wgt*wdysondet(spinref,:), spinref
      end if
    end do ref_beta
    if (spinref/=nelbra) stop 'superdyson%accumulate_dyson - counting error'
    call TimerStop('accumulate_dyson')
  end subroutine accumulate_dyson
  !
  !  Print results in the MO basis
  !
  subroutine print_mo_results
    integer(ik) :: imo
    !
    write (out,"()")
    write (out,"(t10,'-----------------------------------------------------')")
    write (out,"(t10,'Dyson orbital, expanded in MOs of the parent molecule')")
    write (out,"(t10,'-----------------------------------------------------')")
    write (out,"()")
    write (out,"(t5,a5,t12,8(a10,1x))") 'MO', 'C(alpha)', 'C(beta)', 'X(alpha)', 'X(beta)', &
                                              'Y(alpha)', 'Y(beta)', 'Z(alpha)', 'Z(beta)'
    write (out,"(t5,a5,t12,8(a10,1x))") '--', '--------', '-------', '--------', '-------', &
                                              '--------', '-------', '--------', '-------'
    print_coefficients: do imo=1,nmobra
      write (out,"(t5,i5,t12,8(f10.5,1x))") imo, wdyson(imo,:,:)
    end do print_coefficients
    amplitude = sum(wdyson(:,:,1)**2)
    write (out,"()")
    write (out,"('<psid|psid> = ',f14.8)") amplitude
    write (out,"()")
    !
  end subroutine print_mo_results
  !
  !  Print results in the AO basis
  !
  subroutine print_ao_results
    integer(ik) :: ibas
    !
    write (out,"()")
    write (out,"(t10,'------------------------------')")
    write (out,"(t10,'Dyson orbital, expanded in AOs')")
    write (out,"(t10,'------------------------------')")
    write (out,"()")
    write (out,"(t5,a5,1x,a14,1x,t26,8(a11,1x))") 'AO', '', 'C(alpha)', 'C(beta)', 'X(alpha)', 'X(beta)', &
                                                            'Y(alpha)', 'Y(beta)', 'Z(alpha)', 'Z(beta)'
    write (out,"(t5,a5,1x,a14,1x,t26,8(a11,1x))") '--', '', '--------', '-------', '--------', '-------', &
                                                            '--------', '-------', '--------', '-------'
    print_coefficients: do ibas=1,naos_bra
      write (out,"(t5,i5,a14,1x,t26,8(f11.5,1x))") ibas, aolabels(ibas), aodyson(ibas,:,:)
    end do print_coefficients
    write (out,"()")
    !
  end subroutine print_ao_results
  !
  subroutine punch_ao_results
    integer(ik) :: ic
    !
    write (out,"('--- DYSON ORBITALS ---'/1x,a/'<psid|psid> = ',f14.8)") &
           trim(comment), amplitude
    write (out,"(' $VECDYS')")
    !
    punch_terms: do ic=1,1+op_components
      call punch_vector(1+2*(ic-1),aodyson(:,1,ic))
      call punch_vector(2+2*(ic-1),aodyson(:,2,ic))
    end do punch_terms
    !
    write (out,"(' $END   ')")
    !
  end subroutine punch_ao_results
end module superdyson
