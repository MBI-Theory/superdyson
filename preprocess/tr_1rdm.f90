!
!  14 March 2009, SP     - Initial version, based on "superdyson.f90"
!  09 September 2009, SP - Exploit frozen core where practical
!
!  Calculation of reduced 1-particle density matrices for transitions
!  between two general CI wavefunctions. Regular 1-RDM can be done as
!  a degenerate case.
!
!  Note that we calculate total 1-RDM, rather than spin-RDM
!
!  The side efect is that UHF wavefunctions are not acceptable on
!  input - however, R(O)HF wavefunctions are perfectly OK.
!
!  The implementation is _amazingly_ inefficient.
!
!  See read_core_input code in sd_core.f90 for the ad-hoc description of the input
!
module tr_1rdm
  use sd_core
  implicit none
  private
  public call_tr_1rdm
  public rcsid_tr_1rdm
  !
  character(len=clen), save :: rcsid_tr_1rdm = "$Id: tr_1rdm.f90,v 1.21 2021/09/29 13:43:22 ps Exp $"
  !
  !  Local data for dyson
  !
  real(rk)                :: det_included     ! Number of determinant pairs included
  real(rk)                :: det_cutoff       ! Number of determinant pairs cut off
  real(rk)                :: det_distant      ! Number of too distant determinant pairs 
  real(rk)                :: det_earlyzero    ! Number of too distant determinant pairs 
  !
  !  Our results. The last index is spin
  !
  real(rk), allocatable   :: rdm_mo(:,:)      ! Total 1-RDM, MO basis
  integer(ik)             :: rdm_sv_size      ! Number of singular values of the RDM
                                              ! rdm_sv_size = min(nmobra,nmoket)
  real(rk), allocatable   :: rdm_sv(:)        ! Singular values of the RDM - same for MO & AO
  real(rk), allocatable   :: rdm_lv_mo(:,:)   ! Left singular vectors of the MO RDM
  real(rk), allocatable   :: rdm_rv_mo(:,:)   ! Right singular vectors of the MO RDM
  real(rk), allocatable   :: rdm_lv_ao(:,:)   ! Left singular vectors of the AO RDM
  real(rk), allocatable   :: rdm_rv_ao(:,:)   ! Right singular vectors of the AO RDM
  !
contains
  !
  subroutine call_tr_1rdm
    integer(ik)              :: detref, detion
    real(rk), allocatable    :: sdet  (:,:)           ! 1-particle overlaps of spin-orbitals
    real(rk), allocatable    :: detrdm(:,:)           ! Partial spin- 1-RDM for one determinant
    real(rk), allocatable    :: rdm_mo_thread(:,:)    ! Thread-private buffer for 1-RDM accumulation
    integer(ik), allocatable :: indref(:), indion(:)  ! Mapping between determinant MOs and 
                                                      ! the active MO space
    integer(ik)              :: distance              ! Distance between the determinants
    type(bm_structure)       :: sdet_ms               ! Block-descriptor of sdet
    !
    call TimerStart('Transition 1RDM')
    call TimerStart('Initialization')
    !
    !  Done reading input, begin calculations
    !
    call check_input_consistency
    call flush(out)
    !
    if (.not.mode_overlap) then
      write (out,"('Number of electrons in two wavefunctions is not the same')") 
      stop 'tr_1rdm - number of electrons does not match'
    end if
    !
    rdm_sv_size = min(nmobra,nmoket)
    allocate (rdm_mo(nmobra,nmoket))
    allocate (rdm_sv(rdm_sv_size))
    allocate (rdm_lv_mo(nmobra,rdm_sv_size),rdm_rv_mo(nmoket,rdm_sv_size))
    allocate (rdm_lv_ao(naos_bra,rdm_sv_size),rdm_rv_ao(naos_ket,rdm_sv_size))
    !
    !  Transform integrals to the MO basis
    !
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
    !
    call locate_frozen_core
    !
    call TimerStop('Initialization')
    write (out,"(/'Done with preliminaries, beginning the calculation'/)")
    call TimerStart('Calculation')
    call flush(out)
    !
    det_included = 0
    det_cutoff   = 0
    det_distant  = 0
    det_earlyzero= 0
    rdm_mo = 0._rk
    !$omp parallel default(none) &
    !$omp& private(detref,detion,sdet,indref,indion,detrdm,distance,rdm_mo_thread) &
    !$omp& private(sdet_ms) &
    !$omp& shared(nelbra,nelket,ndetket,cdetref,cdetion,eps_cdet) &
    !$omp& shared(occbra,occket,sx,sx_ms,same_mos,nmobra,nmoket,rdm_mo) &
    !$omp& shared(ndetbra,verbose,use_symmetry) &
    !$omp& reduction(+:det_cutoff,det_included,det_distant,det_earlyzero)
    !
    !  Allocate thread-private storage now
    !
    allocate (sdet(nelbra,nelket),indref(nelbra),indion(nelket))
    allocate (detrdm(nelbra,nelket),rdm_mo_thread(nmobra,nmoket))
    rdm_mo_thread    = 0
    sdet_ms%n_blocks = -1  ! This should not be necessary, but gfortran OpenMP is weird ...
    !$omp do schedule(dynamic)
    det_ref: do detref=1,ndetbra
      call report_progress(detref)
      !
      det_ion: do detion=1,ndetket
        if ( abs(cdetion(detion)*cdetref(detref)) <= eps_cdet ) then
          det_cutoff = det_cutoff + 1
          cycle det_ion
        end if
        !
        !  Try to short-circuit most of the zeros for the case of identical
        !  reference MOs. If the distance between the two deteminants is
        !  more than two (meaning a single-orbital substitution), the pair 
        !  won't contribute to the RDM, and should be discarded.
        !
        if (same_mos) then
          distance = determinant_distance(occbra(:,detref),occket(:,detion))
          if (distance>2) then
            det_distant = det_distant + 1
            cycle det_ion
          end if
        end if
        !
        if (use_symmetry) then
          call analyze_spinorbit_integral_blocks(detref,detion,sx_ms,sdet_ms)
          if (verbose>=2) call bm_print(sdet_ms,'sdet')
          if (hopeless_blocks(1_ik,0_ik,1_ik,0_ik,sdet_ms)) then
            !
            !  All our manipulations below will not lead to a non-zero determinant; skip to next determinant
            !
            call bm_destroy(sdet_ms)
            det_earlyzero = det_earlyzero + 1
            cycle det_ion
          end if
        end if
        !
        det_included = det_included + 1
        !
        call fill_spinorbit_integrals(occbra(:,detref),occket(:,detion),sdet,sx)
        call fill_spinorbit_indices(occbra(:,detref),indref)
        call fill_spinorbit_indices(occket(:,detion),indion)
        !
        call determinant_rdm(sdet,sdet_ms,detrdm)
        if (use_symmetry) call bm_destroy(sdet_ms)
        if (verbose>=3) then
          write (out,"(' ref det. = ',i5,' ion det. = ',i5,' weight = ',g16.9)") &
                 detref, detion, cdetion(detion)*cdetref(detref)
          write (out,"(' Partial spin 1-RDM:')")
          call print_matrix(detrdm)
        end if
        !
        call accumulate_rdm(cdetion(detion)*cdetref(detref),indref,indion,detrdm,rdm_mo_thread)
      end do det_ion
    end do det_ref
    !$omp end do
    !$omp critical
    rdm_mo = rdm_mo + rdm_mo_thread
    !$omp end critical
    deallocate (sdet,indref,indion,detrdm,rdm_mo_thread)
    !$omp end parallel
    rdm_mo = sign_bra * sign_ket * rdm_mo
    call TimerStop('Calculation')
    !
    call TimerStart('Output')
    !
    !  Cut-off statistics
    !
    write (out,"('          Determinant pairs processed: ',f13.0)") det_included
    write (out,"('  Determinant pairs cut off by weight: ',f13.0)") det_cutoff
    write (out,"('Determinant pairs cut off by distance: ',f13.0)") det_distant
    write (out,"('Determinant pairs cut off by symmetry: ',f13.0)") det_earlyzero
    !
    write (out,"(/t5,'1-RDM in the MO basis'/)")
    call print_matrix(rdm_mo)
    !
    !  Analyze principal components of the 1-RDM 
    !
    call analyze_rdm
    !
    !  Transform 1-RDM principal components into the AO basis
    !
    rdm_lv_ao = matmul(cbra,rdm_lv_mo)
    rdm_rv_ao = matmul(cket,rdm_rv_mo)
    !
    !  Print our results. We'll do it two ways:
    !  a) Expansion over the MOs, which is good for understanding
    !  b) Principal-component expansion over the AOs, which gives 
    !     a compact expression for numerical use (and a different
    !     way of looking at the data, too)
    !
    call print_mo_results
    call print_ao_results
    call punch_ao_results
    !
    call TimerStop('Output')
    call TimerStop('Transition 1RDM')
    call TimerReport
  end subroutine call_tr_1rdm
  !
  !  Calculate spin 1-RDM for a given determinant
  !
  subroutine determinant_rdm(sdet,sdet_ms,detrdm)
    real(rk), intent(in)           :: sdet(:,:)   ! 1-particle overlaps of spin-orbitals
    type(bm_structure), intent(in) :: sdet_ms     ! Block structure of sdet
    real(rk), intent(out)          :: detrdm(:,:) ! spin 1-RDM
    !
    integer(ik)        :: mo_ref, mo_ion         ! Spin-orbital index
    real(rk)           :: scr(nelbra-1,nelket-1) ! Buffer for the minors
    type(bm_structure) :: dr_ms                  ! "drop row" matrix structure
    type(bm_structure) :: dc_ms                  ! "drop column" matrix structure
    !
    call TimerStart('determinant_rdm')
    !
    !  Handle degenerate case separately
    !
    if (nelbra==1 .and. nelket==1) then
      detrdm = 1
    else
      if (use_symmetry) then
        !
        !  It is not a good idea to screen minors directly: the overhead
        !  of screening is quite large, and we loose the change of killing
        !  off entire rows. Instead, we'll drop a row first, screen, then
        !  drop a column. Obviously, there is no need to actually copy the
        !  integrals - these can be done once we decide a minor is actually
        !  viable.
        !
        detrdm = 0
        sym_mo_ref_: do mo_ref=1,nelbra
          call bm_drop_row(sdet_ms,mo_ref,dr_ms)
          if (.not.hopeless_blocks(0_ik,0_ik,1_ik,0_ik,dr_ms)) then
            sym_mo_ion_: do mo_ion=1,nelket
              call bm_drop_column(dr_ms,mo_ion,dc_ms)
              if (.not.hopeless_blocks(0_ik,0_ik,0_ik,0_ik,dc_ms)) then
                call get_minor(mo_ref,mo_ion,sdet,scr)
                detrdm(mo_ref,mo_ion) = ((-1)**(mo_ref+mo_ion))*block_det(scr,dc_ms)
              end if
              call bm_destroy(dc_ms)
            end do sym_mo_ion_
          end if
          call bm_destroy(dr_ms)
        end do sym_mo_ref_
      else
        mo_ion_: do mo_ion=1,nelket
          mo_ref_: do mo_ref=1,nelbra
            call get_minor(mo_ref,mo_ion,sdet,scr)
            detrdm(mo_ref,mo_ion) = ((-1)**(mo_ref+mo_ion))*determinant(scr)
          end do mo_ref_
        end do mo_ion_
      end if
    end if
    call TimerStop('determinant_rdm')
    !
  end subroutine determinant_rdm
  !
  !  Add determinant's contribution to the total 1-RDM
  !
  subroutine accumulate_rdm(wgt,indref,indion,detrdm,rdm_mo)
    real(rk), intent(in)    :: wgt         ! Weight of this contribution
    integer(ik), intent(in) :: indref(:)   ! Index table for the left-hand side
    integer(ik), intent(in) :: indion(:)   ! Index table for the right-hand side
    real(rk), intent(in)    :: detrdm(:,:) ! Spin-RDM for the current determinant
    real(rk), intent(inout) :: rdm_mo(:,:) ! 1-RDM accumulation buffer; can be thread-private
    !
    integer(ik) :: mo_det_ref, mo_det_ion, mo_ref, mo_ion
    integer(ik) :: aion, aref
    !
    call TimerStart('accumulate_rdm')
    !
    !  Transformation below is actually a reduction - detrdm() is indexed by
    !  spin-orbitals, while rdm() is indexed by spatial MOs alone. Hence an
    !  array expresion would NOT be a good idea here!
    !
    !  Note that only same-spin block (alpha-alpha and beta-beta) should be
    !  accumulated. Since our operators are spin-free, cross-spin blocks
    !  will vanish upon spin integration, and do not contribute
    !
    det_ion: do mo_det_ion=1,size(detrdm,dim=2)
      mo_ion = indion(mo_det_ion)
      aion   = abs(mo_ion)
      det_ref: do mo_det_ref=1,size(detrdm,dim=1)
        mo_ref = indref(mo_det_ref)
        if ( (mo_ion>0.and.mo_ref<0) .or. (mo_ion<0.and.mo_ref>0) ) cycle det_ref
        aref   = abs(mo_ref)
        rdm_mo(aref,aion) = rdm_mo(aref,aion) + wgt*detrdm(mo_det_ref,mo_det_ion)
      end do det_ref
    end do det_ion
    !
    call TimerStop('accumulate_rdm')
  end subroutine accumulate_rdm
  !
  !  A slightly more sofisticated version, which will try to preserve
  !  block-diagonality of the 1-RDM representation. Within the numerical
  !  accuracy, it should produce results identical to the old version
  !
  subroutine analyze_rdm
    integer(ik)           :: ivec, ic, ir
    real(rk)              :: cmax, tmp
    real(rk), allocatable :: a(:,:), vth(:,:)
    !
    allocate (a(nmobra,nmoket),vth(rdm_sv_size,nmoket))
    !
    a = rdm_mo
    ! call lapack_svd(a,rdm_sv,rdm_lv_mo,vth)
      call block_gesvd(a,rdm_sv,rdm_lv_mo,vth)
    rdm_rv_mo = transpose(vth)
    deallocate (a,vth)
    !
    !  Sanity check - we must have a valid decomposition!
    !
    cmax = 1e3_rk*spacing(maxval(abs(rdm_mo)))
    decomposition_check_column: do ic=1,nmoket
      decomposition_check_row: do ir=1,nmobra
        tmp = sum(rdm_lv_mo(ir,:)*rdm_sv(:)*rdm_rv_mo(ic,:))
        if (abs(rdm_mo(ir,ic)-tmp)<=cmax) cycle decomposition_check_row
        write (out,"('1-RDM element ',2i8,' must be ',g24.15,' but is ',g24.15,' error = ',g24.15)") &
               ir, ic, rdm_mo(ir,ic), tmp, rdm_mo(ir,ic) - tmp
        stop 'tr_1rdm%analyze_rdm - faulty decomposition'
      end do decomposition_check_row
    end do decomposition_check_column
    write (out,"(/'Singular value decomposition accuracy check passed, eps = ',g13.5/)") cmax
    !
    !  Simultaneously flipping the signs of both the left and right singular
    !  vectors does not change the overall result; however, it is nice to 
    !  produce consistent results - so let's try to make the (first copy of the) 
    !  largest coeffient of the left singular vector positive.
    !
    scan_for_signs: do ivec=1,rdm_sv_size
      cmax = maxval(abs(rdm_lv_mo(:,ivec)))
      cmax = cmax - 100._rk * spacing(cmax) ! Give a margin for numerical noise
      find_first_coefficient: do ic=1,nmobra
        if (abs(rdm_lv_mo(ic,ivec))<cmax) cycle find_first_coefficient
        !
        if (rdm_lv_mo(ic,ivec)<0) then
          rdm_lv_mo(:,ivec) = -rdm_lv_mo(:,ivec)
          rdm_rv_mo(:,ivec) = -rdm_rv_mo(:,ivec)
        end if
        cycle scan_for_signs
      end do find_first_coefficient
      stop 'tr_1rdm%analyze_rdm - can''t locate the maximum?!'
    end do scan_for_signs
  end subroutine analyze_rdm
  !
  !  Print results in the MO basis
  !
  subroutine print_mo_results
    integer(ik) :: imo, ic
    !
    write (out,"()")
    write (out,"(t10,'--------------------------------------------')")
    write (out,"(t10,'Principal components of the 1-RDM (MO basis)')")
    write (out,"(t10,'--------------------------------------------')")
    write (out,"()")
    component_loop: do ic=1,rdm_sv_size
      if (abs(rdm_sv(ic))<spacing(100._rk)) cycle component_loop ! Skip the null-space
      !
      write (out,"('Component ',i5,' singular value ',g16.10/)") ic, rdm_sv(ic)
      !
      write (out,"(t5,a5,t11,2(a20,1x))") 'MO', 'Left S.V.', 'Right S.V.'
      write (out,"(t5,a5,t11,2(a20,1x))") '--', '---------', '----------'
      print_coefficients: do imo=1,max(nmobra,nmoket)
        if (imo<=min(nmobra,nmoket)) then
          write (out,"(t5,i5,t11,f20.12,1x,f20.12)") imo, rdm_lv_mo(imo,ic), rdm_rv_mo(imo,ic)
        else if (imo<=nmobra) then
          write (out,"(t5,i5,t11,f20.12,1x,a20   )") imo, rdm_lv_mo(imo,ic), ' '
        else
          write (out,"(t5,i5,t11,a20   ,1x,f20.12)") imo, ' '              , rdm_rv_mo(imo,ic)
        end if
      end do print_coefficients
      write (out,"()")
    end do component_loop
    !
    write (out,"(/1x,i5,' components of 1-RDM are in the null-space.'/)") count(rdm_sv<spacing(100._rk))
    !
  end subroutine print_mo_results
  !
  !  Print results in the AO basis
  !
  subroutine print_ao_results
    integer(ik) :: ibas, ic
    !
    write (out,"()")
    write (out,"(t10,'--------------------------------------------')")
    write (out,"(t10,'Principal components of the 1-RDM (AO basis)')")
    write (out,"(t10,'--------------------------------------------')")
    write (out,"()")
    component_loop: do ic=1,rdm_sv_size
      if (abs(rdm_sv(ic))<spacing(100._rk)) cycle component_loop ! Skip the null-space
      !
      write (out,"('Component ',i5,' singular value ',g16.10/)") ic, rdm_sv(ic)
      !
      write (out,"(t5,a5,a12,3x,t26,2(a20,1x))") 'AO', 'Type', 'Left S.V.', 'Right S.V.'
      write (out,"(t5,a5,a12,3x,t26,2(a20,1x))") '--', '----', '---------', '----------'
      print_coefficients: do ibas=1,min(naos_bra,naos_ket)
        write (out,"(t5,i5,a14,1x,t26,2(f20.12,1x))") &
               ibas, aolabels(ibas), rdm_lv_ao(ibas,ic), rdm_rv_ao(ibas,ic)
      end do print_coefficients
      !
      !  If basis sets are not the same, print the rest of the vectors separately
      !
      if (naos_bra>naos_ket) then
        print_coefficients_r0: do ibas=naos_ket+1,naos_bra
          write (out,"(t5,i5,a14,1x,t26,f20.12,1x)") &
                 ibas, aolabels(ibas), rdm_lv_ao(ibas,ic)
        end do print_coefficients_r0
      else ! naos_bra>naos_ket
        print_coefficients_0i: do ibas=naos_bra+1,naos_ket
          write (out,"(t5,i5,a14,1x,t26,21x,f20.12,1x)") &
                 ibas, aolabels(ibas), rdm_rv_ao(ibas,ic)
        end do print_coefficients_0i
      end if
      !
      write (out,"()")
    end do component_loop
    !
    write (out,"(/1x,i5,' components of 1-RDM are in the null-space.'/)") count(rdm_sv<spacing(100._rk))
    !
  end subroutine print_ao_results
  !
  subroutine punch_ao_results
    integer(ik) :: ic, nc
    !
    nc = count(rdm_sv>=spacing(100._rk))
    !
    write (out,"('--- TRANSITION 1-RDM SINGULAR VALUES ---'/1x,a)") trim(comment)
    write (out,"(' $RDMPCE')")
    write (out,"(5(1x,f14.12))") rdm_sv(:nc)
    write (out,"(' $END   ')")
    !
    write (out,"('--- TRANSITION 1-RDM LEFT/RIGHT SINGULAR VECTORS ---'/1x,a)") trim(comment)
    write (out,"(' $VECRDM')")
    !
    punch_terms: do ic=1,nc
      call punch_vector(1+2*(ic-1),rdm_lv_ao(:,ic))
      call punch_vector(2+2*(ic-1),rdm_rv_ao(:,ic))
    end do punch_terms
    !
    write (out,"(' $END   ')")
    !
  end subroutine punch_ao_results
end module tr_1rdm
