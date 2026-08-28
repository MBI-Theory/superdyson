!
!  1-electron integrals for complex operators.
!
!   subroutine gamess_1e_integrals_complex(what,v,bra,ket,op_xyz,op_param,op_index)
!     integer, parameter             :: ok = rk
      character(len=*), intent(in)   :: what         ! What to evaluate
      complex(ok), intent(out)       :: v(:,:)       ! Buffer for the results
      type(gam_structure), intent(in), target, optional :: bra
      type(gam_structure), intent(in), target, optional :: ket
      real(ok), intent(in), optional  :: op_xyz(:)   ! Position of the operator (3c integrals)
      real(ok), intent(in), optional  :: op_param(:) ! Additional parameters for the operator
                                                     ! Required size and allowed values depend on
                                                     ! the operator; see below.
      integer(ik), intent(in), optional :: op_index(:) ! Additional (integer) parameters for the operator.
      !
      character(len=40)            :: lwhat          ! Local copy of "what"
      type(gam_structure), pointer :: l, r
      type(gam_operator_data)      :: op_data        ! Operator data
      !
      call TimerStart('GAMESS '//trim(what))
      l => gam_def ; r => gam_def
      if (present(bra)) l => bra
      if (present(ket)) r => ket
      !
      !  What is it we want to evaluate?
      !
      lwhat = what
      select case (lwhat(1:3))
        case ('AO ')
          if ( (size(v,dim=1)<l%nbasis) .or. (size(v,dim=2)<r%nbasis) ) then
            write (out,"('import_gamess%gamess_1e_integrals_complex: Output buffer is not large enough')")
            stop 'import_gamess%gamess_1e_integrals_complex: Output buffer is not large enough'
          end if
      end select
      !
      !  3-centre integrals need position of the operator; make sure it was supplied
      !  Operator-specific parameters will be handled later.
      !
      if (lwhat(4:5)=='3C') then
        if (.not.present(op_xyz)) then
          write (out,"('import_gamess%gamess_1e_integrals_complex: operator location missing for ',a)") trim(lwhat)
          stop 'import_gamess%gamess_1e_integrals_complex: Missing operator position'
        end if
        op_data%op_xyz = real(op_xyz,kind=kind(op_data%op_xyz))
      end if
      !
      op_data%op_name = lwhat(4:)
      !
      select case (lwhat)
        case default
          write (out,"('import_gamess%gamess_1e_integrals_complex: ',a,' are not implemented')") trim(lwhat)
          stop 'import_gamess%gamess_1e_integrals_complex - integral not implemented'
        !
        !  First the integrals where no additional parmaters are needed:
        !
        !  'AO PLANEWAVE'   - Complex planewave operator [Exp(I K R)]; needs K vector
        !  'AO R PLANEWAVE' - Component of R (X, Y, or Z), times the planewave.
        !                     Neeeds K vector and the X/Y/Z index.
        !
        case ('AO PLANEWAVE')
          if (.not.present(op_param)) then
            write (out,"('import_gamess%gamess_1e_integrals_complex: exponent missing for ',a)") trim(lwhat)
            stop 'import_gamess%gamess_1e_integrals_complex: Missing operator parameter'
          end if
          if (size(op_param)/=3) then
            write (out,"('import_gamess%gamess_1e_integrals_complex: wrong number of parameters for ',a)") trim(lwhat)
            stop 'import_gamess%gamess_1e_integrals_complex: Wrong parameter count'
          end if
          op_data%kvec = real(op_param(1:3),kind=kind(op_data%kvec))
          call os_1e_matrix(op_data,l,r,v)
        case ('AO R PLANEWAVE')
          if (.not.present(op_param)) then
            write (out,"('import_gamess%gamess_1e_integrals_complex: exponent missing for ',a)") trim(lwhat)
            stop 'import_gamess%gamess_1e_integrals_complex: Missing operator parameter'
          end if
          if (size(op_param)/=3) then
            write (out,"('import_gamess%gamess_1e_integrals_complex: wrong number of parameters for ',a)") trim(lwhat)
            stop 'import_gamess%gamess_1e_integrals_complex: Wrong parameter count'
          end if
          op_data%kvec = real(op_param(1:3),kind=kind(op_data%kvec))
          if (.not.present(op_index)) then
            write (out,"('import_gamess%gamess_1e_integrals_complex: index missing for ',a)") trim(lwhat)
            stop 'import_gamess%gamess_1e_integrals_complex: Missing operator parameter'
          end if
          if (size(op_index)/=1) then
            write (out,"('import_gamess%gamess_1e_integrals_complex: wrong number of parameters for ',a)") trim(lwhat)
            stop 'import_gamess%gamess_1e_integrals_complex: Wrong parameter count'
          end if
          op_data%op_i = op_index(1)
          call os_1e_matrix(op_data,l,r,v)
      end select
      call TimerStop('GAMESS '//trim(what))
!   end subroutine gamess_1e_integrals_complex
