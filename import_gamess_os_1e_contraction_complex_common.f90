!
!  Calculate contraction of a primitive integral block (complex version)
!
!   subroutine os_1e_contraction_complex(what,vb,r_l,l_l,z_l,c_l,r_r,l_r,z_r,c_r)
!     integer, parameter :: ok = rk
      type (gam_operator_data), intent(in) :: what    ! Integral to evaluate
      complex(ok), intent(out)             :: vb(:,:) ! Place for the integrals
      real(ok), intent(in)                 :: r_l(:)  ! Coordinate on the left, in Bohr
      integer(ik), intent(in)              :: l_l     ! "Angular momentum" on the left
      real(ok), intent(in)                 :: z_l(:)  ! Primitive exponents on the left
      real(ok), intent(in)                 :: c_l(:)  ! Contractions on the left
      real(ok), intent(in)                 :: r_r(:)  ! Ditto on the right
      integer(ik), intent(in)              :: l_r     ! 
      real(ok), intent(in)                 :: z_r(:)  ! 
      real(ok), intent(in)                 :: c_r(:)  ! 
      !
      integer(ik)              :: ic_l, p1_l, p2_l, p_l, int_l
      integer(ik)              :: ic_r, p1_r, p2_r, p_r, int_r
      complex(ok), allocatable :: vb_r(:,:,:)  ! Buffer for recursions
      real(ok)                 :: wgt
      integer(ik)              :: n_rec        ! Number of additional recursion intermediated
      integer(ik)              :: alloc
      ! real(ok)                 :: ang_c_kind(lbound(ang_c,dim=1):ubound(ang_c,dim=1))
      real(ok)                 :: ang_c_kind(0:gam_nxyz-1)
      !
      !  Prepare angular factors of the right kind
      !
      ang_c_kind = real(ang_c,kind=ok)
      !
      p1_l = ang_loc(l_l)
      p2_l = ang_loc(l_l+1)-1
      p1_r = ang_loc(l_r)
      p2_r = ang_loc(l_r+1)-1
      !
      select case (what%op_name)
        case default
          !
          !  Assume general 3-centre integral; each unit of angular momentum requires
          !  an extra order in recursion.
          !
          n_rec = l_l + l_r
        case ('PLANEWAVE')
          n_rec = 0
        case ('R PLANEWAVE')
          n_rec = 1
      end select
      allocate (vb_r(0:p2_l,0:p2_r,0:n_rec),stat=alloc)
      if (alloc/=0) then
        write (out,"('import_gamess%os_1e_contraction: allocation failed. Parameters are:',3(1x,i8))") p2_l, p2_r, n_rec
        stop 'import_gamess%os_1e_contraction: Buffer allocation failed'
      end if
      !
      vb = 0._ok
      r_contractions_loop: do ic_r=1,size(z_r)
        l_contractions_loop: do ic_l=1,size(z_l)
          !
          !  The only part of the 1e integral code which is property-specific
          !  is the primitive-integral call below. Everything else is identical
          !  for all properties.
          !
          select case (what%op_name)
            case default
              write (out,"('import_gamess%os_1e_contraction_complex: Property ',a,' is not implemented')") trim(what%op_name)
              stop 'import_gamess%os_1e_contraction_complex: Unimplemented property'
            case ('PLANEWAVE')
              call os_1e_planewave(p2_l,p2_r,vb_r(:,:,0),r_l,z_l(ic_l),r_r,z_r(ic_r),real(what%kvec,kind=ok))
            case ('R PLANEWAVE')
              if (what%op_i<1 .or. what%op_i>3) then
                write (out,"('import_gamess%os_1e_contraction_complex: Dipole component ',i0,' is not between 1 and 3 in ',a)") &
                       what%op_i, trim(what%op_name)
                stop 'import_gamess%os_1e_contraction_complex: Bad dipole'
              end if
              call os_1e_planewave(p2_l,p2_r,vb_r(:,:,1),r_l,z_l(ic_l),r_r,z_r(ic_r),real(what%kvec,ok))
              call os_1e_pwdipole (what%op_i,vb_r(:,:,1), &
                                   p2_l,p2_r,vb_r(:,:,0),r_l,z_l(ic_l),r_r,z_r(ic_r),real(what%kvec,ok))
          end select
          !
          wgt = c_l(ic_l) * c_r(ic_r)
          r_accumulate: do p_r=p1_r,p2_r
            int_r = p_r - p1_r + 1
            l_accumulate: do p_l=p1_l,p2_l
              int_l = p_l - p1_l + 1
              vb(int_l,int_r) = vb(int_l,int_r) + wgt*ang_c_kind(p_l)*ang_c_kind(p_r)*vb_r(p_l,p_r,0)
            end do l_accumulate
          end do r_accumulate
          !
        end do l_contractions_loop
      end do r_contractions_loop
      !
      deallocate (vb_r,stat=alloc)
      if (alloc/=0) then
        stop 'import_gamess%os_1e_contraction_complex: Buffer deallocation failed'
      end if
!   end subroutine os_1e_contraction_complex
