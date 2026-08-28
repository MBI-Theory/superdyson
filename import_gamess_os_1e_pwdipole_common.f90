!
!  Note that we are using a simpler, direct expression in terms of overlap integrals. 
!  Effectively, we first recurse over the components of the moment operator, then
!  over the function. Our real dipole expressions are instead in terms of the much 
!  more complex O&S recursion, which recurses over the function first.
!
!   subroutine os_1e_pwdipole(id,vb_s,p2_l,p2_r,vb_p,r_l,z_l,r_r,z_r,kvec)
!     integer, parameter :: ok = rk
      integer(ik), intent(in)  :: id                  ! Dipole component to evaluate
      integer(ik), intent(in)  :: p2_l, p2_r          ! Upper dimensions of the vb_p array
      complex(ok), intent(in)  :: vb_s(0:p2_l,0:p2_r) ! Buffer containing primitive overlaps
      complex(ok), intent(out) :: vb_p(0:p2_l,0:p2_r) ! Final result.
      real(ok), intent(in)     :: r_l(:)              ! Centre of the left b.f.
      real(ok), intent(in)     :: z_l                 ! Orbital exponent of the left b.f.
      real(ok), intent(in)     :: r_r(:)              ! ditto, for the right b.f.
      real(ok), intent(in)     :: z_r                 ! 
      real(ok), intent(in)     :: kvec(3)             ! Planewave momentum
      !
      real(ok)    :: zeta, p(3)
      complex(ok) :: s, scl00
      integer(ik) :: m_l, m_r, id1, id2
      !
      call os_common_primitives(r_l,z_l,r_r,z_r,zeta=zeta,p=p)
      scl00 = p(id) + (0,1)*kvec(id)/(2*zeta)
      !
      right_loop: do m_r=0,p2_r
        left_loop: do m_l=0,p2_l
          s   = scl00 * vb_s(m_l,m_r)
          id1 = drop_xyz(m_l,id)
          id2 = drop_xyz(m_r,id)
          if (id1>=0) s = s + vb_s(id1,m_r) * ang_nxyz(m_l,id)/(2*zeta)
          if (id2>=0) s = s + vb_s(m_l,id2) * ang_nxyz(m_r,id)/(2*zeta)
          vb_p(m_l,m_r) = s
        end do left_loop
      end do right_loop
!   end subroutine os_1e_pwdipole
