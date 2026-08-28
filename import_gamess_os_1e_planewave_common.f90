!
!  Calculate primitive planewave integrals. This is a simple modification of the
!  os_1e_overlap routine
!
!   subroutine os_1e_planewave(p2_l,p2_r,vb_p,r_l,z_l,r_r,z_r,kvec)
!     integer, parameter       :: ok = rk
      integer(ik), intent(in)  :: p2_l, p2_r          ! Upper dimensions of the vb_p array
      complex(ok), intent(out) :: vb_p(0:p2_l,0:p2_r) ! Buffer used for recurrences, and the final result
      real(ok), intent(in)     :: r_l(:)              ! Centre of the left b.f.
      real(ok), intent(in)     :: z_l                 ! Orbital exponent of the left b.f.
      real(ok), intent(in)     :: r_r(:)              ! ditto, for the right b.f.
      real(ok), intent(in)     :: z_r                 ! 
      real(ok), intent(in)     :: kvec(3)             ! Planewave momentum
      !
      real(ok)    :: zeta, xi, p(3), s00
      complex(ok) :: s
      integer(ik) :: m_l, m_r, ic, id1, id2, id3
      !
      call os_common_primitives(r_l,z_l,r_r,z_r,xi=xi,zeta=zeta,p=p,s00=s00)
      !
      !  The primitive overlap has an extra term in the case of planevawe
      !
      vb_p(0,0) = s00 * exp(-sum(kvec**2)/(4*zeta)) * exp((0,1)*sum(kvec*p))
      !
      !  Bootstrap: calculate all overlaps on the right, keeping left at zero. 
      !
      right_bootstrap: do m_r=1,p2_r
        find_right: do ic=1,3
          id1 = drop_xyz(m_r,ic)
          if (id1<0) cycle find_right
          s = (p(ic)-r_r(ic)+(0,1)*kvec(ic)/(2.0_ok*zeta)) * vb_p(0,id1)
          id2 = drop_xyz(id1,ic)
          if (id2>=0) s = s + vb_p(0,id2)*ang_nxyz(id1,ic)/(2.0_ok*zeta)
          vb_p(0,m_r) = s
          exit find_right
        end do find_right
      end do right_bootstrap
      !
      !  Now the general case: always recurse on the left. 
      !
      right_loop: do m_r=0,p2_r
        left_recurrence: do m_l=1,p2_l
          find_left: do ic=1,3
            id1 = drop_xyz(m_l,ic)
            if (id1<0) cycle find_left
            s = (p(ic)-r_l(ic)+(0,1)*kvec(ic)/(2.0_ok*zeta)) * vb_p(id1,m_r)
            id2 = drop_xyz(id1,ic)
            if (id2>=0) s = s + vb_p(id2,m_r)*ang_nxyz(id1,ic)/(2.0_ok*zeta)
            id3 = drop_xyz(m_r,ic)
            if (id3>=0) s = s + vb_p(id1,id3)*ang_nxyz(m_r,ic)/(2.0_ok*zeta)
            vb_p(m_l,m_r) = s
            exit find_left
          end do find_left
        end do left_recurrence
      end do right_loop
!   end subroutine os_1e_planewave
