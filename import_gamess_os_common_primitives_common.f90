!
!  Calculate quantities which appear in (nearly) all Obara-Saika recurrences
!
!   subroutine os_common_primitives(r_l,z_l,r_r,z_r,xi,zeta,p,s00,s00a)
!     integer, parameter              :: ok = rk
      real(ok), intent(in)            :: r_l(:)  ! Centre of the left b.f.
      real(ok), intent(in)            :: z_l     ! Orbital exponent of the left b.f.
      real(ok), intent(in)            :: r_r(:)  ! ditto, for the right b.f.
      real(ok), intent(in)            :: z_r     ! 
      real(ok), intent(out), optional :: xi      ! Eq. 13 of Obara-Saika
      real(ok), intent(out), optional :: zeta    ! Eq. 14
      real(ok), intent(out), optional :: p(:)    ! Eq. 15
      real(ok), intent(out), optional :: s00     ! Eq. 22
      real(ok), intent(out), optional :: s00a    ! Ahlrichs' version of s00 - no prefactor
!     !
      if (present(xi)  ) xi   = z_l*z_r/(z_l+z_r)
      if (present(zeta)) zeta = z_l + z_r
      if (present(p)   ) p    = (z_l*r_l + z_r*r_r)/(z_l+z_r)
      if (present(s00) ) s00  = (real(pi,kind=ok)/(z_l+z_r))**1.5_ok &
                                                * exp(-(z_l*z_r/(z_l+z_r))*sum((r_l-r_r)**2))
      if (present(s00a)) s00a =                   exp(-(z_l*z_r/(z_l+z_r))*sum((r_l-r_r)**2))
!   end subroutine os_common_primitives
