!
!  Mathematical and physical constants not related to numerical precision
!
!  2018.10.26 - Updated to the recommended CODATA 2014 values
!
module constants
  use accuracy
  implicit none
  public
  !
  character(len=clen), save :: rcsid_constants = "$Id: constants.f90,v 1.2 2021/09/29 13:46:03 ps Exp ps $"
  !
  !  Mathematical constants
  !
  real(rk),  parameter  :: pi      = 3.1415926535897932384626433832795028841972_rk  ! pi, of course
  real(xk),  parameter  :: pi_xk   = 3.1415926535897932384626433832795028841972_xk  ! pi, of course
  real(xrk), parameter  :: pi_xrk  = 3.1415926535897932384626433832795028841972_xrk ! pi, of course
  real(rk),  parameter  :: twopi   = 2._rk * pi
  real(rk),  parameter  :: fourpi  = 4._rk * pi
  real(rk),  parameter  :: sqrtpi  = sqrt(pi)
  real(rk),  parameter  :: sqrt2pi = sqrt(twopi)
  !
  !  Powers of imaginary unity
  !
  !  I**N
  complex(rk), parameter ::  ipow(0:3) = (/ (1._rk,0._rk), (0._rk,1._rk), (-1._rk,0._rk), (0._rk,-1._rk) /)
  !  (-I)**N
  complex(rk), parameter :: mipow(0:3) = (/ (1._rk,0._rk), (0._rk,-1._rk), (-1._rk,0._rk), (0._rk,1._rk) /)
  !
  !  Physical constants. Old values are in the commented-out section
  !
! real(rk), parameter :: femtosecond = 41.341373336561364_rk           ! Conversion factor: from femtoseconds to au[t]
! real(rk), parameter :: abohr       = 0.5291772083_rk                 ! Conversion factor: from au[l] (Bohrs) to Angstrom
! real(rk), parameter :: h2ev        = 27.2113845_rk                   ! Conversion factor: from Hartree to electron-volts
! real(rk), parameter :: h2cm        = 219474.6313705_rk               ! Conversion factor: from Hartree to cm^-1 (wavenumbers)
! real(rk), parameter :: k_Boltzmann = 1.0_rk/315774.65_rk             ! Conversion factor: from Kelvin to Hartree 
!                                                                      !                    (aka Boltzmann constant)
! real(rk), parameter :: Hartree     = 4.35974417e-18_rk               ! Conversion factor: from au[e] (Hartree) to Joules
! real(rk), parameter :: vlight      = 137.0359991_rk                  ! Speed of light in atomic units; also the
! real(rk), parameter :: R_gas       = 8.314472_rk                     ! Molar gas constant, in J/mol-K
! real(rk), parameter :: N_Avogadro  = 6.0221415e23_rk                 ! Avogadro number, in particles/mole
! real(rk), parameter :: au2wcm2     = 3.5e16_rk                       ! Conversion factor: from square of peak intensity 
!                                                                      ! (linear polarization) to cycle-average intensity in W/cm^2
! real(rk), parameter :: unified_atomic_mass_unit = 1.660538782e-27_rk ! Unified atomic mass unit, in kg
! real(rk), parameter :: electron_mass_in_amu     = 5.4857990943e-4_rk ! Electron mass in u-amu
! real(rk), parameter :: au2wcm2     = 3.5e16_rk                       ! Conversion factor: from square of peak intensity 
!                                                                      ! (linear polarization) to cycle-average intensity in W/cm^2
! real(rk), parameter :: au2tesla = 2.350517382e5_rk                   ! Conversion factor from "zeeman" to SI Tesla
! real(rk), parameter :: gfactor  = 2.0023193043622_rk                 ! Electron g-factor
  !
  !  CODATA 2014 values
  !
  real(rk), parameter :: aut         = 24.18884326509_rk               ! Conversion factor: atom au[t] to attoseconds
  real(rk), parameter :: femtosecond = 1000._rk / aut                  ! Conversion factor: from femtoseconds to au[t]
  real(rk), parameter :: abohr       = 0.52917721067_rk                ! Conversion factor: from au[l] (Bohrs) to Angstrom
!
  real(rk), parameter :: h2ev        = 27.21138602_rk                  ! Conversion factor: from Hartree to electron-volts
  real(rk), parameter :: h2cm        = 219474.6313702_rk               ! Conversion factor: from Hartree to cm^-1 (wavenumbers)
  real(rk), parameter :: k_Boltzmann = 1.0_rk/315775.13_rk             ! Conversion factor: from Kelvin to Hartree 
                                                                       !                    (aka Boltzmann constant)
  real(rk), parameter :: Hartree     = 4.359744650e-18_rk              ! Conversion factor: from au[e] (Hartree) to Joules
  real(rk), parameter :: bar         = 101325._rk                      ! Conversion factor: from bars to Pascal 
                                                                       !                    (aka standard pressure)
  real(rk), parameter :: fine_alpha  = 7.2973525664e-3_rk              ! Fine-structure constant
  real(rk), parameter :: vlight      = 1._rk / fine_alpha              ! Speed of light in atomic units; also the
                                                                       ! inverse of the fine structure constant
  real(rk), parameter :: R_gas       = 8.3144598_rk                    ! Molar gas constant, in J/mol-K
  real(rk), parameter :: N_Avogadro  = 6.022140857e23_rk               ! Avogadro number, in particles/mole
  real(rk), parameter :: unified_atomic_mass_unit = 1.660539040e-27_rk ! Unified atomic mass unit, in kg
  real(rk), parameter :: electron_mass_in_amu     = 5.48579909070e-4_rk! Electron mass in u-amu
  real(rk), parameter :: au2wcm2     = (vlight/(8._rk*pi)) &           ! Conversion factor: from square of peak intensity 
                                     * Hartree &                       ! (linear polarization) to cycle-average intensity in W/cm^2
                                     * (1._rk/((abohr*1e-8_rk)**2)) &  ! The result is \approx 3.50944546 e+16
                                     * (1._rk/(aut*1e-18_rk))
  real(rk), parameter :: au2tesla    = 2.350517550e5_rk                ! Conversion factor from "zeeman" to SI Tesla
  real(rk), parameter :: gfactor     = 2.00231930436182_rk             ! Electron g-factor
  !
  !  Units
  !
  real(rk), parameter :: electron_mass   =  1._rk                      ! Mass of the electron in our system of units
  real(rk), parameter :: electron_charge = -1._rk                      ! Charhe of the electron in our system of units
end module constants
