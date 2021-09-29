module block_diag
!
!  Simplistic block-diagonalization routines
!
!  Note that this set of routines is not complete: only the data types 
!  We encountered so far are implemented. Sorry.
!
!  2013 Nov 02 - Added block_geev routine
!  2013 Jul 16 - Added sorting of eigenvalues / singular values
!  2014 Jun ?? - Added quad-precision
!
  use accuracy
  use timer
  use lapack
  use printing
  use sort_tools
  implicit none
  private
  public block_syev, block_gesvd, block_geev
  public rcsid_block_diag
  !
  character(len=clen), save :: rcsid_block_diag = "$Id: block_diag.f90,v 1.8 2021/09/29 13:43:22 ps Exp $"
  !
  interface block_syev
     module procedure block_syev_real
!*qd module procedure block_syev_quad
  end interface block_syev
  !
  interface block_geev
     module procedure block_zgeev_real
!*qd module procedure block_zgeev_quad
  end interface block_geev
  !
  interface block_gesvd
     module procedure block_gesvd_real
!*qd module procedure block_gesvd_quad
     module procedure block_gesvd_complex
!*qd module procedure block_gesvd_quad_complex
  end interface block_gesvd
  !
  integer(ik), parameter :: verbose = 0

  contains

  subroutine block_syev_real(h,e,eps)
    integer, parameter :: ok = rk
    include 'block_diag_syev_common.f90'
  end subroutine block_syev_real

  subroutine block_syev_quad(h,e,eps)
    integer, parameter :: ok = xrk
    include 'block_diag_syev_common.f90'
  end subroutine block_syev_quad
  !
  !  This routine is modeled on lapack_zgeev2
  !
  subroutine block_zgeev_real(h,e,eps)
    integer, parameter :: ok = rk
    include "block_diag_zgeev_common.f90"
  end subroutine block_zgeev_real

  subroutine block_zgeev_quad(h,e,eps)
    integer, parameter :: ok = xrk
    include "block_diag_zgeev_common.f90"
  end subroutine block_zgeev_quad

  subroutine block_gesvd_real(a,s,u,vth,eps)
    integer, parameter :: ok = rk
    ! Only the variables where the basic type (not the kind) depend on the 
    ! specific routine are listed below; the remaining variables are defined
    ! in the include file
    real(ok), intent(inout) :: a  (:,:)       ! In:  Matrix to be decomposed
                                              ! Out: Matrix destroyed
    real(ok), intent(out)   :: u  (:,:)       ! Out: Left singular vectors
    real(ok), intent(out)   :: vth(:,:)       ! Out: Right singular vectors, transposed & conjugated
                                              ! The overall result is A = U S VTH
    !
    real(ok), allocatable   :: tmp_a(:,:)     ! Temporary used for decomposing the connected block
    real(ok), allocatable   :: tmp_u(:,:)     ! ditto
    real(ok), allocatable   :: tmp_vth(:,:)   ! ditto
    !
    include 'block_diag_gesvd_common.f90'
  end subroutine block_gesvd_real

  subroutine block_gesvd_quad(a,s,u,vth,eps)
    integer, parameter :: ok = xrk
    ! Only the variables where the basic type (not the kind) depend on the 
    ! specific routine are listed below; the remaining variables are defined
    ! in the include file
    real(ok), intent(inout) :: a  (:,:)       ! In:  Matrix to be decomposed
                                              ! Out: Matrix destroyed
    real(ok), intent(out)   :: u  (:,:)       ! Out: Left singular vectors
    real(ok), intent(out)   :: vth(:,:)       ! Out: Right singular vectors, transposed & conjugated
                                              ! The overall result is A = U S VTH
    !
    real(ok), allocatable   :: tmp_a(:,:)     ! Temporary used for decomposing the connected block
    real(ok), allocatable   :: tmp_u(:,:)     ! ditto
    real(ok), allocatable   :: tmp_vth(:,:)   ! ditto
    !
    include 'block_diag_gesvd_common.f90'
  end subroutine block_gesvd_quad

  subroutine block_gesvd_complex(a,s,u,vth,eps)
    integer, parameter :: ok = rk
    ! Only the variables where the basic type (not the kind) depend on the 
    ! specific routine are listed below; the remaining variables are defined
    ! in the include file
    complex(ok), intent(inout) :: a  (:,:)       ! In:  Matrix to be decomposed
                                                 ! Out: Matrix destroyed
    complex(ok), intent(out)   :: u  (:,:)       ! Out: Left singular vectors
    complex(ok), intent(out)   :: vth(:,:)       ! Out: Right singular vectors, transposed & conjugated
                                                 ! The overall result is A = U S VTH
    !
    complex(ok), allocatable   :: tmp_a(:,:)     ! Temporary used for decomposing the connected block
    complex(ok), allocatable   :: tmp_u(:,:)     ! ditto
    complex(ok), allocatable   :: tmp_vth(:,:)   ! ditto
    !
    include 'block_diag_gesvd_common.f90'
  end subroutine block_gesvd_complex

  subroutine block_gesvd_quad_complex(a,s,u,vth,eps)
    integer, parameter :: ok = xrk
    ! Only the variables where the basic type (not the kind) depend on the 
    ! specific routine are listed below; the remaining variables are defined
    ! in the include file
    complex(ok), intent(inout) :: a  (:,:)       ! In:  Matrix to be decomposed
                                                 ! Out: Matrix destroyed
    complex(ok), intent(out)   :: u  (:,:)       ! Out: Left singular vectors
    complex(ok), intent(out)   :: vth(:,:)       ! Out: Right singular vectors, transposed & conjugated
                                                 ! The overall result is A = U S VTH
    !
    complex(ok), allocatable   :: tmp_a(:,:)     ! Temporary used for decomposing the connected block
    complex(ok), allocatable   :: tmp_u(:,:)     ! ditto
    complex(ok), allocatable   :: tmp_vth(:,:)   ! ditto
    !
    include 'block_diag_gesvd_common.f90'
  end subroutine block_gesvd_quad_complex

end module block_diag
