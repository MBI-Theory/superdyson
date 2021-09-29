!
!  Report versions of all modules in this build
!
 subroutine versions
   use accuracy
   use block_determinant
   use block_diag
   use block_matrices
   use constants
   use gamess_internal
   use import_gamess
   use lapack
   use math
   use os_integral_operators
   use printing
   use sd_core
   use sort_tools
   use superdyson
   use timer
   use tr_1rdm
   !
   write (out,"(t5,a)") trim(rcsid_accuracy)
   write (out,"(t5,a)") trim(rcsid_block_determinant)
   write (out,"(t5,a)") trim(rcsid_block_diag)
   write (out,"(t5,a)") trim(rcsid_block_matrices)
   write (out,"(t5,a)") trim(rcsid_constants)
   write (out,"(t5,a)") trim(rcsid_gamess_internal)
   write (out,"(t5,a)") trim(rcsid_import_gamess)
   write (out,"(t5,a)") trim(rcsid_lapack)
   write (out,"(t5,a)") trim(rcsid_math)
   write (out,"(t5,a)") trim(rcsid_os_integral_operators)
   write (out,"(t5,a)") trim(rcsid_printing)
   write (out,"(t5,a)") trim(rcsid_sd_core)
   write (out,"(t5,a)") trim(rcsid_sort_tools)
   write (out,"(t5,a)") trim(rcsid_superdyson)
   write (out,"(t5,a)") trim(rcsid_timer)
   write (out,"(t5,a)") trim(rcsid_tr_1rdm)
 end subroutine versions
