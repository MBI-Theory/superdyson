!
!  This is a more modern version of the input for superdyson/tr_1rdm
!  In addition to the namelist-based, extensible input, it also
!  offers support for greater choice of operators in the matrix 
!  elements.
!
!  The input is still handled in sd_core.f90; depending on the
!  input data, we'll chose either superdyson or tr_1rdm branch.
!
program sd_v2
  use sd_core
  use superdyson
  use tr_1rdm
  !
  external versions
  !
  write (out,"('Version: ',a/)") __BUILD_ID__
  call versions
  !
  call TimerStart('start')
  call read_input_new
  call initialize_core_data
  !
  select case (task)
    case default
      write (out,"(' Task ',a,' is not recognized')") trim(task)
      stop 'sd_v2 - bad task'
    case ('braket','BRAKET')
      call call_superdyson
    case ('1-rdm','1-RDM')
      call call_tr_1rdm
  end select
  call TimerStop('start')
  call TimerReport
end program sd_v2
