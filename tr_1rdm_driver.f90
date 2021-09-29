program rdm
  use sd_core
  use tr_1rdm

  call read_core_input
  task = '1-rdm'
  call initialize_core_data
  call call_tr_1rdm
end program rdm
