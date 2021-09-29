program dyson
  use sd_core
  use superdyson

  call read_core_input
  task = 'braket'
  call initialize_core_data
  call call_superdyson
end program dyson
