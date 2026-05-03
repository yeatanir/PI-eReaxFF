program test_crystx_periodic_loader
  use ereaxff_kinds, only: dp
  use ereaxff_system_state, only: system_state, cell_is_defined, cell_lengths_A
  use ereaxff_legacy_geo_loader, only: load_legacy_geo_file
  use ereaxff_reaxff_compat_output, only: write_reaxff_compat_snapshot
  implicit none

  type(system_state) :: state
  real(dp), allocatable :: charges(:)
  character(len=512) :: geo_path, out_dir
  integer :: nparticles

  call get_command_argument(1, geo_path)
  call get_command_argument(2, out_dir)
  if (len_trim(geo_path) == 0 .or. len_trim(out_dir) == 0) error stop &
      "test_crystx_periodic_loader: missing arguments"

  call load_legacy_geo_file(trim(geo_path), state)
  if (.not. cell_is_defined(state)) error stop "test_crystx_periodic_loader: cell not defined"
  nparticles = state%nparticles
  allocate(charges(nparticles))
  charges = 0.0_dp

  call execute_command_line('mkdir -p "' // trim(out_dir) // '"')
  call write_reaxff_compat_snapshot(trim(out_dir), 1, state%position_A, charges, 0.0_dp, 0.0_dp, 0.0_dp, 1, state)

  write(*,'(a,l1)') "cell_defined=", cell_is_defined(state)
  write(*,'(a,l1)') "cell_periodic=", state%cell%periodic
  write(*,'(a,3(1x,f12.5))') "cell_lengths_A:", cell_lengths_A(state)
  write(*,'(a,i0)') "nparticles=", nparticles
end program test_crystx_periodic_loader
