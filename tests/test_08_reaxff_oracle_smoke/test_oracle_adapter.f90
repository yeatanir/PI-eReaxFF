program test_oracle_adapter
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use ereaxff_kinds, only: dp
  use ereaxff_system_state, only: system_state, init_system_state, set_particle
  use ereaxff_oracle_adapter, only: oracle_adapter_state, oracle_adapter_result, &
      oracle_adapter_init, oracle_adapter_eval
  implicit none

  integer, parameter :: ndim = 3
  integer, parameter :: nparticles = 61
  real(dp), parameter :: kcal2eV = 0.043364_dp
  real(dp), parameter :: tol_energy = 1.0e-10_dp
  real(dp), parameter :: tol_force = 1.0e-12_dp
  real(dp), parameter :: tol_charge = 1.0e-12_dp
  real(dp), parameter :: tol_force_entry = 1.0e-12_dp
  type(system_state) :: template
  type(oracle_adapter_state) :: adapter
  type(oracle_adapter_result) :: result_first, result_second
  real(dp) :: positions(ndim, nparticles)
  real(dp) :: charges(nparticles)
  real(dp) :: solver_guess(nparticles)
  character(len=512) :: v1_root, v2_root, work_root
  character(len=32) :: mode
  integer :: arg_count

  call get_environment_variable("EREAXFF_V1_ROOT", v1_root)
  call get_environment_variable("EREAXFF_V2_ROOT", v2_root)
  call get_environment_variable("EREAXFF_ADAPTER_WORK_ROOT", work_root)
  if (len_trim(v1_root) == 0) v1_root = "/Users/yeatanir/PhD/PIMD/ereaxff_pimd/V1"
  if (len_trim(v2_root) == 0) v2_root = "/Users/yeatanir/PhD/PIMD/ereaxff_pimd/V2"
  if (len_trim(work_root) == 0) work_root = "build/test_adapter_work"

  arg_count = command_argument_count()
  mode = "smoke"
  if (arg_count >= 1) call get_command_argument(1, mode)

  call build_phase4_template(template, positions)
  call oracle_adapter_init(adapter, template, trim(v1_root), trim(v2_root), trim(work_root))

  charges = 0.0_dp
  solver_guess = 0.0_dp

  if (trim(mode) == "nan") then
    call exercise_nan_rejection(adapter, positions, charges, solver_guess)
  else
    call oracle_adapter_eval(adapter, positions, charges, solver_guess, bead_id=1, &
        charge_policy="bead_local_acks2", result=result_first)
    call oracle_adapter_eval(adapter, positions, charges, solver_guess, bead_id=1, &
        charge_policy="bead_local_acks2", result=result_second)

    call assert_true("legacy status zero", result_first%force%status == 0)
    call assert_true("repeat status zero", result_second%force%status == 0)
    call assert_close("energy reproducible", result_first%force%energy_eV, &
        result_second%force%energy_eV, tol_energy)
    call assert_close("force reproducible", maxval(abs(result_first%force%forces_eV_per_A - &
        result_second%force%forces_eV_per_A)), 0.0_dp, tol_force)
    call assert_close("charge reproducible", maxval(abs(result_first%force%charges - &
        result_second%force%charges)), 0.0_dp, tol_charge)
    call assert_close("solver reproducible", maxval(abs(result_first%solver_guess - &
        result_second%solver_guess)), 0.0_dp, tol_charge)

    call assert_close("energy matches legacy kcal->eV", result_first%force%energy_eV, &
        -10916.85308542530_dp * kcal2eV, 1.0e-10_dp)
    call assert_close("max force matches legacy", maxval(abs(result_first%force%forces_eV_per_A)), &
        9.903670910000796e-02_dp, tol_force_entry)
    call assert_close("charge sum matches legacy", sum(result_first%force%charges), &
        -5.717690199427694e-08_dp, tol_charge)
    call assert_close("atom 1 charge matches legacy", result_first%force%charges(1), &
        -1.000000259853050_dp, tol_charge)
    call assert_close("atom 1 solver matches legacy", result_first%solver_guess(1), &
        -1.000000259853050_dp, tol_charge)

    call write_summary(result_first)
  end if

contains

  subroutine build_phase4_template(state, positions_out)
    type(system_state), intent(inout) :: state
    real(dp), intent(out) :: positions_out(ndim, nparticles)
    integer :: ix, iy, unit_index, idx
    real(dp), parameter :: a_cc = 1.4210_dp
    real(dp), parameter :: z_c = 10.0_dp
    real(dp), parameter :: z_el = 10.001_dp
    real(dp), parameter :: el_x = 7.38300_dp
    real(dp), parameter :: el_y = 5.68433_dp
    real(dp), parameter :: el_z = z_el
    real(dp) :: ax, ay, x, y
    real(dp), dimension(4, 2) :: unit

    unit = reshape([0.0_dp, 0.0_dp, 0.5_dp, 1.0_dp/6.0_dp, 0.5_dp, 3.0_dp/6.0_dp, &
        0.0_dp, 4.0_dp/6.0_dp], shape(unit))
    ax = sqrt(3.0_dp) * a_cc
    ay = 3.0_dp * a_cc

    call init_system_state(state, ndim, nparticles)
    call set_particle(state, 1, 1, "El", 0.00054858_dp, [el_x, el_y, el_z])
    positions_out(:, 1) = [el_x, el_y, el_z]

    idx = 2
    do ix = 0, 4
      do iy = 0, 2
        do unit_index = 1, 4
          x = (real(ix, dp) + unit(unit_index, 1)) * ax
          y = (real(iy, dp) + unit(unit_index, 2)) * ay
          call set_particle(state, idx, idx, "C", 12.011_dp, [x, y, z_c])
          positions_out(:, idx) = [x, y, z_c]
          idx = idx + 1
        end do
      end do
    end do
  end subroutine build_phase4_template

  subroutine exercise_nan_rejection(adapter_in, positions_in, charges_in, solver_in)
    type(oracle_adapter_state), intent(inout) :: adapter_in
    real(dp), intent(in) :: positions_in(ndim, nparticles)
    real(dp), intent(in) :: charges_in(nparticles)
    real(dp), intent(in) :: solver_in(nparticles)
    type(oracle_adapter_result) :: result_nan
    real(dp) :: bad_positions(ndim, nparticles)

    bad_positions = positions_in
    bad_positions(1, 1) = ieee_value(0.0_dp, ieee_quiet_nan)
    call oracle_adapter_eval(adapter_in, bad_positions, charges_in, solver_in, bead_id=1, &
        charge_policy="bead_local_acks2", result=result_nan)
  end subroutine exercise_nan_rejection

  subroutine write_summary(result_in)
    type(oracle_adapter_result), intent(in) :: result_in

    write(*,'(a)') "test_08_reaxff_oracle_smoke: PASS"
    write(*,'(a,i0)') "  natom = ", size(result_in%force%charges)
    write(*,'(a,es24.15)') "  energy_eV = ", result_in%force%energy_eV
    write(*,'(a,es24.15)') "  max_force = ", maxval(abs(result_in%force%forces_eV_per_A))
    write(*,'(a,es24.15)') "  charge_sum = ", sum(result_in%force%charges)
    write(*,'(a,es24.15)') "  solver_atom1 = ", result_in%solver_guess(1)
  end subroutine write_summary

  subroutine assert_true(label, condition)
    character(len=*), intent(in) :: label
    logical, intent(in) :: condition

    if (.not. condition) then
      write(*,'(a)') "FAIL: " // trim(label)
      error stop "oracle adapter regression failed"
    end if
  end subroutine assert_true

  subroutine assert_close(label, value, expected, threshold)
    character(len=*), intent(in) :: label
    real(dp), intent(in) :: value, expected, threshold

    if (abs(value - expected) > threshold) then
      write(*,'(a)') "FAIL: " // trim(label)
      write(*,'(a,es24.16)') "  value    = ", value
      write(*,'(a,es24.16)') "  expected = ", expected
      write(*,'(a,es24.16)') "  tol      = ", threshold
      error stop "oracle adapter regression failed"
    end if
  end subroutine assert_close

end program test_oracle_adapter
