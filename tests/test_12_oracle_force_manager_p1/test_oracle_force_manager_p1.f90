program test_oracle_force_manager_p1
  use ereaxff_kinds, only: dp
  use ereaxff_system_state, only: system_state, init_system_state, set_particle
  use ereaxff_charge_coupling_policy, only: charge_memory_state, init_charge_memory, &
      charge_sum_for_bead, charge_policy_bead_local_acks2
  use ereaxff_oracle_adapter, only: oracle_adapter_state, oracle_adapter_init
  use ereaxff_bead_force_manager, only: oracle_bead_summary, evaluate_oracle_beads
  implicit none

  integer, parameter :: ndim = 3
  integer, parameter :: nparticles = 61
  integer, parameter :: p = 1
  real(dp), parameter :: kcal2eV = 0.043364_dp
  type(system_state) :: template
  type(oracle_adapter_state) :: adapter
  type(charge_memory_state) :: charge_memory
  type(oracle_bead_summary) :: summary
  real(dp) :: positions(ndim, nparticles, p)
  real(dp) :: solver_guess(nparticles, p)
  character(len=512) :: v1_root, v2_root, work_root

  call get_environment_variable("EREAXFF_V1_ROOT", v1_root)
  call get_environment_variable("EREAXFF_V2_ROOT", v2_root)
  call get_environment_variable("EREAXFF_ADAPTER_WORK_ROOT", work_root)
  if (len_trim(v1_root) == 0) v1_root = "/Users/yeatanir/PhD/PIMD/ereaxff_pimd/V1"
  if (len_trim(v2_root) == 0) v2_root = "/Users/yeatanir/PhD/PIMD/ereaxff_pimd/V2"
  if (len_trim(work_root) == 0) work_root = "build/test_adapter_work"

  call build_phase4_template(template, positions(:, :, 1))
  call oracle_adapter_init(adapter, template, trim(v1_root), trim(v2_root), trim(work_root))
  call init_charge_memory(charge_memory, nparticles, p, charge_policy_bead_local_acks2)
  solver_guess = 0.0_dp

  call evaluate_oracle_beads(adapter, positions, charge_memory, solver_guess, &
      charge_policy_bead_local_acks2, summary)

  call assert_true("manager status zero", summary%bead_status(1) == 0)
  call assert_close("manager energy matches legacy", summary%pimd_potential_eV, &
      -10916.85308542530_dp * kcal2eV, 1.0e-10_dp)
  call assert_close("manager max force matches legacy", &
      maxval(abs(summary%scaled_forces_eV_per_A(:, :, 1))), 9.903670910000796e-02_dp, 1.0e-12_dp)
  call assert_close("manager charge sum matches legacy", sum(summary%bead_charges(:, 1)), &
      -5.717690199427694e-08_dp, 1.0e-12_dp)
  call assert_close("manager solver memory matches legacy", summary%bead_solver_guess(1, 1), &
      -1.000000259853050_dp, 1.0e-12_dp)
  call assert_close("charge memory updated", charge_sum_for_bead(charge_memory, 1), &
      -5.717690199427694e-08_dp, 1.0e-12_dp)

  write(*,'(a)') "test_12_oracle_force_manager_p1: PASS"
  write(*,'(a,es24.15)') "  energy_eV = ", summary%pimd_potential_eV
  write(*,'(a,es24.15)') "  max_force = ", maxval(abs(summary%scaled_forces_eV_per_A(:, :, 1)))
  write(*,'(a,es24.15)') "  charge_sum = ", sum(summary%bead_charges(:, 1))
  write(*,'(a,es24.15)') "  solver_atom1 = ", summary%bead_solver_guess(1, 1)

contains

  subroutine build_phase4_template(state, positions_out)
    type(system_state), intent(inout) :: state
    real(dp), intent(out) :: positions_out(ndim, nparticles)
    integer :: ix, iy, unit_index, idx
    real(dp), parameter :: a_cc = 1.4210_dp
    real(dp), parameter :: z_c = 10.0_dp
    real(dp), parameter :: el_x = 7.38300_dp
    real(dp), parameter :: el_y = 5.68433_dp
    real(dp), parameter :: el_z = 10.001_dp
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

  subroutine assert_true(label, condition)
    character(len=*), intent(in) :: label
    logical, intent(in) :: condition

    if (.not. condition) then
      write(*,'(a)') "FAIL: " // trim(label)
      error stop "oracle force-manager regression failed"
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
      error stop "oracle force-manager regression failed"
    end if
  end subroutine assert_close

end program test_oracle_force_manager_p1
