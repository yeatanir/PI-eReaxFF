program test_oracle_charge_policy_modes
  use ereaxff_kinds, only: dp
  use ereaxff_system_state, only: system_state, init_system_state, set_particle
  use ereaxff_charge_coupling_policy, only: charge_memory_state, init_charge_memory, &
      charge_sum_for_bead, max_charge_difference_across_beads, &
      charge_policy_bead_local_acks2, charge_policy_centroid_acks2, &
      charge_policy_frozen_charges, charge_policy_no_charge_coupling, &
      charge_policy_supports_oracle_manager, charge_policy_uses_centroid_contraction
  use ereaxff_oracle_adapter, only: oracle_adapter_state, oracle_adapter_init
  use ereaxff_bead_force_manager, only: oracle_bead_summary, evaluate_oracle_beads
  implicit none

  integer, parameter :: ndim = 3
  integer, parameter :: nparticles = 61
  integer, parameter :: p = 2
  real(dp), parameter :: kcal2eV = 0.043364_dp
  type(system_state) :: template
  type(oracle_adapter_state) :: adapter
  type(charge_memory_state) :: local_memory, centroid_memory
  type(oracle_bead_summary) :: local_summary, centroid_summary
  real(dp) :: positions(ndim, nparticles, p)
  real(dp) :: local_solver(nparticles, p)
  real(dp) :: centroid_solver(nparticles, p)
  character(len=512) :: v1_root, v2_root, work_root
  integer :: bead

  call get_environment_variable("EREAXFF_V1_ROOT", v1_root)
  call get_environment_variable("EREAXFF_V2_ROOT", v2_root)
  call get_environment_variable("EREAXFF_ADAPTER_WORK_ROOT", work_root)
  if (len_trim(v1_root) == 0) v1_root = "/Users/yeatanir/PhD/PIMD/ereaxff_pimd/V1"
  if (len_trim(v2_root) == 0) v2_root = "/Users/yeatanir/PhD/PIMD/ereaxff_pimd/V2"
  if (len_trim(work_root) == 0) work_root = "build/test_adapter_work"

  call build_phase4_template(template, positions(:, :, 1))
  positions(:, :, 2) = positions(:, :, 1)
  call oracle_adapter_init(adapter, template, trim(v1_root), trim(v2_root), trim(work_root))

  call init_charge_memory(local_memory, nparticles, p, charge_policy_bead_local_acks2)
  call init_charge_memory(centroid_memory, nparticles, p, charge_policy_centroid_acks2)
  local_solver = 0.0_dp
  centroid_solver = 0.0_dp

  call evaluate_oracle_beads(adapter, positions, local_memory, local_solver, &
      charge_policy_bead_local_acks2, local_summary)
  call evaluate_oracle_beads(adapter, positions, centroid_memory, centroid_solver, &
      charge_policy_centroid_acks2, centroid_summary)

  call assert_true("bead-local policy supports oracle manager", &
      charge_policy_supports_oracle_manager(charge_policy_bead_local_acks2))
  call assert_true("centroid policy supports oracle manager", &
      charge_policy_supports_oracle_manager(charge_policy_centroid_acks2))
  call assert_true("centroid policy uses centroid contraction", &
      charge_policy_uses_centroid_contraction(charge_policy_centroid_acks2))
  call assert_true("frozen policy does not support oracle manager", &
      .not. charge_policy_supports_oracle_manager(charge_policy_frozen_charges))
  call assert_true("no coupling policy does not support oracle manager", &
      .not. charge_policy_supports_oracle_manager(charge_policy_no_charge_coupling))

  call assert_close("local and centroid energy agree", local_summary%pimd_potential_eV, &
      centroid_summary%pimd_potential_eV, 1.0e-12_dp)
  call assert_close("local and centroid energy match legacy", local_summary%pimd_potential_eV, &
      -10916.85308542530_dp * kcal2eV, 1.0e-10_dp)
  call assert_close("local and centroid force fields agree", &
      maxval(abs(local_summary%scaled_forces_eV_per_A - centroid_summary%scaled_forces_eV_per_A)), &
      0.0_dp, 1.0e-12_dp)
  call assert_close("local and centroid charges agree", &
      maxval(abs(local_summary%bead_charges - centroid_summary%bead_charges)), &
      0.0_dp, 1.0e-12_dp)
  call assert_close("local and centroid solver guesses agree", &
      maxval(abs(local_summary%bead_solver_guess - centroid_summary%bead_solver_guess)), &
      0.0_dp, 1.0e-12_dp)
  call assert_close("local charge memory remains bead-local", &
      max_charge_difference_across_beads(local_memory), 0.0_dp, 1.0e-12_dp)
  call assert_close("centroid charge memory is shared", &
      max_charge_difference_across_beads(centroid_memory), 0.0_dp, 1.0e-12_dp)

  do bead = 1, p
    call assert_close("local charge sum", charge_sum_for_bead(local_memory, bead), &
        charge_sum_for_bead(centroid_memory, bead), 1.0e-12_dp)
  end do

  write(*,'(a)') "test_13_oracle_charge_policy_modes: PASS"
  write(*,'(a,es24.15)') "  energy_eV = ", local_summary%pimd_potential_eV
  write(*,'(a,es24.15)') "  max_force = ", maxval(abs(local_summary%scaled_forces_eV_per_A))
  write(*,'(a,es24.15)') "  charge_sum = ", charge_sum_for_bead(local_memory, 1)
  write(*,'(a,es24.15)') "  centroid_charge_sum = ", charge_sum_for_bead(centroid_memory, 1)

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
      error stop "oracle charge-policy regression failed"
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
      error stop "oracle charge-policy regression failed"
    end if
  end subroutine assert_close

end program test_oracle_charge_policy_modes
