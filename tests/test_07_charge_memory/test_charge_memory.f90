program test_charge_memory
  use ereaxff_kinds, only: dp
  use ereaxff_charge_coupling_policy, only: charge_memory_state, charge_policy_is_valid, &
      charge_policy_uses_bead_local_memory, init_charge_memory, set_bead_charges, &
      charge_sum_for_bead, max_charge_difference_across_beads, &
      charge_policy_bead_local_acks2, charge_policy_centroid_acks2, &
      charge_policy_frozen_charges, charge_policy_no_charge_coupling
  implicit none

  integer, parameter :: nparticles = 3
  integer, parameter :: p = 4
  type(charge_memory_state) :: memory
  real(dp) :: charges(nparticles)
  integer :: bead

  call assert_true("bead_local policy valid", charge_policy_is_valid(charge_policy_bead_local_acks2))
  call assert_true("centroid policy valid", charge_policy_is_valid(charge_policy_centroid_acks2))
  call assert_true("frozen policy valid", charge_policy_is_valid(charge_policy_frozen_charges))
  call assert_true("no charge policy valid", charge_policy_is_valid(charge_policy_no_charge_coupling))
  call assert_true("invalid policy rejected", .not. charge_policy_is_valid("shared_accidental_memory"))
  call assert_true("bead local uses bead memory", charge_policy_uses_bead_local_memory(charge_policy_bead_local_acks2))
  call assert_true("centroid does not use bead-local memory", &
      .not. charge_policy_uses_bead_local_memory(charge_policy_centroid_acks2))

  call init_charge_memory(memory, nparticles, p, charge_policy_bead_local_acks2)
  do bead = 1, p
    charges = [0.1_dp * real(bead, dp), -0.2_dp, 0.05_dp * real(bead * bead, dp)]
    call set_bead_charges(memory, bead, charges)
  end do

  call assert_close("charge sum bead 1", charge_sum_for_bead(memory, 1), -0.05_dp, 1.0e-14_dp)
  call assert_close("charge sum bead 4", charge_sum_for_bead(memory, 4), 1.0_dp, 1.0e-14_dp)
  call assert_close("max charge difference across beads", max_charge_difference_across_beads(memory), &
      0.75_dp, 1.0e-14_dp)
  call assert_close("bead 1 memory remains independent", memory%charges_by_bead(1, 1), 0.1_dp, 1.0e-14_dp)
  call assert_close("bead 4 memory remains independent", memory%charges_by_bead(1, 4), 0.4_dp, 1.0e-14_dp)

  write(*,'(a)') "test_07_charge_memory: PASS"
  write(*,'(a,es12.4)') "  charge_sum_bead_1 = ", charge_sum_for_bead(memory, 1)
  write(*,'(a,es12.4)') "  charge_sum_bead_4 = ", charge_sum_for_bead(memory, 4)
  write(*,'(a,es12.4)') "  max_charge_diff   = ", max_charge_difference_across_beads(memory)

contains

  subroutine assert_true(label, condition)
    character(len=*), intent(in) :: label
    logical, intent(in) :: condition

    if (.not. condition) then
      write(*,'(a)') "FAIL: " // trim(label)
      error stop "charge policy regression failed"
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
      error stop "charge memory regression failed"
    end if
  end subroutine assert_close

end program test_charge_memory

