program harmonic_pimd_main
  use ereaxff_kinds, only: dp, i8
  use ereaxff_rng, only: rng_state, rng_seed
  use ereaxff_ring_polymer, only: ring_polymer_state, init_ring_polymer
  use ereaxff_normal_modes, only: nm_to_bead
  use ereaxff_spring_energy, only: compute_spring_energy_nm
  use ereaxff_harmonic_quantum_variance, only: harmonic_finite_p_bead_variance_A2
  use ereaxff_thermostat_base, only: kinetic_energy_2d_eV
  use ereaxff_thermostat_pile_l, only: pile_l_thermostat_state, pile_l_init
  use ereaxff_baoab_integrator, only: baoab_harmonic_pimd_step_1d
  use ereaxff_pimd_output, only: write_centroid_xyz_frame, compute_centroid_1d, compute_bead_rms_spread_1d
  use ereaxff_output_manager, only: join_path
  use ereaxff_diagnostics, only: nvt_conserved_like
  implicit none

  integer, parameter :: p = 16
  integer, parameter :: steps = 2000
  integer, parameter :: output_stride = 100
  real(dp), parameter :: temperature_K = 300.0_dp
  real(dp), parameter :: mass_amu = 1.0_dp
  real(dp), parameter :: k_eV_per_A2 = 1.0_dp
  real(dp), parameter :: dt_fs = 0.2_dp
  real(dp), parameter :: centroid_gamma_fs_inv = 0.02_dp
  integer(i8), parameter :: seed = 13579_i8

  type(ring_polymer_state) :: rp
  type(pile_l_thermostat_state) :: thermostat
  type(rng_state) :: rng
  real(dp) :: q_nm(p), v_nm(p), q_bead(p), q2d(1, p), v2d(1, p), masses(1)
  real(dp) :: centroid_A, spread_A, kinetic_eV, spring_eV, physical_potential_eV
  real(dp) :: conserved_like, finite_p_variance
  character(len=512) :: output_dir
  character(len=1100) :: mkdir_command
  integer :: step, energy_unit, spread_unit, heat_unit, xyz_unit

  call get_command_argument(1, output_dir)
  if (len_trim(output_dir) == 0) output_dir = "harmonic_pimd_output"

  mkdir_command = 'mkdir -p "' // trim(output_dir) // '"'
  call execute_command_line(trim(mkdir_command))

  write(*,'(a)') "PI-eReaxFF harmonic PIMD driver"
  write(*,'(a)') "build/version: development"
  write(*,'(a)') "ring-polymer convention: standard_ipi"
  write(*,'(a,i0)') "P: ", p
  write(*,'(a,es12.4)') "T_K: ", temperature_K
  write(*,'(a)') "thermostat mode: pile_l"
  write(*,'(a)') "charge coupling mode: no_charge_coupling"
  write(*,'(a)') "force provider: harmonic"
  write(*,'(a)') "output directory: " // trim(output_dir)

  call init_ring_polymer(rp, p, temperature_K, "standard_ipi")
  call pile_l_init(thermostat, temperature_K * real(p, dp), centroid_gamma_fs_inv, 1.0_dp, dt_fs, rp%omega_k)
  call rng_seed(rng, seed)

  q_nm = 0.0_dp
  v_nm = 0.0_dp
  masses = [mass_amu]
  finite_p_variance = harmonic_finite_p_bead_variance_A2(k_eV_per_A2, mass_amu, temperature_K, rp%omega_k)

  open(newunit=energy_unit, file=join_path(output_dir, "pimd_energy.out"), status="replace", action="write")
  open(newunit=spread_unit, file=join_path(output_dir, "bead_spread.out"), status="replace", action="write")
  open(newunit=heat_unit, file=join_path(output_dir, "thermostat_heat.out"), status="replace", action="write")
  open(newunit=xyz_unit, file=join_path(output_dir, "centroid.xyz"), status="replace", action="write")

  write(energy_unit,'(a)') "# step kinetic_eV physical_potential_eV spring_eV " // &
      "thermostat_heat_eV conserved_like_eV finite_P_variance_A2"
  write(spread_unit,'(a)') "# step centroid_A bead_rms_spread_A"
  write(heat_unit,'(a)') "# step thermostat_heat_eV"

  do step = 0, steps
    if (mod(step, output_stride) == 0) then
      call nm_to_bead(rp%cmat, q_nm, q_bead)
      centroid_A = compute_centroid_1d(q_bead)
      spread_A = compute_bead_rms_spread_1d(q_bead)
      q2d(1, :) = q_nm
      v2d(1, :) = v_nm
      kinetic_eV = kinetic_energy_2d_eV(v2d, masses)
      spring_eV = compute_spring_energy_nm(q2d, masses, rp%omega_k)
      physical_potential_eV = 0.5_dp * k_eV_per_A2 * sum(q_bead**2) / real(p, dp)
      conserved_like = nvt_conserved_like(kinetic_eV, physical_potential_eV, spring_eV, thermostat%heat_eV)

      write(energy_unit,'(i8,6(1x,es22.14))') step, kinetic_eV, physical_potential_eV, spring_eV, &
          thermostat%heat_eV, conserved_like, finite_p_variance
      write(spread_unit,'(i8,2(1x,es22.14))') step, centroid_A, spread_A
      write(heat_unit,'(i8,1x,es22.14)') step, thermostat%heat_eV
      call write_centroid_xyz_frame(xyz_unit, step, centroid_A)
    end if

    if (step < steps) then
      call baoab_harmonic_pimd_step_1d(q_nm, v_nm, rp%cmat, rp%omega_k, mass_amu, &
          k_eV_per_A2, dt_fs, thermostat, rng, physical_force_scale=1.0_dp)
    end if
  end do

  close(energy_unit)
  close(spread_unit)
  close(heat_unit)
  close(xyz_unit)
end program harmonic_pimd_main
