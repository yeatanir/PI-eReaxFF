program oracle_single_point_force_export
  implicit none
  integer, parameter :: dp = selected_real_kind(15, 300)
  integer, parameter :: natmax = 2600
  integer :: natom, status, i
  real(dp) :: x(3,natmax), force(3,natmax), charge_in(natmax), solver_in(natmax)
  real(dp) :: charge_out(natmax), solver_out(natmax), energy
  real(dp) :: max_force, charge_sum

  interface
    subroutine ereaxff_force_oracle_init()
    end subroutine ereaxff_force_oracle_init
    subroutine ereaxff_oracle_get_geometry(x_all, natom_eval)
      integer :: natom_eval
      double precision :: x_all(3,*)
    end subroutine ereaxff_oracle_get_geometry
    subroutine ereaxff_force_eval(x_all, natom_eval, charge_in, &
        acks2_solver_guess_in, energy_kcal, forces_work, charge_out, &
        acks2_solver_guess_out, status)
      integer :: natom_eval, status
      double precision :: x_all(3,*), charge_in(*)
      double precision :: acks2_solver_guess_in(*), energy_kcal
      double precision :: forces_work(3,*), charge_out(*)
      double precision :: acks2_solver_guess_out(*)
    end subroutine ereaxff_force_eval
  end interface

  call ereaxff_force_oracle_init()
  call ereaxff_oracle_get_geometry(x, natom)
  charge_in = 0.0d0
  solver_in = 0.0d0
  call ereaxff_force_eval(x, natom, charge_in, solver_in, energy, force, &
      charge_out, solver_out, status)
  if (status /= 0) stop 'oracle_single_point_force_export: oracle failure'

  max_force = 0.0d0
  charge_sum = 0.0d0
  do i = 1, natom
    max_force = max(max_force, abs(force(1,i)))
    max_force = max(max_force, abs(force(2,i)))
    max_force = max(max_force, abs(force(3,i)))
    charge_sum = charge_sum + charge_out(i)
  end do

  open(20, file='phase4_single_point_oracle.out', status='replace')
  write(20,'(a)') '# energy_kcal max_force charge_sum status natom'
  write(20,'(es24.15,1x,es24.15,1x,es24.15,1x,i8,1x,i8)') &
      energy, max_force, charge_sum, status, natom
  close(20)

  open(21, file='single_point_forces.out', status='replace')
  write(21,'(a)') '# atom Fx Fy Fz charge solver_guess'
  do i = 1, natom
    write(21,'(i8,5es24.15)') i, force(1,i), force(2,i), force(3,i), charge_out(i), solver_out(i)
  end do
  close(21)

  write(*,'(a,1x,es12.4,1x,es12.4,1x,es12.4,1x,i0)') &
      'oracle_single_point_force_export PASS', energy, max_force, charge_sum, status
end program oracle_single_point_force_export
