program oracle_single_point_force_export
  implicit none
  integer, parameter :: dp = selected_real_kind(15, 300)
  integer, parameter :: natmax = 2600
  integer, parameter :: nbondmax = 400000
  integer, parameter :: nbucketmax = 256
  integer :: natom, status, i, nbond, nbuckets
  integer :: bond_atoms(2,nbondmax), bo_natom
  integer :: bucket_left(nbucketmax), bucket_right(nbucketmax), bucket_self(nbucketmax), bucket_count(nbucketmax)
  real(dp) :: x(3,natmax), force(3,natmax), charge_in(natmax), solver_in(natmax)
  real(dp) :: charge_out(natmax), solver_out(natmax), energy
  real(dp) :: energy_terms(14), energy_terms_total
  real(dp) :: abo(natmax), bo(nbondmax), bopi(nbondmax), bopi2(nbondmax)
  real(dp) :: bucket_vdw(nbucketmax), bucket_coulomb(nbucketmax)
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
    subroutine ereaxff_oracle_get_energy_terms(terms_kcal, total_kcal)
      double precision :: terms_kcal(14), total_kcal
    end subroutine ereaxff_oracle_get_energy_terms
    subroutine ereaxff_oracle_get_bond_orders(natom_eval, nbond_eval, abo_out, bond_atoms_out, bo_out, bopi_out, bopi2_out)
      integer :: natom_eval, nbond_eval
      integer :: bond_atoms_out(2,*)
      double precision :: abo_out(*), bo_out(*), bopi_out(*), bopi2_out(*)
    end subroutine ereaxff_oracle_get_bond_orders
    subroutine ereaxff_oracle_get_nonbonded_buckets(max_rows, nrows, left_type_out, right_type_out, self_out, count_out, vdw_out, coulomb_out)
      integer :: max_rows, nrows
      integer :: left_type_out(*), right_type_out(*), self_out(*), count_out(*)
      double precision :: vdw_out(*), coulomb_out(*)
    end subroutine ereaxff_oracle_get_nonbonded_buckets
  end interface

  call ereaxff_force_oracle_init()
  call ereaxff_oracle_get_geometry(x, natom)
  charge_in = 0.0d0
  solver_in = 0.0d0
  call ereaxff_force_eval(x, natom, charge_in, solver_in, energy, force, &
      charge_out, solver_out, status)
  if (status /= 0) stop 'oracle_single_point_force_export: oracle failure'
  call ereaxff_oracle_get_energy_terms(energy_terms, energy_terms_total)
  call ereaxff_oracle_get_bond_orders(bo_natom, nbond, abo, bond_atoms, bo, bopi, bopi2)
  call ereaxff_oracle_get_nonbonded_buckets(nbucketmax, nbuckets, bucket_left, bucket_right, bucket_self, &
      bucket_count, bucket_vdw, bucket_coulomb)
  if (bo_natom /= natom) stop 'oracle_single_point_force_export: bond-order atom count mismatch'
  if (nbond > nbondmax) stop 'oracle_single_point_force_export: too many bond-order rows'

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
  call write_energy_terms('single_point_energy_terms.out', energy_terms, energy_terms_total)
  call write_bond_orders('single_point_bond_orders.out', natom, nbond, abo, bond_atoms, bo, bopi, bopi2)
  call write_nonbonded_buckets('single_point_nonbonded_buckets.out', nbuckets, bucket_left, bucket_right, &
      bucket_self, bucket_count, bucket_vdw, bucket_coulomb)

  write(*,'(a,1x,es12.4,1x,es12.4,1x,es12.4,1x,i0)') &
      'oracle_single_point_force_export PASS', energy, max_force, charge_sum, status

contains

  subroutine write_energy_terms(path, terms_kcal, total_kcal)
    character(len=*), intent(in) :: path
    real(dp), intent(in) :: terms_kcal(14), total_kcal
    character(len=32), parameter :: names(14) = [ &
        character(len=32) :: &
        'bond', 'overcoordination', 'lone_pair', 'molecule', &
        'valence_penalty', 'angle_conjugation', 'hydrogen_bond', &
        'torsion', 'torsion_conjugation', 'vdw', 'coulomb', &
        'charge', 'external_field', 'electron' &
    ]
    integer :: unit, i

    open(newunit=unit, file=trim(path), status='replace', action='write')
    write(unit,'(a)') '# term energy_kcal'
    do i = 1, size(names)
      write(unit,'(a,1x,es24.15)') trim(names(i)), terms_kcal(i)
    end do
    write(unit,'(a,1x,es24.15)') 'total', total_kcal
    close(unit)
  end subroutine write_energy_terms

  subroutine write_bond_orders(path, natom_in, nbond_in, abo_in, bond_atoms_in, bo_in, bopi_in, bopi2_in)
    character(len=*), intent(in) :: path
    integer, intent(in) :: natom_in, nbond_in, bond_atoms_in(2,*)
    real(dp), intent(in) :: abo_in(*), bo_in(*), bopi_in(*), bopi2_in(*)
    integer :: unit, i

    open(newunit=unit, file=trim(path), status='replace', action='write')
    write(unit,'(a)') '# kind index atom_i atom_j value bopi bopi2'
    do i = 1, natom_in
      write(unit,'(a,1x,i8,1x,es24.15)') 'ATOM', i, abo_in(i)
    end do
    do i = 1, nbond_in
      write(unit,'(a,1x,i8,2(1x,i8),3es24.15)') 'BOND', i, bond_atoms_in(1,i), bond_atoms_in(2,i), &
          bo_in(i), bopi_in(i), bopi2_in(i)
    end do
    close(unit)
  end subroutine write_bond_orders

  subroutine write_nonbonded_buckets(path, nrows, left_type, right_type, self_image, bucket_count, vdw_kcal, coulomb_kcal)
    character(len=*), intent(in) :: path
    integer, intent(in) :: nrows, left_type(*), right_type(*), self_image(*), bucket_count(*)
    real(dp), intent(in) :: vdw_kcal(*), coulomb_kcal(*)
    integer :: unit, i

    open(newunit=unit, file=trim(path), status='replace', action='write')
    write(unit,'(a)') '# left_type right_type self_image count vdw_kcal coulomb_kcal'
    do i = 1, nrows
      write(unit,'(4i8,2es24.15)') left_type(i), right_type(i), self_image(i), bucket_count(i), vdw_kcal(i), coulomb_kcal(i)
    end do
    close(unit)
  end subroutine write_nonbonded_buckets
end program oracle_single_point_force_export
