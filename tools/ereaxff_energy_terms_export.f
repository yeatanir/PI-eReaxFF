      subroutine ereaxff_oracle_get_energy_terms(terms_kcal,
     $  total_kcal)
      include 'cbka_v2.blk'
      double precision terms_kcal(14), total_kcal

      terms_kcal(1) = eb
      terms_kcal(2) = ea
      terms_kcal(3) = elp
      terms_kcal(4) = emol
      terms_kcal(5) = ev + epen
      terms_kcal(6) = ecoa
      terms_kcal(7) = ehb
      terms_kcal(8) = et
      terms_kcal(9) = eco
      terms_kcal(10) = ew
      terms_kcal(11) = ep
      terms_kcal(12) = ech
      terms_kcal(13) = efi
      terms_kcal(14) = eelect
      total_kcal = estrc

      return
      end

      subroutine ereaxff_oracle_get_bond_orders(natom_eval,
     $  nbond_eval, abo_out, bond_atoms_out, bo_out, bopi_out,
     $  bopi2_out)
      include 'cbka_v2.blk'
      integer natom_eval, nbond_eval
      integer bond_atoms_out(2,*)
      double precision abo_out(*), bo_out(*), bopi_out(*), bopi2_out(*)
      integer i1

      natom_eval = na
      nbond_eval = nbon
      do i1 = 1, na
        abo_out(i1) = abo(i1)
      end do
      do i1 = 1, nbon
        bond_atoms_out(1,i1) = ib(i1,2)
        bond_atoms_out(2,i1) = ib(i1,3)
        bo_out(i1) = bo(i1)
        bopi_out(i1) = bopi(i1)
        bopi2_out(i1) = bopi2(i1)
      end do

      return
      end
