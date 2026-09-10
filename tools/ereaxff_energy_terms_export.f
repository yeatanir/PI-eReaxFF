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
