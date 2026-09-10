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

      subroutine ereaxff_oracle_get_nonbonded_buckets(max_rows,
     $  nrows, left_type_out, right_type_out, self_out, count_out,
     $  vdw_out, coulomb_out)
      include 'cbka_v2.blk'
      integer max_rows, nrows
      integer left_type_out(*), right_type_out(*), self_out(*)
      integer count_out(*)
      double precision vdw_out(*), coulomb_out(*)
      integer ivl, i1, i2, ix, iy, iz, ity1, ity2, row, iself
      double precision ax, ay, az, rr, rr2, p1, p2, p3
      double precision hulpw, rrw, h1, h2, ewh, q1q2, hulp1, eph
      double precision c1c_local, third_local

      nrows = 0
      c1c_local = 332.0638d0
      third_local = 1.0d0 / 3.0d0
      do ivl = 1, max_rows
        left_type_out(ivl) = 0
        right_type_out(ivl) = 0
        self_out(ivl) = 0
        count_out(ivl) = 0
        vdw_out(ivl) = 0.0d0
        coulomb_out(ivl) = 0.0d0
      end do

      do 10 ivl = 1, nvpair - nvlself
        i1 = nvl1(ivl)
        i2 = nvl2(ivl)
        ix = nvlx(ivl)
        iy = nvly(ivl)
        iz = nvlz(ivl)
        ax = c(i1,1) - c(i2,1) + ix * tm11
        ay = c(i1,2) - c(i2,2) + ix * tm21 + iy * tm22
        az = c(i1,3) - c(i2,3) + ix * tm31 + iy * tm32
     $    + iz * tm33
        rr = sqrt(ax * ax + ay * ay + az * az)
        if (rr.gt.swb.or.rr.lt.0.000000000001d0) goto 10
        ity1 = ia(i1,1)
        ity2 = ia(i2,1)
        iself = 0
        call nbucket_find(max_rows, nrows, left_type_out,
     $    right_type_out, self_out, count_out, vdw_out, coulomb_out,
     $    ity1, ity2, iself, row)
        rr2 = rr * rr
        sw = 1.0d0
        sw1 = 0.0d0
        call taper(rr, rr2)
        p1 = p1co(ity1,ity2)
        p2 = p2co(ity1,ity2)
        p3 = p3co(ity1,ity2)
        hulpw = rr**vpar(29) + gamwco(ity1,ity2)
        rrw = hulpw**(1.0d0 / vpar(29))
        h1 = exp(p3 * (1.0d0 - rrw / p1))
        h2 = exp(0.50d0 * p3 * (1.0d0 - rrw / p1))
        ewh = p2 * (h1 - 2.0d0 * h2)
        q1q2 = ch(i1) * ch(i2)
        hulp1 = rr2 * rr + gamcco(ity1,ity2)
        eph = c1c_local * q1q2 / (hulp1**third_local)
        count_out(row) = count_out(row) + 1
        vdw_out(row) = vdw_out(row) + ewh * sw
        coulomb_out(row) = coulomb_out(row) + eph * sw
   10 continue

      do 20 ivl = nvpair - nvlself, nvpair
        i1 = nvl1(ivl)
        ix = nvlx(ivl)
        iy = nvly(ivl)
        iz = nvlz(ivl)
        ax = ix * tm11
        ay = ix * tm21 + iy * tm22
        az = ix * tm31 + iy * tm33
        rr = sqrt(ax * ax + ay * ay + az * az)
        if (rr.gt.swb.or.rr.lt.0.000000000001d0) goto 20
        ity1 = ia(i1,1)
        ity2 = ity1
        iself = 1
        call nbucket_find(max_rows, nrows, left_type_out,
     $    right_type_out, self_out, count_out, vdw_out, coulomb_out,
     $    ity1, ity2, iself, row)
        rr2 = rr * rr
        sw = 1.0d0
        sw1 = 0.0d0
        call taper(rr, rr2)
        p1 = p1co(ity1,ity1)
        p2 = p2co(ity1,ity1)
        p3 = p3co(ity1,ity1)
        hulpw = rr**vpar(29) + gamwco(ity1,ity1)
        rrw = hulpw**(1.0d0 / vpar(29))
        h1 = exp(p3 * (1.0d0 - rrw / p1))
        h2 = exp(0.50d0 * p3 * (1.0d0 - rrw / p1))
        ewh = 0.50d0 * p2 * (h1 - 2.0d0 * h2)
        q1q2 = ch(i1) * ch(i1)
        hulp1 = rr * rr2 + gamcco(ity1,ity1)
        eph = 0.50d0 * c1c_local * q1q2
     $    / (hulp1**third_local)
        count_out(row) = count_out(row) + 1
        vdw_out(row) = vdw_out(row) + ewh * sw
        coulomb_out(row) = coulomb_out(row) + eph * sw
   20 continue

      return
      end

      subroutine nbucket_find(max_rows, nrows, left_type_out,
     $  right_type_out, self_out, count_out, vdw_out, coulomb_out,
     $  ity1, ity2, iself, row)
      integer max_rows, nrows, ity1, ity2, iself, row
      integer left_type_out(*), right_type_out(*), self_out(*)
      integer count_out(*)
      double precision vdw_out(*), coulomb_out(*)
      integer left_type, right_type, i

      left_type = ity1
      right_type = ity2
      if (iself.eq.0.and.left_type.gt.right_type) then
        left_type = ity2
        right_type = ity1
      end if
      do i = 1, nrows
        if (left_type_out(i).eq.left_type.and.
     $      right_type_out(i).eq.right_type.and.
     $      self_out(i).eq.iself) then
          row = i
          return
        end if
      end do
      if (nrows.ge.max_rows) stop 'Too many nonbonded buckets'
      nrows = nrows + 1
      row = nrows
      left_type_out(row) = left_type
      right_type_out(row) = right_type
      self_out(row) = iself
      count_out(row) = 0
      vdw_out(row) = 0.0d0
      coulomb_out(row) = 0.0d0
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
