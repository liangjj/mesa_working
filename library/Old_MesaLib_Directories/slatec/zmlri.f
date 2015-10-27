*deck zmlri
      subroutine zmlri (zr, zi, fnu, kode, n, yr, yi, nz, tol)
c***begin prologue  zmlri
c***subsidiary
c***purpose  subsidiary to zbesi and zbesk
c***library   slatec
c***type      all (cmlri-a, zmlri-a)
c***author  amos, d. e., (snl)
c***description
c
c     zmlri computes the i bessel function for re(z).ge.0.0 by the
c     miller algorithm normalized by a neumann series.
c
c***see also  zbesi, zbesk
c***routines called  d1mach, dgamln, zabs, zexp, zlog, zmlt
c***revision history  (yymmdd)
c   830501  date written
c   910415  prologue converted to version 4.0 format.  (bab)
c   930122  added zexp and zlog to external statement.  (rwc)
c***end prologue  zmlri
c     complex ck,cnorm,cone,ctwo,czero,pt,p1,p2,rz,sum,y,z
      double precision ack, ak, ap, at, az, bk, cki, ckr, cnormi,
     * cnormr, conei, coner, fkap, fkk, flam, fnf, fnu, pti, ptr, p1i,
     * p1r, p2i, p2r, raz, rho, rho2, rzi, rzr, scle, sti, str, sumi,
     * sumr, tfnf, tol, tst, yi, yr, zeroi, zeror, zi, zr, dgamln,
     * d1mach, zabs
      integer i, iaz, idum, ifnu, inu, itime, k, kk, km, kode, m, n, nz
      dimension yr(n), yi(n)
      external zabs, zexp, zlog
      data zeror,zeroi,coner,conei / 0.0d0, 0.0d0, 1.0d0, 0.0d0 /
c***first executable statement  zmlri
      scle = d1mach(1)/tol
      nz=0
      az = zabs(zr,zi)
      iaz = az
      ifnu = fnu
      inu = ifnu + n - 1
      at = iaz + 1.0d0
      raz = 1.0d0/az
      str = zr*raz
      sti = -zi*raz
      ckr = str*at*raz
      cki = sti*at*raz
      rzr = (str+str)*raz
      rzi = (sti+sti)*raz
      p1r = zeror
      p1i = zeroi
      p2r = coner
      p2i = conei
      ack = (at+1.0d0)*raz
      rho = ack + sqrt(ack*ack-1.0d0)
      rho2 = rho*rho
      tst = (rho2+rho2)/((rho2-1.0d0)*(rho-1.0d0))
      tst = tst/tol
c-----------------------------------------------------------------------
c     compute relative truncation error index for series
c-----------------------------------------------------------------------
      ak = at
      do 10 i=1,80
        ptr = p2r
        pti = p2i
        p2r = p1r - (ckr*ptr-cki*pti)
        p2i = p1i - (cki*ptr+ckr*pti)
        p1r = ptr
        p1i = pti
        ckr = ckr + rzr
        cki = cki + rzi
        ap = zabs(p2r,p2i)
        if (ap.gt.tst*ak*ak) go to 20
        ak = ak + 1.0d0
   10 continue
      go to 110
   20 continue
      i = i + 1
      k = 0
      if (inu.lt.iaz) go to 40
c-----------------------------------------------------------------------
c     compute relative truncation error for ratios
c-----------------------------------------------------------------------
      p1r = zeror
      p1i = zeroi
      p2r = coner
      p2i = conei
      at = inu + 1.0d0
      str = zr*raz
      sti = -zi*raz
      ckr = str*at*raz
      cki = sti*at*raz
      ack = at*raz
      tst = sqrt(ack/tol)
      itime = 1
      do 30 k=1,80
        ptr = p2r
        pti = p2i
        p2r = p1r - (ckr*ptr-cki*pti)
        p2i = p1i - (ckr*pti+cki*ptr)
        p1r = ptr
        p1i = pti
        ckr = ckr + rzr
        cki = cki + rzi
        ap = zabs(p2r,p2i)
        if (ap.lt.tst) go to 30
        if (itime.eq.2) go to 40
        ack = zabs(ckr,cki)
        flam = ack + sqrt(ack*ack-1.0d0)
        fkap = ap/zabs(p1r,p1i)
        rho = min(flam,fkap)
        tst = tst*sqrt(rho/(rho*rho-1.0d0))
        itime = 2
   30 continue
      go to 110
   40 continue
c-----------------------------------------------------------------------
c     backward recurrence and sum normalizing relation
c-----------------------------------------------------------------------
      k = k + 1
      kk = max(i+iaz,k+inu)
      fkk = kk
      p1r = zeror
      p1i = zeroi
c-----------------------------------------------------------------------
c     scale p2 and sum by scle
c-----------------------------------------------------------------------
      p2r = scle
      p2i = zeroi
      fnf = fnu - ifnu
      tfnf = fnf + fnf
      bk = dgamln(fkk+tfnf+1.0d0,idum) - dgamln(fkk+1.0d0,idum) -
     * dgamln(tfnf+1.0d0,idum)
      bk = exp(bk)
      sumr = zeror
      sumi = zeroi
      km = kk - inu
      do 50 i=1,km
        ptr = p2r
        pti = p2i
        p2r = p1r + (fkk+fnf)*(rzr*ptr-rzi*pti)
        p2i = p1i + (fkk+fnf)*(rzi*ptr+rzr*pti)
        p1r = ptr
        p1i = pti
        ak = 1.0d0 - tfnf/(fkk+tfnf)
        ack = bk*ak
        sumr = sumr + (ack+bk)*p1r
        sumi = sumi + (ack+bk)*p1i
        bk = ack
        fkk = fkk - 1.0d0
   50 continue
      yr(n) = p2r
      yi(n) = p2i
      if (n.eq.1) go to 70
      do 60 i=2,n
        ptr = p2r
        pti = p2i
        p2r = p1r + (fkk+fnf)*(rzr*ptr-rzi*pti)
        p2i = p1i + (fkk+fnf)*(rzi*ptr+rzr*pti)
        p1r = ptr
        p1i = pti
        ak = 1.0d0 - tfnf/(fkk+tfnf)
        ack = bk*ak
        sumr = sumr + (ack+bk)*p1r
        sumi = sumi + (ack+bk)*p1i
        bk = ack
        fkk = fkk - 1.0d0
        m = n - i + 1
        yr(m) = p2r
        yi(m) = p2i
   60 continue
   70 continue
      if (ifnu.le.0) go to 90
      do 80 i=1,ifnu
        ptr = p2r
        pti = p2i
        p2r = p1r + (fkk+fnf)*(rzr*ptr-rzi*pti)
        p2i = p1i + (fkk+fnf)*(rzr*pti+rzi*ptr)
        p1r = ptr
        p1i = pti
        ak = 1.0d0 - tfnf/(fkk+tfnf)
        ack = bk*ak
        sumr = sumr + (ack+bk)*p1r
        sumi = sumi + (ack+bk)*p1i
        bk = ack
        fkk = fkk - 1.0d0
   80 continue
   90 continue
      ptr = zr
      pti = zi
      if (kode.eq.2) ptr = zeror
      call zlog(rzr, rzi, str, sti, idum)
      p1r = -fnf*str + ptr
      p1i = -fnf*sti + pti
      ap = dgamln(1.0d0+fnf,idum)
      ptr = p1r - ap
      pti = p1i
c-----------------------------------------------------------------------
c     the division cexp(pt)/(sum+p2) is altered to avoid overflow
c     in the denominator by squaring large quantities
c-----------------------------------------------------------------------
      p2r = p2r + sumr
      p2i = p2i + sumi
      ap = zabs(p2r,p2i)
      p1r = 1.0d0/ap
      call zexp(ptr, pti, str, sti)
      ckr = str*p1r
      cki = sti*p1r
      ptr = p2r*p1r
      pti = -p2i*p1r
      call zmlt(ckr, cki, ptr, pti, cnormr, cnormi)
      do 100 i=1,n
        str = yr(i)*cnormr - yi(i)*cnormi
        yi(i) = yr(i)*cnormi + yi(i)*cnormr
        yr(i) = str
  100 continue
      return
  110 continue
      nz=-2
      return
      end
