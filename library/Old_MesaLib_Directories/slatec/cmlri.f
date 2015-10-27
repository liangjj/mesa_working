*deck cmlri
      subroutine cmlri (z, fnu, kode, n, y, nz, tol)
c***begin prologue  cmlri
c***subsidiary
c***purpose  subsidiary to cbesi and cbesk
c***library   slatec
c***type      all (cmlri-a, zmlri-a)
c***author  amos, d. e., (snl)
c***description
c
c     cmlri computes the i bessel function for re(z).ge.0.0 by the
c     miller algorithm normalized by a neumann series.
c
c***see also  cbesi, cbesk
c***routines called  gamln, r1mach
c***revision history  (yymmdd)
c   830501  date written
c   910415  prologue converted to version 4.0 format.  (bab)
c***end prologue  cmlri
      complex ck, cnorm, cone, ctwo, czero, pt, p1, p2, rz, sum, y, z
      real ack, ak, ap, at, az, bk, fkap, fkk, flam, fnf, fnu, rho,
     * rho2, scle, tfnf, tol, tst, x, gamln, r1mach
      integer i, iaz, idum, ifnu, inu, itime, k, kk, km, kode, m, n, nz
      dimension y(n)
      data czero,cone,ctwo /(0.0e0,0.0e0),(1.0e0,0.0e0),(2.0e0,0.0e0)/
      scle = 1.0e+3*r1mach(1)/tol
c***first executable statement  cmlri
      nz=0
      az = abs(z)
      x = real(z)
      iaz = az
      ifnu = fnu
      inu = ifnu + n - 1
      at = iaz + 1.0e0
      ck = cmplx(at,0.0e0)/z
      rz = ctwo/z
      p1 = czero
      p2 = cone
      ack = (at+1.0e0)/az
      rho = ack + sqrt(ack*ack-1.0e0)
      rho2 = rho*rho
      tst = (rho2+rho2)/((rho2-1.0e0)*(rho-1.0e0))
      tst = tst/tol
c-----------------------------------------------------------------------
c     compute relative truncation error index for series
c-----------------------------------------------------------------------
      ak = at
      do 10 i=1,80
        pt = p2
        p2 = p1 - ck*p2
        p1 = pt
        ck = ck + rz
        ap = abs(p2)
        if (ap.gt.tst*ak*ak) go to 20
        ak = ak + 1.0e0
   10 continue
      go to 110
   20 continue
      i = i + 1
      k = 0
      if (inu.lt.iaz) go to 40
c-----------------------------------------------------------------------
c     compute relative truncation error for ratios
c-----------------------------------------------------------------------
      p1 = czero
      p2 = cone
      at = inu + 1.0e0
      ck = cmplx(at,0.0e0)/z
      ack = at/az
      tst = sqrt(ack/tol)
      itime = 1
      do 30 k=1,80
        pt = p2
        p2 = p1 - ck*p2
        p1 = pt
        ck = ck + rz
        ap = abs(p2)
        if (ap.lt.tst) go to 30
        if (itime.eq.2) go to 40
        ack = abs(ck)
        flam = ack + sqrt(ack*ack-1.0e0)
        fkap = ap/abs(p1)
        rho = min(flam,fkap)
        tst = tst*sqrt(rho/(rho*rho-1.0e0))
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
      p1 = czero
c-----------------------------------------------------------------------
c     scale p2 and sum by scle
c-----------------------------------------------------------------------
      p2 = cmplx(scle,0.0e0)
      fnf = fnu - ifnu
      tfnf = fnf + fnf
      bk = gamln(fkk+tfnf+1.0e0,idum) - gamln(fkk+1.0e0,idum)
     *     -gamln(tfnf+1.0e0,idum)
      bk = exp(bk)
      sum = czero
      km = kk - inu
      do 50 i=1,km
        pt = p2
        p2 = p1 + cmplx(fkk+fnf,0.0e0)*rz*p2
        p1 = pt
        ak = 1.0e0 - tfnf/(fkk+tfnf)
        ack = bk*ak
        sum = sum + cmplx(ack+bk,0.0e0)*p1
        bk = ack
        fkk = fkk - 1.0e0
   50 continue
      y(n) = p2
      if (n.eq.1) go to 70
      do 60 i=2,n
        pt = p2
        p2 = p1 + cmplx(fkk+fnf,0.0e0)*rz*p2
        p1 = pt
        ak = 1.0e0 - tfnf/(fkk+tfnf)
        ack = bk*ak
        sum = sum + cmplx(ack+bk,0.0e0)*p1
        bk = ack
        fkk = fkk - 1.0e0
        m = n - i + 1
        y(m) = p2
   60 continue
   70 continue
      if (ifnu.le.0) go to 90
      do 80 i=1,ifnu
        pt = p2
        p2 = p1 + cmplx(fkk+fnf,0.0e0)*rz*p2
        p1 = pt
        ak = 1.0e0 - tfnf/(fkk+tfnf)
        ack = bk*ak
        sum = sum + cmplx(ack+bk,0.0e0)*p1
        bk = ack
        fkk = fkk - 1.0e0
   80 continue
   90 continue
      pt = z
      if (kode.eq.2) pt = pt - cmplx(x,0.0e0)
      p1 = -cmplx(fnf,0.0e0)*clog(rz) + pt
      ap = gamln(1.0e0+fnf,idum)
      pt = p1 - cmplx(ap,0.0e0)
c-----------------------------------------------------------------------
c     the division cexp(pt)/(sum+p2) is altered to avoid overflow
c     in the denominator by squaring large quantities
c-----------------------------------------------------------------------
      p2 = p2 + sum
      ap = abs(p2)
      p1 = cmplx(1.0e0/ap,0.0e0)
      ck = cexp(pt)*p1
      pt = conjg(p2)*p1
      cnorm = ck*pt
      do 100 i=1,n
        y(i) = y(i)*cnorm
  100 continue
      return
  110 continue
      nz=-2
      return
      end
