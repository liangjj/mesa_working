*deck crati
      subroutine crati (z, fnu, n, cy, tol)
c***begin prologue  crati
c***subsidiary
c***purpose  subsidiary to cbesh, cbesi and cbesk
c***library   slatec
c***type      all (crati-a, zrati-a)
c***author  amos, d. e., (snl)
c***description
c
c     crati computes ratios of i bessel functions by backward
c     recurrence.  the starting index is determined by forward
c     recurrence as described in j. res. of nat. bur. of standards-b,
c     mathematical sciences, vol 77b, p111-114, september, 1973,
c     bessel functions i and j of complex argument and integer order,
c     by d. j. sookne.
c
c***see also  cbesh, cbesi, cbesk
c***routines called  (none)
c***revision history  (yymmdd)
c   830501  date written
c   910415  prologue converted to version 4.0 format.  (bab)
c***end prologue  crati
      complex cdfnu, cone, cy, czero, pt, p1, p2, rz, t1, z
      real ak, amagz, ap1, ap2, arg, az, dfnu, fdnu, flam, fnu, fnup,
     * rap1, rho, test, test1, tol
      integer i, id, idnu, inu, itime, k, kk, magz, n
      dimension cy(n)
      data czero, cone / (0.0e0,0.0e0), (1.0e0,0.0e0) /
c***first executable statement  crati
      az = abs(z)
      inu = fnu
      idnu = inu + n - 1
      fdnu = idnu
      magz = az
      amagz = magz+1
      fnup = max(amagz,fdnu)
      id = idnu - magz - 1
      itime = 1
      k = 1
      rz = (cone+cone)/z
      t1 = cmplx(fnup,0.0e0)*rz
      p2 = -t1
      p1 = cone
      t1 = t1 + rz
      if (id.gt.0) id = 0
      ap2 = abs(p2)
      ap1 = abs(p1)
c-----------------------------------------------------------------------
c     the overflow test on k(fnu+i-1,z) before the call to cbknx
c     guarantees that p2 is on scale. scale test1 and all subsequent
c     p2 values by ap1 to ensure that an overflow does not occur
c     prematurely.
c-----------------------------------------------------------------------
      arg = (ap2+ap2)/(ap1*tol)
      test1 = sqrt(arg)
      test = test1
      rap1 = 1.0e0/ap1
      p1 = p1*cmplx(rap1,0.0e0)
      p2 = p2*cmplx(rap1,0.0e0)
      ap2 = ap2*rap1
   10 continue
      k = k + 1
      ap1 = ap2
      pt = p2
      p2 = p1 - t1*p2
      p1 = pt
      t1 = t1 + rz
      ap2 = abs(p2)
      if (ap1.le.test) go to 10
      if (itime.eq.2) go to 20
      ak = abs(t1)*0.5e0
      flam = ak + sqrt(ak*ak-1.0e0)
      rho = min(ap2/ap1,flam)
      test = test1*sqrt(rho/(rho*rho-1.0e0))
      itime = 2
      go to 10
   20 continue
      kk = k + 1 - id
      ak = kk
      dfnu = fnu + (n-1)
      cdfnu = cmplx(dfnu,0.0e0)
      t1 = cmplx(ak,0.0e0)
      p1 = cmplx(1.0e0/ap2,0.0e0)
      p2 = czero
      do 30 i=1,kk
        pt = p1
        p1 = rz*(cdfnu+t1)*p1 + p2
        p2 = pt
        t1 = t1 - cone
   30 continue
      if (real(p1).ne.0.0e0 .or. aimag(p1).ne.0.0e0) go to 40
      p1 = cmplx(tol,tol)
   40 continue
      cy(n) = p2/p1
      if (n.eq.1) return
      k = n - 1
      ak = k
      t1 = cmplx(ak,0.0e0)
      cdfnu = cmplx(fnu,0.0e0)*rz
      do 60 i=2,n
        pt = cdfnu + t1*rz + cy(k+1)
        if (real(pt).ne.0.0e0 .or. aimag(pt).ne.0.0e0) go to 50
        pt = cmplx(tol,tol)
   50   continue
        cy(k) = cone/pt
        t1 = t1 - cone
        k = k - 1
   60 continue
      return
      end
