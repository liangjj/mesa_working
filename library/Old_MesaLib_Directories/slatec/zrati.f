*deck zrati
      subroutine zrati (zr, zi, fnu, n, cyr, cyi, tol)
c***begin prologue  zrati
c***subsidiary
c***purpose  subsidiary to zbesh, zbesi and zbesk
c***library   slatec
c***type      all (crati-a, zrati-a)
c***author  amos, d. e., (snl)
c***description
c
c     zrati computes ratios of i bessel functions by backward
c     recurrence.  the starting index is determined by forward
c     recurrence as described in j. res. of nat. bur. of standards-b,
c     mathematical sciences, vol 77b, p111-114, september, 1973,
c     bessel functions i and j of complex argument and integer order,
c     by d. j. sookne.
c
c***see also  zbesh, zbesi, zbesk
c***routines called  zabs, zdiv
c***revision history  (yymmdd)
c   830501  date written
c   910415  prologue converted to version 4.0 format.  (bab)
c***end prologue  zrati
      double precision ak, amagz, ap1, ap2, arg, az, cdfnui, cdfnur,
     * conei, coner, cyi, cyr, czeroi, czeror, dfnu, fdnu, flam, fnu,
     * fnup, pti, ptr, p1i, p1r, p2i, p2r, rak, rap1, rho, rt2, rzi,
     * rzr, test, test1, tol, tti, ttr, t1i, t1r, zi, zr, zabs
      integer i, id, idnu, inu, itime, k, kk, magz, n
      dimension cyr(n), cyi(n)
      external zabs
      data czeror,czeroi,coner,conei,rt2/
     1 0.0d0, 0.0d0, 1.0d0, 0.0d0, 1.41421356237309505d0 /
c***first executable statement  zrati
      az = zabs(zr,zi)
      inu = fnu
      idnu = inu + n - 1
      magz = az
      amagz = magz+1
      fdnu = idnu
      fnup = max(amagz,fdnu)
      id = idnu - magz - 1
      itime = 1
      k = 1
      ptr = 1.0d0/az
      rzr = ptr*(zr+zr)*ptr
      rzi = -ptr*(zi+zi)*ptr
      t1r = rzr*fnup
      t1i = rzi*fnup
      p2r = -t1r
      p2i = -t1i
      p1r = coner
      p1i = conei
      t1r = t1r + rzr
      t1i = t1i + rzi
      if (id.gt.0) id = 0
      ap2 = zabs(p2r,p2i)
      ap1 = zabs(p1r,p1i)
c-----------------------------------------------------------------------
c     the overflow test on k(fnu+i-1,z) before the call to cbknu
c     guarantees that p2 is on scale. scale test1 and all subsequent
c     p2 values by ap1 to ensure that an overflow does not occur
c     prematurely.
c-----------------------------------------------------------------------
      arg = (ap2+ap2)/(ap1*tol)
      test1 = sqrt(arg)
      test = test1
      rap1 = 1.0d0/ap1
      p1r = p1r*rap1
      p1i = p1i*rap1
      p2r = p2r*rap1
      p2i = p2i*rap1
      ap2 = ap2*rap1
   10 continue
      k = k + 1
      ap1 = ap2
      ptr = p2r
      pti = p2i
      p2r = p1r - (t1r*ptr-t1i*pti)
      p2i = p1i - (t1r*pti+t1i*ptr)
      p1r = ptr
      p1i = pti
      t1r = t1r + rzr
      t1i = t1i + rzi
      ap2 = zabs(p2r,p2i)
      if (ap1.le.test) go to 10
      if (itime.eq.2) go to 20
      ak = zabs(t1r,t1i)*0.5d0
      flam = ak + sqrt(ak*ak-1.0d0)
      rho = min(ap2/ap1,flam)
      test = test1*sqrt(rho/(rho*rho-1.0d0))
      itime = 2
      go to 10
   20 continue
      kk = k + 1 - id
      ak = kk
      t1r = ak
      t1i = czeroi
      dfnu = fnu + (n-1)
      p1r = 1.0d0/ap2
      p1i = czeroi
      p2r = czeror
      p2i = czeroi
      do 30 i=1,kk
        ptr = p1r
        pti = p1i
        rap1 = dfnu + t1r
        ttr = rzr*rap1
        tti = rzi*rap1
        p1r = (ptr*ttr-pti*tti) + p2r
        p1i = (ptr*tti+pti*ttr) + p2i
        p2r = ptr
        p2i = pti
        t1r = t1r - coner
   30 continue
      if (p1r.ne.czeror .or. p1i.ne.czeroi) go to 40
      p1r = tol
      p1i = tol
   40 continue
      call zdiv(p2r, p2i, p1r, p1i, cyr(n), cyi(n))
      if (n.eq.1) return
      k = n - 1
      ak = k
      t1r = ak
      t1i = czeroi
      cdfnur = fnu*rzr
      cdfnui = fnu*rzi
      do 60 i=2,n
        ptr = cdfnur + (t1r*rzr-t1i*rzi) + cyr(k+1)
        pti = cdfnui + (t1r*rzi+t1i*rzr) + cyi(k+1)
        ak = zabs(ptr,pti)
        if (ak.ne.czeror) go to 50
        ptr = tol
        pti = tol
        ak = tol*rt2
   50   continue
        rak = coner/ak
        cyr(k) = rak*ptr*rak
        cyi(k) = -rak*pti*rak
        t1r = t1r - coner
        k = k - 1
   60 continue
      return
      end
