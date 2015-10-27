*deck zbinu
      subroutine zbinu (zr, zi, fnu, kode, n, cyr, cyi, nz, rl, fnul,
     +   tol, elim, alim)
c***begin prologue  zbinu
c***subsidiary
c***purpose  subsidiary to zairy, zbesh, zbesi, zbesj, zbesk and zbiry
c***library   slatec
c***type      all (cbinu-a, zbinu-a)
c***author  amos, d. e., (snl)
c***description
c
c     zbinu computes the i function in the right half z plane
c
c***see also  zairy, zbesh, zbesi, zbesj, zbesk, zbiry
c***routines called  zabs, zasyi, zbuni, zmlri, zseri, zuoik, zwrsk
c***revision history  (yymmdd)
c   830501  date written
c   910415  prologue converted to version 4.0 format.  (bab)
c***end prologue  zbinu
      double precision alim, az, cwi, cwr, cyi, cyr, dfnu, elim, fnu,
     * fnul, rl, tol, zeroi, zeror, zi, zr, zabs
      integer i, inw, kode, n, nlast, nn, nui, nw, nz
      dimension cyr(n), cyi(n), cwr(2), cwi(2)
      external zabs
      data zeror,zeroi / 0.0d0, 0.0d0 /
c***first executable statement  zbinu
      nz = 0
      az = zabs(zr,zi)
      nn = n
      dfnu = fnu + (n-1)
      if (az.le.2.0d0) go to 10
      if (az*az*0.25d0.gt.dfnu+1.0d0) go to 20
   10 continue
c-----------------------------------------------------------------------
c     power series
c-----------------------------------------------------------------------
      call zseri(zr, zi, fnu, kode, nn, cyr, cyi, nw, tol, elim, alim)
      inw = abs(nw)
      nz = nz + inw
      nn = nn - inw
      if (nn.eq.0) return
      if (nw.ge.0) go to 120
      dfnu = fnu + (nn-1)
   20 continue
      if (az.lt.rl) go to 40
      if (dfnu.le.1.0d0) go to 30
      if (az+az.lt.dfnu*dfnu) go to 50
c-----------------------------------------------------------------------
c     asymptotic expansion for large z
c-----------------------------------------------------------------------
   30 continue
      call zasyi(zr, zi, fnu, kode, nn, cyr, cyi, nw, rl, tol, elim,
     * alim)
      if (nw.lt.0) go to 130
      go to 120
   40 continue
      if (dfnu.le.1.0d0) go to 70
   50 continue
c-----------------------------------------------------------------------
c     overflow and underflow test on i sequence for miller algorithm
c-----------------------------------------------------------------------
      call zuoik(zr, zi, fnu, kode, 1, nn, cyr, cyi, nw, tol, elim,
     * alim)
      if (nw.lt.0) go to 130
      nz = nz + nw
      nn = nn - nw
      if (nn.eq.0) return
      dfnu = fnu+(nn-1)
      if (dfnu.gt.fnul) go to 110
      if (az.gt.fnul) go to 110
   60 continue
      if (az.gt.rl) go to 80
   70 continue
c-----------------------------------------------------------------------
c     miller algorithm normalized by the series
c-----------------------------------------------------------------------
      call zmlri(zr, zi, fnu, kode, nn, cyr, cyi, nw, tol)
      if(nw.lt.0) go to 130
      go to 120
   80 continue
c-----------------------------------------------------------------------
c     miller algorithm normalized by the wronskian
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     overflow test on k functions used in wronskian
c-----------------------------------------------------------------------
      call zuoik(zr, zi, fnu, kode, 2, 2, cwr, cwi, nw, tol, elim,
     * alim)
      if (nw.ge.0) go to 100
      nz = nn
      do 90 i=1,nn
        cyr(i) = zeror
        cyi(i) = zeroi
   90 continue
      return
  100 continue
      if (nw.gt.0) go to 130
      call zwrsk(zr, zi, fnu, kode, nn, cyr, cyi, nw, cwr, cwi, tol,
     * elim, alim)
      if (nw.lt.0) go to 130
      go to 120
  110 continue
c-----------------------------------------------------------------------
c     increment fnu+nn-1 up to fnul, compute and recur backward
c-----------------------------------------------------------------------
      nui = fnul-dfnu + 1
      nui = max(nui,0)
      call zbuni(zr, zi, fnu, kode, nn, cyr, cyi, nw, nui, nlast, fnul,
     * tol, elim, alim)
      if (nw.lt.0) go to 130
      nz = nz + nw
      if (nlast.eq.0) go to 120
      nn = nlast
      go to 60
  120 continue
      return
  130 continue
      nz = -1
      if(nw.eq.(-2)) nz=-2
      return
      end
