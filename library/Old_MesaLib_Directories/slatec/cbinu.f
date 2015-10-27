*deck cbinu
      subroutine cbinu (z, fnu, kode, n, cy, nz, rl, fnul, tol, elim,
     +   alim)
c***begin prologue  cbinu
c***subsidiary
c***purpose  subsidiary to cairy, cbesh, cbesi, cbesj, cbesk and cbiry
c***library   slatec
c***type      all (cbinu-a, zbinu-a)
c***author  amos, d. e., (snl)
c***description
c
c     cbinu computes the i function in the right half z plane
c
c***see also  cairy, cbesh, cbesi, cbesj, cbesk, cbiry
c***routines called  casyi, cbuni, cmlri, cseri, cuoik, cwrsk
c***revision history  (yymmdd)
c   830501  date written
c   910415  prologue converted to version 4.0 format.  (bab)
c***end prologue  cbinu
      complex cw, cy, czero, z
      real alim, az, dfnu, elim, fnu, fnul, rl, tol
      integer i, inw, kode, n, nlast, nn, nui, nw, nz
      dimension cy(n), cw(2)
      data czero / (0.0e0,0.0e0) /
c***first executable statement  cbinu
      nz = 0
      az = abs(z)
      nn = n
      dfnu = fnu + (n-1)
      if (az.le.2.0e0) go to 10
      if (az*az*0.25e0.gt.dfnu+1.0e0) go to 20
   10 continue
c-----------------------------------------------------------------------
c     power series
c-----------------------------------------------------------------------
      call cseri(z, fnu, kode, nn, cy, nw, tol, elim, alim)
      inw = abs(nw)
      nz = nz + inw
      nn = nn - inw
      if (nn.eq.0) return
      if (nw.ge.0) go to 120
      dfnu = fnu + (nn-1)
   20 continue
      if (az.lt.rl) go to 40
      if (dfnu.le.1.0e0) go to 30
      if (az+az.lt.dfnu*dfnu) go to 50
c-----------------------------------------------------------------------
c     asymptotic expansion for large z
c-----------------------------------------------------------------------
   30 continue
      call casyi(z, fnu, kode, nn, cy, nw, rl, tol, elim, alim)
      if (nw.lt.0) go to 130
      go to 120
   40 continue
      if (dfnu.le.1.0e0) go to 70
   50 continue
c-----------------------------------------------------------------------
c     overflow and underflow test on i sequence for miller algorithm
c-----------------------------------------------------------------------
      call cuoik(z, fnu, kode, 1, nn, cy, nw, tol, elim, alim)
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
      call cmlri(z, fnu, kode, nn, cy, nw, tol)
      if(nw.lt.0) go to 130
      go to 120
   80 continue
c-----------------------------------------------------------------------
c     miller algorithm normalized by the wronskian
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     overflow test on k functions used in wronskian
c-----------------------------------------------------------------------
      call cuoik(z, fnu, kode, 2, 2, cw, nw, tol, elim, alim)
      if (nw.ge.0) go to 100
      nz = nn
      do 90 i=1,nn
        cy(i) = czero
   90 continue
      return
  100 continue
      if (nw.gt.0) go to 130
      call cwrsk(z, fnu, kode, nn, cy, nw, cw, tol, elim, alim)
      if (nw.lt.0) go to 130
      go to 120
  110 continue
c-----------------------------------------------------------------------
c     increment fnu+nn-1 up to fnul, compute and recur backward
c-----------------------------------------------------------------------
      nui = fnul-dfnu + 1
      nui = max(nui,0)
      call cbuni(z, fnu, kode, nn, cy, nw, nui, nlast, fnul, tol, elim,
     * alim)
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
