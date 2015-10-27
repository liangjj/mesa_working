*deck zacai
      subroutine zacai (zr, zi, fnu, kode, mr, n, yr, yi, nz, rl, tol,
     +   elim, alim)
c***begin prologue  zacai
c***subsidiary
c***purpose  subsidiary to zairy
c***library   slatec
c***type      all (cacai-a, zacai-a)
c***author  amos, d. e., (snl)
c***description
c
c     zacai applies the analytic continuation formula
c
c         k(fnu,zn*exp(mp))=k(fnu,zn)*exp(-mp*fnu) - mp*i(fnu,zn)
c                 mp=pi*mr*cmplx(0.0,1.0)
c
c     to continue the k function from the right half to the left
c     half z plane for use with zairy where fnu=1/3 or 2/3 and n=1.
c     zacai is the same as zacon with the parts for larger orders and
c     recurrence removed. a recursive call to zacon can result if zacon
c     is called from zairy.
c
c***see also  zairy
c***routines called  d1mach, zabs, zasyi, zbknu, zmlri, zs1s2, zseri
c***revision history  (yymmdd)
c   830501  date written
c   910415  prologue converted to version 4.0 format.  (bab)
c***end prologue  zacai
c     complex csgn,cspn,c1,c2,y,z,zn,cy
      double precision alim, arg, ascle, az, csgnr, csgni, cspnr,
     * cspni, c1r, c1i, c2r, c2i, cyr, cyi, dfnu, elim, fmr, fnu, pi,
     * rl, sgn, tol, yy, yr, yi, zr, zi, znr, zni, d1mach, zabs
      integer inu, iuf, kode, mr, n, nn, nw, nz
      dimension yr(n), yi(n), cyr(2), cyi(2)
      external zabs
      data pi / 3.14159265358979324d0 /
c***first executable statement  zacai
      nz = 0
      znr = -zr
      zni = -zi
      az = zabs(zr,zi)
      nn = n
      dfnu = fnu + (n-1)
      if (az.le.2.0d0) go to 10
      if (az*az*0.25d0.gt.dfnu+1.0d0) go to 20
   10 continue
c-----------------------------------------------------------------------
c     power series for the i function
c-----------------------------------------------------------------------
      call zseri(znr, zni, fnu, kode, nn, yr, yi, nw, tol, elim, alim)
      go to 40
   20 continue
      if (az.lt.rl) go to 30
c-----------------------------------------------------------------------
c     asymptotic expansion for large z for the i function
c-----------------------------------------------------------------------
      call zasyi(znr, zni, fnu, kode, nn, yr, yi, nw, rl, tol, elim,
     * alim)
      if (nw.lt.0) go to 80
      go to 40
   30 continue
c-----------------------------------------------------------------------
c     miller algorithm normalized by the series for the i function
c-----------------------------------------------------------------------
      call zmlri(znr, zni, fnu, kode, nn, yr, yi, nw, tol)
      if(nw.lt.0) go to 80
   40 continue
c-----------------------------------------------------------------------
c     analytic continuation to the left half plane for the k function
c-----------------------------------------------------------------------
      call zbknu(znr, zni, fnu, kode, 1, cyr, cyi, nw, tol, elim, alim)
      if (nw.ne.0) go to 80
      fmr = mr
      sgn = -dsign(pi,fmr)
      csgnr = 0.0d0
      csgni = sgn
      if (kode.eq.1) go to 50
      yy = -zni
      csgnr = -csgni*sin(yy)
      csgni = csgni*cos(yy)
   50 continue
c-----------------------------------------------------------------------
c     calculate cspn=exp(fnu*pi*i) to minimize losses of significance
c     when fnu is large
c-----------------------------------------------------------------------
      inu = fnu
      arg = (fnu-inu)*sgn
      cspnr = cos(arg)
      cspni = sin(arg)
      if (mod(inu,2).eq.0) go to 60
      cspnr = -cspnr
      cspni = -cspni
   60 continue
      c1r = cyr(1)
      c1i = cyi(1)
      c2r = yr(1)
      c2i = yi(1)
      if (kode.eq.1) go to 70
      iuf = 0
      ascle = 1.0d+3*d1mach(1)/tol
      call zs1s2(znr, zni, c1r, c1i, c2r, c2i, nw, ascle, alim, iuf)
      nz = nz + nw
   70 continue
      yr(1) = cspnr*c1r - cspni*c1i + csgnr*c2r - csgni*c2i
      yi(1) = cspnr*c1i + cspni*c1r + csgnr*c2i + csgni*c2r
      return
   80 continue
      nz = -1
      if(nw.eq.(-2)) nz=-2
      return
      end
