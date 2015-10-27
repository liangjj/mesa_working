*deck cacai
      subroutine cacai (z, fnu, kode, mr, n, y, nz, rl, tol, elim, alim)
c***begin prologue  cacai
c***subsidiary
c***purpose  subsidiary to cairy
c***library   slatec
c***type      all (cacai-a, zacai-a)
c***author  amos, d. e., (snl)
c***description
c
c     cacai applies the analytic continuation formula
c
c         k(fnu,zn*exp(mp))=k(fnu,zn)*exp(-mp*fnu) - mp*i(fnu,zn)
c                 mp=pi*mr*cmplx(0.0,1.0)
c
c     to continue the k function from the right half to the left
c     half z plane for use with cairy where fnu=1/3 or 2/3 and n=1.
c     cacai is the same as cacon with the parts for larger orders and
c     recurrence removed. a recursive call to cacon can result if cacon
c     is called from cairy.
c
c***see also  cairy
c***routines called  casyi, cbknu, cmlri, cs1s2, cseri, r1mach
c***revision history  (yymmdd)
c   830501  date written
c   910415  prologue converted to version 4.0 format.  (bab)
c***end prologue  cacai
      complex csgn, cspn, c1, c2, y, z, zn, cy
      real alim, arg, ascle, az, cpn, dfnu, elim, fmr, fnu, pi, rl,
     * sgn, spn, tol, yy, r1mach
      integer inu, iuf, kode, mr, n, nn, nw, nz
      dimension y(n), cy(2)
      data pi / 3.14159265358979324e0 /
c***first executable statement  cacai
      nz = 0
      zn = -z
      az = abs(z)
      nn = n
      dfnu = fnu + (n-1)
      if (az.le.2.0e0) go to 10
      if (az*az*0.25e0.gt.dfnu+1.0e0) go to 20
   10 continue
c-----------------------------------------------------------------------
c     power series for the i function
c-----------------------------------------------------------------------
      call cseri(zn, fnu, kode, nn, y, nw, tol, elim, alim)
      go to 40
   20 continue
      if (az.lt.rl) go to 30
c-----------------------------------------------------------------------
c     asymptotic expansion for large z for the i function
c-----------------------------------------------------------------------
      call casyi(zn, fnu, kode, nn, y, nw, rl, tol, elim, alim)
      if (nw.lt.0) go to 70
      go to 40
   30 continue
c-----------------------------------------------------------------------
c     miller algorithm normalized by the series for the i function
c-----------------------------------------------------------------------
      call cmlri(zn, fnu, kode, nn, y, nw, tol)
      if(nw.lt.0) go to 70
   40 continue
c-----------------------------------------------------------------------
c     analytic continuation to the left half plane for the k function
c-----------------------------------------------------------------------
      call cbknu(zn, fnu, kode, 1, cy, nw, tol, elim, alim)
      if (nw.ne.0) go to 70
      fmr = mr
      sgn = -sign(pi,fmr)
      csgn = cmplx(0.0e0,sgn)
      if (kode.eq.1) go to 50
      yy = -aimag(zn)
      cpn = cos(yy)
      spn = sin(yy)
      csgn = csgn*cmplx(cpn,spn)
   50 continue
c-----------------------------------------------------------------------
c     calculate cspn=exp(fnu*pi*i) to minimize losses of significance
c     when fnu is large
c-----------------------------------------------------------------------
      inu = fnu
      arg = (fnu-inu)*sgn
      cpn = cos(arg)
      spn = sin(arg)
      cspn = cmplx(cpn,spn)
      if (mod(inu,2).eq.1) cspn = -cspn
      c1 = cy(1)
      c2 = y(1)
      if (kode.eq.1) go to 60
      iuf = 0
      ascle = 1.0e+3*r1mach(1)/tol
      call cs1s2(zn, c1, c2, nw, ascle, alim, iuf)
      nz = nz + nw
   60 continue
      y(1) = cspn*c1 + csgn*c2
      return
   70 continue
      nz = -1
      if(nw.eq.(-2)) nz=-2
      return
      end
