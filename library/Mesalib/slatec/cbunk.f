*deck cbunk
      subroutine cbunk (z, fnu, kode, mr, n, y, nz, tol, elim, alim)
c***begin prologue  cbunk
c***subsidiary
c***purpose  subsidiary to cbesh and cbesk
c***library   slatec
c***type      all (cbunk-a, zbunk-a)
c***author  amos, d. e., (snl)
c***description
c
c     cbunk computes the k bessel function for fnu.gt.fnul.
c     according to the uniform asymptotic expansion for k(fnu,z)
c     in cunk1 and the expansion for h(2,fnu,z) in cunk2
c
c***see also  cbesh, cbesk
c***routines called  cunk1, cunk2
c***revision history  (yymmdd)
c   830501  date written
c   910415  prologue converted to version 4.0 format.  (bab)
c***end prologue  cbunk
      complex y, z
      real alim, ax, ay, elim, fnu, tol, xx, yy
      integer kode, mr, n, nz
      dimension y(n)
c***first executable statement  cbunk
      nz = 0
      xx = real(z)
      yy = aimag(z)
      ax = abs(xx)*1.7321e0
      ay = abs(yy)
      if (ay.gt.ax) go to 10
c-----------------------------------------------------------------------
c     asymptotic expansion for k(fnu,z) for large fnu applied in
c     -pi/3.le.arg(z).le.pi/3
c-----------------------------------------------------------------------
      call cunk1(z, fnu, kode, mr, n, y, nz, tol, elim, alim)
      go to 20
   10 continue
c-----------------------------------------------------------------------
c     asymptotic expansion for h(2,fnu,z*exp(m*hpi)) for large fnu
c     applied in pi/3.lt.abs(arg(z)).le.pi/2 where m=+i or -i
c     and hpi=pi/2
c-----------------------------------------------------------------------
      call cunk2(z, fnu, kode, mr, n, y, nz, tol, elim, alim)
   20 continue
      return
      end
