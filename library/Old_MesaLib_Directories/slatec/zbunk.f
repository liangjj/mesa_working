*deck zbunk
      subroutine zbunk (zr, zi, fnu, kode, mr, n, yr, yi, nz, tol, elim,
     +   alim)
c***begin prologue  zbunk
c***subsidiary
c***purpose  subsidiary to zbesh and zbesk
c***library   slatec
c***type      all (cbuni-a, zbuni-a)
c***author  amos, d. e., (snl)
c***description
c
c     zbunk computes the k bessel function for fnu.gt.fnul.
c     according to the uniform asymptotic expansion for k(fnu,z)
c     in zunk1 and the expansion for h(2,fnu,z) in zunk2
c
c***see also  zbesh, zbesk
c***routines called  zunk1, zunk2
c***revision history  (yymmdd)
c   830501  date written
c   910415  prologue converted to version 4.0 format.  (bab)
c***end prologue  zbunk
c     complex y,z
      double precision alim, ax, ay, elim, fnu, tol, yi, yr, zi, zr
      integer kode, mr, n, nz
      dimension yr(n), yi(n)
c***first executable statement  zbunk
      nz = 0
      ax = abs(zr)*1.7321d0
      ay = abs(zi)
      if (ay.gt.ax) go to 10
c-----------------------------------------------------------------------
c     asymptotic expansion for k(fnu,z) for large fnu applied in
c     -pi/3.le.arg(z).le.pi/3
c-----------------------------------------------------------------------
      call zunk1(zr, zi, fnu, kode, mr, n, yr, yi, nz, tol, elim, alim)
      go to 20
   10 continue
c-----------------------------------------------------------------------
c     asymptotic expansion for h(2,fnu,z*exp(m*hpi)) for large fnu
c     applied in pi/3.lt.abs(arg(z)).le.pi/2 where m=+i or -i
c     and hpi=pi/2
c-----------------------------------------------------------------------
      call zunk2(zr, zi, fnu, kode, mr, n, yr, yi, nz, tol, elim, alim)
   20 continue
      return
      end
