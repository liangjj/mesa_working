*deck @(#)mkcoul.f	1.1 9/8/91
c***begin prologue     mkcoul
c***date written       xxxxxx   (yymmdd)
c***revision date      890427   (yymmdd)
c***keywords           m6005, link 6005, bessel functions
c***author             rescigno, tom (llnl)
c***source             m6005
c***purpose            coulomb functions from spline fit
c***description        calculation of coulomb functions from spline
c***                   fits in link m6004
c***routines called    iosys, util and mdutil
c***end prologue       mkcoul
      subroutine mkcoul(hs,cj,hp,rvec,x,krvec,kchan,rmin,rdel,
     1                  rd26,l,npt,nr)
      implicit integer (a-z)
      real *8  rvec, x, krvec, kchan, rmin, rdel, rd26, ksq, k52
      real *8 a, b, hs, cj, hp
      dimension rvec(npt), krvec(npt), hs(nr), cj(nr), hp(npt), x(nr)
      common /io/ inp, iout
      ipow=1
      if (l.eq.0) then
          ipow=0
      endif
      do 10 i=1,npt
         krvec(i)=kchan*rvec(i)
   10 continue
      ksq = sqrt(kchan)
      k52 = ksq*kchan**2
      do 20 i=1,npt
c********************
c the expression for klo is split into two fortran
c statements so that klo is never less than 1
c*****************
         klo=(krvec(i)-rmin)/rdel
         klo=klo+1
         a=(x(klo+1)-krvec(i))/rdel
         b=(krvec(i)-x(klo))/rdel
         hp(i)=(a*hs(klo)+b*hs(klo+1)+(a*(a*a-1.)*
     1                    cj(klo)+b*(b*b-1.)*cj(klo+1))*rd26)*ksq
   20 continue
      return
      end

