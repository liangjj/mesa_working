*deck @(#)mkbes.f	1.1 9/8/91
c***begin prologue     mkbes
c***date written       xxxxxx   (yymmdd)
c***revision date      890427   (yymmdd)
c***keywords           m6005, link 6005, bessel functions
c***author             rescigno, tom (llnl)
c***source             m6005
c***purpose            bessel functions from spline fit
c***description        calculation of bessel functions from spline
c***                   fits in link m6004
c***routines called    iosys, util and mdutil
c***end prologue       mkbes
      subroutine mkbes(hs,hsder,cj,cy,hp,hd,rvec,x,krvec,kchan,
     1                         rmin,rdel,rd26,l,npt,nr)
      implicit integer (a-z)
      real *8  rvec, x, krvec, kchan, rmin, rdel, rd26, ksq, k52
      real *8 a, b, hsder, cy, hd
      complex *16 hs, cj, hp
      dimension rvec(npt), krvec(npt), hs(nr), hsder(nr), cj(nr)
      dimension cy(nr), hp(npt), hd(npt), x(nr)
      common /io/ inp, iout
c----------------------------------------------------------------------c
c             this is set up to compute the free function              c
c                              and                                     c
c             the effect of ( H0 - E ) on the free function.           c
c             this is done in order to make the integrand vanish       c
c             at large values of r.                                    c
c----------------------------------------------------------------------c  
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
         hd(i)=a*hsder(klo)+b*hsder(klo+1)+(a*(a*a-1.)*
     1                     cy(klo)+b*(b*b-1.)*cy(klo+1))*rd26
         hd(i)=hd(i)*k52/krvec(i)**(1-ipow)
   20 continue
      return
      end

