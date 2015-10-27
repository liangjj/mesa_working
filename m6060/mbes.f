*deck @(#)mbes.f
c***begin prologue     mbes
c***date written       xxxxxx   (yymmdd)
c***revision date      890427   (yymmdd)
c***keywords           m6060, link 6060, bessel functions
c***author             rescigno, tom (llnl)
c***source             m6060
c***purpose            bessel functions from spline fit
c***description        calculation of bessel functions from spline
c***                   no derivatives
c***routines called    iosys, util and mdutil
c***end prologue       mbes
      subroutine mbes(hs,cj,hp,grid,rvec,x,krvec,kchan,
     1                 rmin,rdel,rd26,ipow,nlm,lch,npt,nr,
     2                 lmax,maxlm,nchan,dimlm,dimc)
      implicit integer (a-z)
      real *8 grid, rvec, x, krvec, kchan, rmin, rdel, rd26, ksq
      real *8 a, b
      complex *16 hs, cj, hp
      dimension grid(4,npt), rvec(npt), krvec(npt), kchan(nchan)
      dimension hs(nr,0:lmax), cj(nr,0:lmax), hp(npt,maxlm,nchan)
      dimension ipow(0:lmax), x(nr), nlm(dimc), lch(dimlm,dimc)
      do 10 i = 1,npt
         rvec(i)=sqrt(grid(1,i)**2 + grid(2,i)**2 +grid(3,i)**2)
   10 continue
      do 20 ch1=1,nchan
         do 30 i=1,npt
            krvec(i)=kchan(ch1)*rvec(i)
   30    continue
         ksq = sqrt(kchan(ch1))
         k52 = ksq*kchan(ch1)**2
         do 40 lv=1,nlm(ch1)
               l=lch(lv,ch1)
            do 50 i=1,npt
c********************
c the expression for klo is split into two fortran
c statements so that klo is never less than 1
c*****************
               klo=(krvec(i)-rmin)/rdel
               klo=klo+1
               a=(x(klo+1)-krvec(i))/rdel
               b=(krvec(i)-x(klo))/rdel
               hp(i,lv,ch1)=(a*hs(klo,l)+b*hs(klo+1,l)+(a*(a*a-1.)*
     1                    cj(klo,l)+b*(b*b-1.)*cj(klo+1,l))*rd26)*ksq
   50       continue
   40    continue
   20 continue
      return
      end
