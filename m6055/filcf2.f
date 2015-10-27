*deck @(#)filcf2.f	1.1 9/8/91
c***begin prologue     filcf2
c***date written       880423   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           utility, matrix
c***author             schneider, barry (lanl)
c***source             m6005
c***purpose            fill matrix in packed form
c***description
c***
c***
c***references         none
c
c***routines called
c***end prologue       filcf2
      subroutine filcf2(c,a,b,v,nc,l,m,maxl,npt,lmax,dim)
      implicit integer (a-z)
      complex *16 b, c
      real *8 a, v
      dimension c(nc,npt), a(npt,0:lmax,0:2*lmax), b(npt,maxl)
      dimension l(dim), m(dim), v(npt)
      do 10 nolm1=1,nc
         l1=l(nolm1)
         m1=m(nolm1)
         do 20 grpt=1,npt
            c(nolm1,grpt)=a(grpt,l1,m1)*v(grpt)*b(grpt,nolm1)
   20    continue
   10 continue
      return
      end
