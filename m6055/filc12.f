*deck @(#)filc12.f	1.1 9/8/91
c***begin prologue     filc12
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
c***end prologue       filc12
      subroutine filc12(cc,cr,a,b,h1,h2,nc,l,m,maxl,npt,lmax,dim)
      implicit integer (a-z)
      complex*16 cc, h1
      real *8 a, b, h2, cr
      dimension cc(nc,npt), a(npt,0:lmax,0:2*lmax), b(4,npt)
      dimension h1(npt,maxl), h2(npt,maxl), l(dim), m(dim)
      dimension cr(nc,npt)
      do 10 nolm1=1,nc
         l1=l(nolm1)
         m1=m(nolm1)
         do 20 grpt=1,npt
            cc(nolm1,grpt)=a(grpt,l1,m1)*b(4,grpt)*h1(grpt,nolm1)
            cr(nolm1,grpt)=a(grpt,l1,m1)*b(4,grpt)*h2(grpt,nolm1)
   20    continue
   10 continue
      return
      end
