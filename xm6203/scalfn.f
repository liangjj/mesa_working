*deck @(#)scalfn.f	1.2  10/27/94
c***begin prologue     scalfn
c***date written       921221   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           scalfn, link m6203, legendre integrals
c***author             schneider, barry (nsf)
c***source             m6203
c***purpose            scale a matrix of function values
c***references         none
c
c***routines called
c***end prologue       scalfn
      subroutine scalfn (f,wt,n,l,m)
      implicit integer (a-z)
      real*8 f, wt
      dimension f(n,m:l), wt(n)
      common/io/inp,iout
      do 10 i=m,l
         do 20 j=1,n
            f(j,i) = wt(j)*f(j,i)
   20    continue
   10 continue 
      return
      end

