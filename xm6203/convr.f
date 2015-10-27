*deck @(#)convr.f	1.2  10/27/94
c***begin prologue     convr
c***date written       930922   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           convr, link m6203, legendre quadrature
c***author             schneider, barry (nsf)
c***source             m6203
c***purpose            convert legendre quadrature from (-1,1) to (a,b)
c***references         none
c
c***routines called
c***end prologue       convr
      subroutine convr (a,b,r,wtr,npt)
      implicit integer (a-z)
      real*8 a, b, r, wtr, afac, bfac
      dimension r(npt), wtr(npt)
      afac=(b-a)*.5d0
      bfac=(b+a)*.5d0      
      do 10 i=1,npt
         r(i)=afac*r(i)+bfac
         wtr(i)=afac*wtr(i)
   10 continue 
      return
      end

