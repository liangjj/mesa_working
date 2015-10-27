*deck @(#)rtocm.f	1.1 9/8/91
c***begin prologue     rtocm
c***date written       xxxxxx   (yymmdd)
c***revision date      890427   (yymmdd)
c***keywords           m7000, bessel functions
c***author             schneider, barry (nsf)
c***source             m7000
c***purpose            copy real to complex array
c***routines called    iosys, util and mdutil
c***end prologue       rtocm
      subroutine rtocm(ar,ac,n)
      implicit integer (a-z)
      real *8  ar
      complex *16 ac
      dimension ar(n), ac(n)
      do 10 i=1,n
         ac(i)=ar(i)   
   10 continue     
      return
      end

