*deck @(#)tstzro.f	1.2  10/27/94
c***begin prologue     tstzro
c***date written       930922   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           tstzro, link m6203
c***author             schneider, barry (nsf)
c***source             m6203
c***purpose            
c***references         none
c
c***routines called
c***end prologue       tstzro
      subroutine tstzro (r,n)
      implicit integer (a-z)
      real*8 r
      dimension r(n)
      do 10 i=1,n
         if(r(i).lt.1.d-10) then
            r(i)=1.d-12
         endif   
   10 continue
      return
      end

