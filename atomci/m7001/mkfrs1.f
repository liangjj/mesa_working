*deck @(#)mkfrs1.f	1.1 9/8/91
c***begin prologue     mkfrs1
c***date written       920417   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m7001, bessel functions
c***author             schneider, barry (nsf)
c***source             m7000
c***purpose            real bessel function from universal complex function
c***routines called    iosys, util and mdutil
c***end prologue       mkfrs1
      subroutine mkfrs1(frec,frer,npt,sze)
      implicit integer (a-z)
      real *8  frer
      complex *16 frec
      dimension frec(npt,sze), frer(npt,sze)
      do 10 i=1,sze
         do 20 j=1,npt
            frer(j,i) = imag(frec(j,i))
   20    continue  
   10 continue     
      return
      end

