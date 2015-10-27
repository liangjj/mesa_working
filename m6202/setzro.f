*deck setzro.f
c***begin prologue     setzro
c***date written       930922   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           setzro, link m6201
c***author             schneider, barry (nsf)
c***source             m6201
c***purpose            
c***references         none
c
c***routines called
c***end prologue       setzro
      subroutine setzro (f,n)
      implicit integer (a-z)
      real*8 f
      dimension f(n)
      do 10 i=1,n
         if(abs(f(i)).lt.1.d-15) then
            f(i)=0.d0
         endif   
   10 continue
      return
      end

