*deck stnzro.f
c***begin prologue     stnzro
c***date written       930922   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           stnzro, link m6201
c***author             schneider, barry (nsf)
c***source             m6201
c***purpose            
c***references         none
c
c***routines called
c***end prologue       stnzro
      subroutine stnzro (f,n)
      implicit integer (a-z)
      real*8 f
      dimension f(n)
      do 10 i=1,n
         if( f(i).le.1.d-10 ) then
             f(i)=1.d-10
         endif   
   10 continue
      return
      end

