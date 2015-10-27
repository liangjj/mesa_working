*deck mke.f
c***begin prologue     mke
c***date written       930623   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           energy grid
c***author             schneider, barry (nsf)
c***source             
c***purpose            make a grid                     
c***                   
c
c***references         
c
c***routines called    
c***end prologue      mke
      subroutine mke(e,e0,dele,n)
c
      implicit integer (a-z)
      real*8  e, e0, dele
      dimension e(n)
c
      e(1)=e0
      do 10 i=2,n
         e(i)=e(i-1)+dele
   10 continue     
      return
      end
