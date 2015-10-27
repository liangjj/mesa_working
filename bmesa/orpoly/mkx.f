*deck mkx.f
c***begin prologue     mkx
c***date written       930623   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           root finder
c***author             schneider, barry (nsf)
c***source             m6202
c***purpose            make a grid                     
c***                   
c
c***references         
c
c***routines called    
c***end prologue      mkx
      subroutine mkx(x,rl,del,n)
c
      implicit integer (a-z)
      real*8  del, rl, x
      dimension x(n)
c
      x(1)=rl
      do 10 i=2,n
         x(i)=x(i-1)+del
   10 continue     
      return
      end
