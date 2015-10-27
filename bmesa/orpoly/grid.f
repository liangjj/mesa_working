*deck grid.f
c***begin prologue     grid
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             math
c***purpose            
c***                   
c***                   
c
c***references         
c
c***routines called    
c***end prologue       grid
      subroutine grid(x,npts)
c
      implicit integer (a-z)
      real*8 x, stp
      common/io/inp, iout 
      dimension x(npts)
      stp=2.d0/(npts-1)
      x(1)=-1.d0
      do 10 i=2,npts
         x(i)=x(i-1)+stp
 10   continue
      return
      end       
