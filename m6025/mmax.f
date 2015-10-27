      subroutine mmax(x,nd,xmax)
      implicit real*8 ( a - h, o - z )
      dimension x(1)
      xmax=0.d0
      do 1 i=1,nd
      if(xmax-abs(x(i))) 2,1,1
    2 xmax=abs(x(i))
    1 continue
      return
      end
