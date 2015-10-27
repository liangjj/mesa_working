*deck max
      subroutine max(x,nd,xmax)
      real *8 x, xmax
      dimension x(1)
      xmax=0.d0
      do 1 i=1,nd
      if(xmax-abs(x(i))) 2,1,1
    2 xmax=abs(x(i))
    1 continue
      return
      end
