c----------------------------------------------------------------------c
c                lots of utility routines to handle                    c
c                complex or mixture of real and complex                c
c                         matrices                                     c
c----------------------------------------------------------------------c
*dk csxpy
      subroutine csxpy (n,sa,sx,sy)
      real*8 sx
      complex*16 sa, sy
      dimension sx(n), sy(n)
      do 10 i=1,n
   10 sy(i)=sy(i)+sa*sx(i)
      return
      end
