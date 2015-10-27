*deck %W%  %G%
      subroutine brzero(v,n)
      implicit real*8(a-h,o-z)
      dimension v(n)
c
      do 10 i=1,n
        v(i)=0.d0
  10  continue
c
      return
      end
