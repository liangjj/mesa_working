*deck %W%  %G%
      subroutine astv(v,xx,n)
      implicit real*8(a-h,o-z)
      dimension v(n)
c
      do 10 i=1,n
        v(i)=v(i)+xx
  10  continue
c
      return
      end
