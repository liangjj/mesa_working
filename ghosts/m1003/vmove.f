*deck %W%  %G%
      subroutine vmove(v,w,n)
      implicit real*8(a-h,o-z)
      dimension v(n),w(n)
c
      do 10 i=1,n
        v(i)=w(i)
  10  continue
c
      return
      end
