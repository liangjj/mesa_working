*deck %W%  %G%
      subroutine ssdot(v,w,xx,n)
      implicit real*8(a-h,o-z)
      dimension v(n),w(n)
c
      xx=0.d0
      do 10 i=1,n
        xx=xx+v(i)*w(i)
  10  continue
c
      return
      end
