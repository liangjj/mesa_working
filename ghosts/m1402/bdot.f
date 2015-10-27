*deck %W%  %G%
      subroutine bdot(v,w,xx,n)
      implicit real*8(a-h,o-z)
      dimension v(n,n),w(n,n)
c
      xx=0.d0
      do 20 i=1,n
      do 10 j=1,n
        xx=xx+v(j,i)*w(i,j)
  10  continue
 20   continue
c
      return
      end
