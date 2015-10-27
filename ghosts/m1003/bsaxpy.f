*deck %W%  %G%
      subroutine bsaxpy(v,w,xx,n)
      implicit real*8(a-h,o-z)
      dimension v(n),w(n)
c
      do 10 i=1,n
        v(i)=v(i)+xx*w(i)
  10  continue
c
      return
      end
