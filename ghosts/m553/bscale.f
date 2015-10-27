      subroutine bscale(v,n,xx)
      implicit real*8(a-h,o-z)
      dimension v(n)
c
c  scale a vector so the largest element has abs. val. one
c
      xx=0.d0
      do 10 i=1,n
      if(dabs(v(i)).lt.xx)go to 10
       xx=dabs(v(i))
  10  continue
c
      xfac=1.d0/xx
c
      call vscale(v,xfac,n)
c
      return
      end
