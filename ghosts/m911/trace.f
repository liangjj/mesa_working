*deck @(#)trace.f	4.1  7/7/93
      subroutine trace(dm,n)
      implicit real*8(a-h,o-z)
      real*8 dm(n,n)
c
c
      trac=0.d0
      do 10 i=1,n
         trac=trac+dm(i,i)
  10  continue
c
c     write(iout,20) trac
  20  format(/,'  trace of square 1-e density = ',f10.7)
c
c
      return
      end
