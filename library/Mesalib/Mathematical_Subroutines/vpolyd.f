*deck @(#)vpolyd.f	5.1  11/6/94
      subroutine vpolyd(t,x,c,nc,n)
c***begin prologue     vpolyd.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             unknown
c***source             @(#)vpolyd.f	5.1   11/6/94
c***purpose            
c***description
c    to evaluate the -(nc-1)th order polynomial t(x)
c    using the coefficients c.  c(1) is the coefficient of
c    x**-(nc-1); c(nc) is the constant in the polynomial.
c     
c***references
c
c***routines called
c
c***end prologue       vpolyd.f
c
c
      implicit integer (a-z)
c
      character*4 itoc
      real*8 t(n),x(n),c(nc)
c
      if (nc.lt.1) then
         call lnkerr('vpolyd: error in polynomial degree:'//itoc(nc))
      end if
c
c
      do 2 i=1,n
         t(i)=c(1)/x(i)
    2 continue
      do 4 coef=2,nc-1
         do 3 i=1,n
            t(i)=(t(i)+c(coef))/x(i)
    3    continue
    4 continue
      do 5 i=1,n
         t(i)=t(i)+c(nc)
    5 continue
c
c
      return
      end
