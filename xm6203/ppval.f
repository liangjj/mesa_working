*deck @(#)ppval.f	1.2  10/27/94
      subroutine ppval (x,y,break,coef,n,nbreak,order,ind,jderiv)
calculates value at  x  of  jderiv-th derivative of spline from pp-repr
      implicit real *8 (a-h,o-z)
      integer order
      common /io/ inp, iout
      save
      dimension break(nbreak), coef(order,nbreak), x(n), y(n), ind(n)
      call rzero(y,n)
      fmmjdr = order - jderiv
c  derivatives of order  k  or higher are identically zero.
      if (fmmjdr .gt. 0.d0) then
c  index of largest breakpoint to the left of  x is in ind.
c  evaluate  jderiv -th derivative of  i -th polynomial piece at  x .
          lim=order-jderiv
          do 10 mm=1,lim
             m=order-mm+1
             do 20 i=1,n
                h = x(i) - break(ind(i))
                y(i) = (y(i)/fmmjdr)*h + coef(m,ind(i))
   20        continue
             fmmjdr = fmmjdr - 1.d0
   10     continue     
      endif
      return
      end
