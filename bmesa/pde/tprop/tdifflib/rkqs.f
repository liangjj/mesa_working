*deck rkqs
      subroutine rkqs(y,dy,yscal,yerr,ytemp,temp,x,htry,eps,
     1                hdid,hnext,n)
      implicit integer(a-z)
      real*8 y, dy, yscal, eps, hdid, hnext, htry, x
      real*8 errmax, h, htemp, xnew, yerr, ytemp, temp
      real*8 safety, pgrow, pshrnk, errcon
      dimension y(n), dy(n), yscal(n), yerr(n), ytemp(n)
      dimension temp(n,6)
      parameter (safety=0.9d0, pgrow=-.2d0, pshrnk=-.25d0,
     1           errcon=1.89d-04)
      h=htry
 1    call rkck(y,dy,ytemp,yerr,temp,x,h,n)
      errmax=0.d0
      do 10 i=1,n
        errmax=max(errmax,abs(yerr(i)/yscal(i)))
 10   continue
      errmax=errmax/eps
      if(errmax.gt.1.d0) then
         htemp=safety*h*(errmax**pshrnk)
         h=sign(max(abs(htemp),0.1d0*abs(h)),h)
         xnew=x+h
         if(xnew.eq.x) then
            call lnkerr('stepsize underflow in rkqs')
         endif
      else
         if(errmax.gt.errcon)then
            hnext=safety*h*(errmax**pgrow)
         else
            hnext=5.d0*h
        endif
        hdid=h
        x=x+h
        do 20 i=1,n
           y(i)=ytemp(i)
 20     continue
        return
      endif
      end


