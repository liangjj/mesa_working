      subroutine rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
      integer n,nmax
      real eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
      external derivs
      parameter (nmax=50)
cu    uses derivs,rkck
      integer i
      real errmax,h,htemp,xnew,yerr(nmax),ytemp(nmax),safety,pgrow,
     *pshrnk,errcon
      parameter (safety=0.9,pgrow=-.2,pshrnk=-.25,errcon=1.89e-4)
      h=htry
1     call rkck(y,dydx,n,x,h,ytemp,yerr,derivs)
      errmax=0.
      do 11 i=1,n
        errmax=max(errmax,abs(yerr(i)/yscal(i)))
11    continue
      errmax=errmax/eps
      if(errmax.gt.1.)then
        htemp=safety*h*(errmax**pshrnk)
        h=sign(max(abs(htemp),0.1*abs(h)),h)
        xnew=x+h
        if(xnew.eq.x)pause 'stepsize underflow in rkqs'
        goto 1
      else
        if(errmax.gt.errcon)then
          hnext=safety*h*(errmax**pgrow)
        else
          hnext=5.*h
        endif
        hdid=h
        x=x+h
        do 12 i=1,n
          y(i)=ytemp(i)
12      continue
        return
      endif
      end


