      subroutine rkck(y,dydx,n,x,h,yout,yerr,derivs)
      integer n,nmax
      real h,x,dydx(n),y(n),yerr(n),yout(n)
      external derivs
      parameter (nmax=50)
cu    uses derivs
      integer i
      real ak2(nmax),ak3(nmax),ak4(nmax),ak5(nmax),ak6(nmax),
     *ytemp(nmax),a2,a3,a4,a5,a6,b21,b31,b32,b41,b42,b43,b51,b52,b53,
     *b54,b61,b62,b63,b64,b65,c1,c3,c4,c6,dc1,dc3,dc4,dc5,dc6
      parameter (a2=.2,a3=.3,a4=.6,a5=1.,a6=.875,b21=.2,b31=3./40.,
     *b32=9./40.,b41=.3,b42=-.9,b43=1.2,b51=-11./54.,b52=2.5,
     *b53=-70./27.,b54=35./27.,b61=1631./55296.,b62=175./512.,
     *b63=575./13824.,b64=44275./110592.,b65=253./4096.,c1=37./378.,
     *c3=250./621.,c4=125./594.,c6=512./1771.,dc1=c1-2825./27648.,
     *dc3=c3-18575./48384.,dc4=c4-13525./55296.,dc5=-277./14336.,
     *dc6=c6-.25)
      do 11 i=1,n
        ytemp(i)=y(i)+b21*h*dydx(i)
11    continue
      call derivs(x+a2*h,ytemp,ak2)
      do 12 i=1,n
        ytemp(i)=y(i)+h*(b31*dydx(i)+b32*ak2(i))
12    continue
      call derivs(x+a3*h,ytemp,ak3)
      do 13 i=1,n
        ytemp(i)=y(i)+h*(b41*dydx(i)+b42*ak2(i)+b43*ak3(i))
13    continue
      call derivs(x+a4*h,ytemp,ak4)
      do 14 i=1,n
        ytemp(i)=y(i)+h*(b51*dydx(i)+b52*ak2(i)+b53*ak3(i)+b54*ak4(i))
14    continue
      call derivs(x+a5*h,ytemp,ak5)
      do 15 i=1,n
        ytemp(i)=y(i)+h*(b61*dydx(i)+b62*ak2(i)+b63*ak3(i)+b64*ak4(i)+
     *b65*ak5(i))
15    continue
      call derivs(x+a6*h,ytemp,ak6)
      do 16 i=1,n
        yout(i)=y(i)+h*(c1*dydx(i)+c3*ak3(i)+c4*ak4(i)+c6*ak6(i))
16    continue
      do 17 i=1,n
        yerr(i)=h*(dc1*dydx(i)+dc3*ak3(i)+dc4*ak4(i)+dc5*ak5(i)+dc6*
     *ak6(i))
17    continue
      return
      end


