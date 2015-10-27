      subroutine outward(xx,yy,index,g)
      implicit real*8 (a-h,o-z)
      real*8 k1,k2,k3
      real*8 mu
      dimension xx(1),yy(1)
      common /integ/ xstart,ystart,ypstart,al,mu,h,xmax
      common /energy/ e,eta,ak
c
      f(z)=((2.*eta/z-1.)+al*(al+1.)/z/z)*t-g(z)*2.
      tmu=2.*mu
c initialize solution
      itest=0
      h8=h/8.
      h2=h/2.
      yp=ypstart
      y=ystart
      x=xstart
      index=0
c begin integration
    1 t=y
      k1=h*f(x)
      index=index+1
      xx(index)=x
      yy(index)=y
      itest=itest+1
      t=y+h2*yp+h8*k1
      k2=h*f(x+h2)
      t=y+h*yp+h2*k2
      k3=h*f(x+h)
      ypp=yp+k1/6.+2.*k2/3.+k3/6.
      yn=y+h*(yp+(k1+2.*k2)/6.)
      xn=x+h
      if(xn.gt.xmax)go to 23
      yp=ypp
      y=yn
      x=xn
      go to 1
   23 continue
      return
      end
