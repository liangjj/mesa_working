      program coul
c
c bessel fuction spline code with new prescription for
c the outgoing wave continuum function.
c
      real mu
      parameter (nbig=15000,lbig=11,nsmall=4015)
      parameter (nchmx=20,nemx=100)
      parameter (nsmall2=2*nsmall)
      dimension ipow(0:lbig),echan(nchmx),energy(nemx)
      common /integ/ xstart,ystart,ypstart,al,mu,h,xmax
      common /energy/ e,eta,ak
      common/spl/rd,r(nsmall),w(nsmall,0:lbig),cs(nsmall),rd26,l
      complex csc(nsmall,0:lbig),scc(nsmall)
      real aj(0:(lbig+1)),ay(0:(lbig+1))
      real aj1(0:(lbig+1)),ay1(0:(lbig+1))
      real csl(nsmall,0:lbig),derl(nsmall,0:lbig)
      complex x1,x2,f1,f2,f,fp,gp,sig(0:lbig),tv
      complex yc1,yc2,zl
      real css(nsmall)
      complex c1,c2,ai,ww(nsmall),wwl(nsmall,0:lbig)
      external zero,gg
      dimension der(nsmall),xx(nbig),yy(nbig),zz(nbig)
      dimension scr(nsmall)
      dimension wwlr(nsmall2,0:lbig),cscr(nsmall2,0:lbig)
      equivalence (wwl(1),wwlr(1)),(csc(1),cscr(1))
      ai=(0.,1.)
      ipr=1
c
c..unicos
c      call link ("unit5=infreeg,unit6=(outfree,hc,create)
c..unicos     1,unit8=bessplr//")
c
      open(5,file='infree')
      open(6,file='outfree')
      open(8,file='bessplr',form='unformatted')
c
      ipow(0)=0
      do 66 i=1,lbig
   66 ipow(i)=1
c     call keep80(4hwave)
c     call frid(10hxerox+film,1,4,1)
c     call dders(-1)
      eta=0.
      znuc=0.
      mu=1.
c
c lmax=max. l value desired
c nper=no. of points per unit a.u. for splined functions
c xstart = first r-point (.00001 is o.k.)
c xmax = last r-point(r=k*x!)
c nint = multiple of nper for integration(3 or 4 is reasonable)
c alpha = exponential parameter in test function
c
      read(5,*)lmax,nper,xstart,xmax,nint,alpha
      write(6,*)' Greens Function Code '
      write(6,101)lmax,nper,xstart,xmax,nint,alpha
 101  format(/,' L-Max  = ',i5,/,
     2         ' Nper   = ',i5,/,
     3         ' XStart = ',f20.12,/,
     4         ' XMax   = ',f20.12,/,
     5         ' Nint   = ',i5,/,
     6         ' Alpha  = ',f20.12)
      np=nper*xmax+1
      if(np.gt.nsmall)then
      write(6,100)
100   format(" np.gt.nsmall")
      stop
      endif
      rd=(xmax-xstart)/(np-1)
      rd26=rd*rd/6.
      do 1 i=1,np
1     r(i)=rd*(i-1.)+xstart
      aint=nint
      h=rd/aint
      nl=lmax+1
      do 668 i=1,np
      x=r(i)
      call sjymec(x,aj,ay,ncomp,lmax)
      do 668 j=0,lmax
668   w(i,j)=exp(-x*alpha)*aj(j)*x
667   continue
      do 33 l=0,lmax
      al=l
      ystart=xstart**(l+1)
      ypstart=(al+1.)*xstart**l
      yp1=(w(2,l)-w(1,l))/rd
      yp2=(w(np,l)-w(np-1,l))/rd
      call spline(r,w(1,l),np,yp1,yp2,scr,cs)
      call outward(xx,yy,index,zero)
      call outward(xx,zz,index,gg)
      i1=index
      i2=index-1./h
      x1=xx(i1)
      x2=xx(i2)
      det=yy(i1)*zz(i2)-yy(i2)*zz(i1)
      call sjymec(x1,aj,ay,ncomp,l)
      call sjymec(x2,aj1,ay1,ncomp,l)
      f1=x1*(-ay(l)+ai*aj(l))
      f2=x2*(-ay1(l)+ai*aj1(l))
      c1=f1*zz(i2)-f2*zz(i1)
      c1=c1/det
      c2=f2*yy(i1)-f1*yy(i2)
      c2=c2/det
      j=0
      do 3 i=1,index,nint
      j=j+1
      ww(j)=c1*yy(i)+c2*zz(i)
      ww(j)=ww(j)/xx(i)
      wwl(j,l)=ww(j)
3     continue
      yc1=(ww(2)-ww(1))/rd
      yc2=(ww(np)-ww(np-1))/rd
      call splinec(r,ww,np,yc1,yc2,scc,csc(1,l))
      c2r=c2
      do 5 i=1,np
      der(i)=w(i,l)*c2r
      der(i)=der(i)/r(i)**ipow(l)
      derl(i,l)=der(i)
5     continue
      yp1=(der(2)-der(1))/rd
      yp2=(der(np)-der(np-1))/rd
      call spline(r,der,np,yp1,yp2,scr,csl(1,l))
      if(ipr.ne.0)then
        write(6,*)' l=',l
        write(6,*)' ww'
        write(6,200)(ww(jj),jj=1,np)
200     format(5x,4e15.8)
        write(6,*)' der'
        write(6,200)(der(jj),jj=1,np)
        write(6,*)' cj'
        write(6,200)(csc(jj,l),jj=1,np)
        write(6,*)' cy'
        write(6,200)(csl(jj,l),jj=1,np)
      endif
33    continue
      write(8)lmax,np,xstart,rd,alpha
      write(8)(r(i),i=1,np)
      np2=2*np
      write(8)((wwlr(i,k),i=1,np2),k=0,lmax)
      write(8)((derl(i,k),i=1,np),k=0,lmax)
      write(8)((cscr(i,k),i=1,np2),k=0,lmax)
      write(8)((csl(i,k),i=1,np),k=0,lmax)
      write(6,*)' BessplG Finished '
 
      stop
      end
      subroutine splinec(x,y,n,yp1,ypn,scr,y2)
      complex y,yp1,ypn,scr,y2,p,qn,un,sig
      dimension x(n),y(n),scr(n),y2(n)
      y2(1)=-.5
      scr(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      do 11 i=2,n-1
      sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
      p=sig*y2(i-1)+2.
      y2(i)=(sig-1.)/p
11    scr(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     1 /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*scr(i-1))/p
      qn=.5
      un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      y2(n)=(un-qn*scr(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
      y2(k)=y2(k)*y2(k+1)+scr(k)
12    continue
      return
      end
      function zero(x)
      zero=0.
      return
      end
      subroutine outward(xx,yy,index,g)
      real k1,k2,k3
      real mu
      dimension xx(1),yy(1)
      common /integ/ xstart,ystart,ypstart,al,mu,h,xmax
      common /energy/ e,eta,ak
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
      subroutine spline(x,y,n,yp1,ypn,scr,y2)
      real y,yp1,ypn,scr,y2,p,qn,un,sig
      dimension x(n),y(n),scr(n),y2(n)
      y2(1)=-.5
      scr(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      do 11 i=2,n-1
      sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
      p=sig*y2(i-1)+2.
      y2(i)=(sig-1.)/p
11    scr(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     1 /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*scr(i-1))/p
      qn=.5
      un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      y2(n)=(un-qn*scr(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
      y2(k)=y2(k)*y2(k+1)+scr(k)
12    continue
      return
      end
      function gg(x)
      parameter (nsmall=4015,lbig=11)
      common /integ/ xstart,ystart,ypstart,al,mu,h,xmax
      common /energy/ e,eta,ak
      common/spl/rd,r(nsmall),w(nsmall,0:lbig),cs(nsmall),rd26,l
      klo=(x-xstart)/rd
      klo=klo+1
      a=(r(klo+1)-x)/rd
      b=(x-r(klo))/rd
      gg=a*w(klo,l)+b*w(klo+1,l)+(a*(a*a-1.)*cs(klo)+b*
     1 (b*b-1.)*cs(klo+1))*rd26
      return
      end
      subroutine sjymec(x,j,y,mn,nout)
c
c
c
c     jsubn(x),ysubn(x)      spherical bessel functions
c     do not use x=0.
c     jsub0(0)=1.   jsub n (0)=0.   for all n.gt.0
c     ysub n (0)= infinity  for all n.ge.0.
c     y sub n (x)=-j sub (-(n+1))  (x)        for all n.ge.0
c     mn  maximum number of subscript of bessel functions computed
c     nout   maximum subscript of desired bessel functions
c     recurrence techniques for bessel functions
c     fr.  mechel   math. of comp.  vol. 22   no. 101  jan. 1968
c
c
c
      real j,jx,jnmp,jnm,jnm1
      dimension j(0:1),y(0:1)
c
c     the dimension of j and y may be varied to suit the users needs
c     set dimension in calling routine  j((0,nout+5)),y((0,nout+5))
c
c
c
      data  bound /1.e13/, nadd /10/,  nxadd /10/
c
c     bound  nadd   nxadd   may be varied
c     nxadd controls the number of terms in the continued fraction
c
      mn=1
      j(0)=sin(x)/x
      y(0)=-cos(x)/x
      y(1)=y(0)/x-j(0)
      j(1)=j(0)/x+y(0)
      nmin=1
      nmax= nout
c
c     upward recursion valid for n<x for both bessel functions
c
      ix=x
      nx=min0(nout,ix)
      do 66 n=1,nx
      mn=n+1
      twon1=n+mn
      j(mn)=twon1*j(n)/x-j(n-1)
   66 y(mn)=twon1*y(n)/x-y(n-1)
      if(nout.le.mn) return
      jx=j(mn)
      mnout=mn
      nmin=mn
c
c     upward recursion for ysubn(x) to obtain good guess for jsubn(x)
c
      do 1 n=nmin,nout
      mn=n+1
      twon1=n+mn
    1 y(mn)=twon1*y(n)/x-y(n-1)
      ynm1=y(nout-1)
      ynm=y(nout)
      nmin=nout
      nmax=nout+nadd
   33 do 11 n=nmin,nmax
      mn=n+1
      twon1=n+mn
      ynmp=twon1*ynm/x-ynm1
      ynm1=ynm
      ynm=ynmp
   11 continue
      if(abs(ynmp).gt.bound) 2,22
   22 nmin=nmax+1
      nmax=nmax+nadd
      go to 33
    2 x2=x*x
      jnm=-1.0/(x2*ynmp)
      pmin1=1.
      qmin1=0.
      q0=1.0
c
c     develop continued fraction to obtain good guess for jsubn+1(x)
c
      b0=mn+mn+1
      p0=mn+mn+1
      do 4 n=1,nxadd
      b0=b0+2.
      p1=b0*p0-x2*pmin1
      q1=b0*q0-x2*qmin1
      pmin1=p0
      p0=p1
      qmin1=q0
      q0=q1
    4 continue
c
      gnmax1=p1/q1
      jnmp=x*jnm/gnmax1
      do 55 n=mn-1,nout+1,-1
      nm=n+1
      twon1=n+nm
      jnm1=twon1*jnm/x-jnmp
      jnmp=jnm
      jnm=jnm1
   55 continue
      j(nout+1)=jnmp
      j(nout)=jnm1
      do 5 n=nout,mnout+1,-1
      nm=n+1
      twon1=n+nm
    5 j(n-1)=twon1*j(n)/x-j(n+1)
      c=jx/j(mnout)
c
c     correct bessel functions with weighting factor
c
      if (x.lt.1.e-3) return
      do 6 n=mnout,nout
    6 j(n)=c*j(n)
c
      return
      end
