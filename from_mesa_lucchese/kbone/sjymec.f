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
      implicit real*8(a-h,o-z)
      real*8 j,jx,jnmp,jnm,jnm1
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
c      if(abs(ynmp).gt.bound) 2,22
      if(abs(ynmp).gt.bound)go to 2
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
