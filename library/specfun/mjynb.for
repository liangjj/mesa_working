	program mjynb
c
c       ====================================================
c       purpose: this program computes bessel functions 
c                jn(x) and yn(x), and their derivatives 
c                using subroutine jynb
c       input :  x --- argument of jn(x) & yn(x)  ( x ע 0 )
c                n --- order of jn(x) & yn(x)
c                      ( n = 0,1,2,תתת, n ף 250 )
c       output:  bj(n) --- jn(x)
c                dj(n) --- jn'(x)
c                by(n) --- yn(x)
c                dy(n) --- yn'(x)
c       example:
c                x = 10.0
c
c                n        jn(x)           jn'(x)
c              -------------------------------------
c                0    -.2459358d+00   -.4347275d-01
c               10     .2074861d+00    .8436958d-01
c               20     .1151337d-04    .2011954d-04
c               30     .1551096d-11    .4396479d-11
c
c                n        yn(x)           yn'(x)
c              -------------------------------------
c                0     .5567117d-01   -.2490154d+00
c               10    -.3598142d+00    .1605149d+00
c               20    -.1597484d+04    .2737803d+04
c               30    -.7256142d+10    .2047617d+11
c       ====================================================
c
	implicit double precision (a-h,o-z)
	dimension bj(0:250),by(0:250),dj(0:250),dy(0:250)
	write(*,*)'  please enter n, x'
	read(*,*)n,x
	write(*,30)n,x
	if (n.le.8) then
	   ns=1
	else
	   write(*,*)'  please enter order step ns'
	   read(*,*)ns
	endif
	write(*,*)
	call jynb(n,x,nm,bj,dj,by,dy)
	write(*,*)'  n        jn(x)           jn''(x)'
	write(*,*)'--------------------------------------'
	do 10 k=0,nm,ns
10         write(*,40)k,bj(k),dj(k)
	write(*,*)
	write(*,*)'  n        yn(x)           yn''(x)'
	write(*,*)'--------------------------------------'
	do 20 k=0,nm,ns
20         write(*,40)k,by(k),dy(k)
30      format(3x,6hnmax =,i3,',    ',3hx =,f6.1)
40      format(1x,i3,1x,2d16.7)
	end


	subroutine jynb(n,x,nm,bj,dj,by,dy)
c
c       =====================================================
c       purpose: compute bessel functions jn(x), yn(x) and
c                their derivatives
c       input :  x --- argument of jn(x) and yn(x) ( x ע 0 )
c                n --- order of jn(x) and yn(x)
c       output:  bj(n) --- jn(x)
c                dj(n) --- jn'(x)
c                by(n) --- yn(x)
c                dy(n) --- yn'(x)
c                nm --- highest order computed
c       routines called:
c                msta1 and msta2 to calculate the starting 
c                point for backward recurrence
c       =====================================================
c
	implicit double precision (a-h,o-z)
	dimension bj(0:n),dj(0:n),by(0:n),dy(0:n),
     &            a(4),b(4),a1(4),b1(4)
	pi=3.141592653589793d0
	r2p=.63661977236758d0
	nm=n
	if (x.lt.1.0d-100) then
	   do 10 k=0,n
	      bj(k)=0.0d0
	      dj(k)=0.0d0
	      by(k)=-1.0d+300
10            dy(k)=1.0d+300
	   bj(0)=1.0d0
	   dj(1)=0.5d0
	   return
	endif
	if (x.le.300.0.or.n.gt.int(0.9*x)) then
	   if (n.eq.0) nm=1
	   m=msta1(x,200)
	   if (m.lt.nm) then
	      nm=m
	   else
	      m=msta2(x,nm,15)
	   endif
	   bs=0.0d0
	   su=0.0d0
	   sv=0.0d0
	   f2=0.0d0
	   f1=1.0d-100
	   do 15 k=m,0,-1
	      f=2.0d0*(k+1.0d0)/x*f1-f2
	      if (k.le.nm) bj(k)=f
	      if (k.eq.2*int(k/2).and.k.ne.0) then
		 bs=bs+2.0d0*f
		 su=su+(-1)**(k/2)*f/k
	      else if (k.gt.1) then
		 sv=sv+(-1)**(k/2)*k/(k*k-1.0)*f
	      endif
	      f2=f1
15            f1=f
	   s0=bs+f
	   do 20 k=0,nm
20            bj(k)=bj(k)/s0
	   ec=dlog(x/2.0d0)+0.5772156649015329d0
	   by0=r2p*(ec*bj(0)-4.0d0*su/s0)
	   by(0)=by0
	   by1=r2p*((ec-1.0d0)*bj(1)-bj(0)/x-4.0d0*sv/s0)
	   by(1)=by1
	else
	   data a/-.7031250000000000d-01,.1121520996093750d+00,
     &            -.5725014209747314d+00,.6074042001273483d+01/
	   data b/ .7324218750000000d-01,-.2271080017089844d+00,
     &             .1727727502584457d+01,-.2438052969955606d+02/
	   data a1/.1171875000000000d+00,-.1441955566406250d+00,
     &             .6765925884246826d+00,-.6883914268109947d+01/
	   data b1/-.1025390625000000d+00,.2775764465332031d+00,
     &             -.1993531733751297d+01,.2724882731126854d+02/
	   t1=x-0.25d0*pi
	   p0=1.0d0
	   q0=-0.125d0/x
	   do 25 k=1,4
	      p0=p0+a(k)*x**(-2*k)
25            q0=q0+b(k)*x**(-2*k-1)
	   cu=dsqrt(r2p/x)
	   bj0=cu*(p0*dcos(t1)-q0*dsin(t1))
	   by0=cu*(p0*dsin(t1)+q0*dcos(t1))
	   bj(0)=bj0
	   by(0)=by0
	   t2=x-0.75d0*pi
	   p1=1.0d0
	   q1=0.375d0/x
	   do 30 k=1,4
	      p1=p1+a1(k)*x**(-2*k)
30            q1=q1+b1(k)*x**(-2*k-1)
	   bj1=cu*(p1*dcos(t2)-q1*dsin(t2))
	   by1=cu*(p1*dsin(t2)+q1*dcos(t2))
	   bj(1)=bj1
	   by(1)=by1
	   do 35 k=2,nm
	      bjk=2.0d0*(k-1.0d0)/x*bj1-bj0
	      bj(k)=bjk
	      bj0=bj1
35            bj1=bjk
	endif
	dj(0)=-bj(1)
	do 40 k=1,nm
40         dj(k)=bj(k-1)-k/x*bj(k)
	do 45 k=2,nm
	   byk=2.0d0*(k-1.0d0)*by1/x-by0
	   by(k)=byk
	   by0=by1
45         by1=byk
	dy(0)=-by(1)
	do 50 k=1,nm
50         dy(k)=by(k-1)-k*by(k)/x
	return
	end


	integer function msta1(x,mp)
c
c       ===================================================
c       purpose: determine the starting point for backward  
c                recurrence such that the magnitude of    
c                jn(x) at that point is about 10^(-mp)
c       input :  x     --- argument of jn(x)
c                mp    --- value of magnitude
c       output:  msta1 --- starting point   
c       ===================================================
c
	implicit double precision (a-h,o-z)
	a0=dabs(x)
	n0=int(1.1*a0)+1
	f0=envj(n0,a0)-mp
	n1=n0+5
	f1=envj(n1,a0)-mp
	do 10 it=1,20             
	   nn=n1-(n1-n0)/(1.0d0-f0/f1)                  
	   f=envj(nn,a0)-mp
	   if(abs(nn-n1).lt.1) go to 20
	   n0=n1
	   f0=f1
	   n1=nn
 10        f1=f
 20     msta1=nn
	return
	end


	integer function msta2(x,n,mp)
c
c       ===================================================
c       purpose: determine the starting point for backward
c                recurrence such that all jn(x) has mp
c                significant digits
c       input :  x  --- argument of jn(x)
c                n  --- order of jn(x)
c                mp --- significant digit
c       output:  msta2 --- starting point
c       ===================================================
c
	implicit double precision (a-h,o-z)
	a0=dabs(x)
	hmp=0.5d0*mp
	ejn=envj(n,a0)
	if (ejn.le.hmp) then
	   obj=mp
	   n0=int(1.1*a0)
	else
	   obj=hmp+ejn
	   n0=n
	endif
	f0=envj(n0,a0)-obj
	n1=n0+5
	f1=envj(n1,a0)-obj
	do 10 it=1,20
	   nn=n1-(n1-n0)/(1.0d0-f0/f1)
	   f=envj(nn,a0)-obj
	   if (abs(nn-n1).lt.1) go to 20
	   n0=n1
	   f0=f1
	   n1=nn
10         f1=f
20      msta2=nn+10
	return
	end

	real*8 function envj(n,x)
	double precision x
	envj=0.5d0*dlog10(6.28d0*n)-n*dlog10(1.36d0*x/n)
	return
	end
