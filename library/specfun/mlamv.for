	program mlamv
c
c       =======================================================
c       purpose: this program computes the lambda functions 
c                for an arbitrary order, and their derivative 
c                using subroutine lamv
c       input :  x --- argument of lambda function
c                v --- order of lambda function
c                      ( v = n+v0, 0 ó n ó 250, 0 ó v0 < 1 )
c       output:  vl(n) --- lambda function of order n+v0
c                dl(n) --- derivative of lambda function 
c       example: x = 10.0
c
c                   v         lambda(x)        lambda'(x)
c                 ------------------------------------------
c                  0.25    -.12510515d+00    -.78558916d-01
c                  0.50    -.54402111d-01    -.78466942d-01
c                  0.75    -.13657787d-01    -.66234027d-01
c                  1.00     .86945492d-02    -.50926063d-01
c                  1.25     .19639729d-01    -.36186221d-01
c                  1.50     .23540083d-01    -.23382658d-01
c                  1.75     .23181910d-01    -.12893894d-01
c                  2.00     .20370425d-01    -.46703503d-02
c                  2.25     .16283799d-01     .15101684d-02
c                  2.50     .11691329d-01     .59243767d-02
c       =======================================================
c
	implicit double precision (a-h,o-z)
	dimension vl(0:250),dl(0:250)
	write(*,*)'  please enter v and x '
	read(*,*)v,x
	write(*,20)v,x
	if (v.le.8) then
	   ns=1
	else
	   write(*,*)'  please enter order step ns'
	   read(*,*)ns
	endif
	write(*,*)
	write(*,*) '   v         lambda(x)        lambda''(x)'
	write(*,*)'-------------------------------------------'
	call lamv(v,x,vm,vl,dl)
	nm=int(vm)
	v0=vm-nm
	do 10 k=0,nm,ns
	   vk=k+v0
10         write(*,15)vk,vl(k),dl(k)
15      format(1x,f6.2,2d18.8)
20      format(1x,'v =',f6.2,'    ','x =',f8.2)
	end


	subroutine lamv(v,x,vm,vl,dl)
c
c       =========================================================
c       purpose: compute lambda function with arbitrary order v,
c                and their derivative
c       input :  x --- argument of lambda function
c                v --- order of lambda function 
c       output:  vl(n) --- lambda function of order n+v0
c                dl(n) --- derivative of lambda function 
c                vm --- highest order computed
c       routines called:
c            (1) msta1 and msta2 for computing the starting 
c                point for backward recurrence
c            (2) gam0 for computing gamma function (|x| ó 1)
c       =========================================================
c
	implicit double precision (a-h,o-z)
	dimension vl(0:*),dl(0:*)
	pi=3.141592653589793d0
	rp2=0.63661977236758d0
	x=dabs(x)
	x2=x*x
	n=int(v)
	v0=v-n
	vm=v
	if (x.le.12.0d0) then
	   do 25 k=0,n
	      vk=v0+k
	      bk=1.0d0
	      r=1.0d0
	      do 10 i=1,50
		 r=-0.25d0*r*x2/(i*(i+vk))
		 bk=bk+r
		 if (dabs(r).lt.dabs(bk)*1.0d-15) go to 15
10            continue
15            vl(k)=bk
	      uk=1.0d0
	      r=1.0d0
	      do 20 i=1,50
		 r=-0.25d0*r*x2/(i*(i+vk+1.0d0))
		 uk=uk+r
		 if (dabs(r).lt.dabs(uk)*1.0d-15) go to 25
20            continue
25            dl(k)=-0.5d0*x/(vk+1.0d0)*uk
	   return
	endif
	k0=11
	if (x.ge.35.0d0) k0=10
	if (x.ge.50.0d0) k0=8
	do 40 j=0,1
	   vv=4.0d0*(j+v0)*(j+v0)
	   px=1.0d0
	   rp=1.0d0
	   do 30 k=1,k0
	      rp=-0.78125d-2*rp*(vv-(4.0*k-3.0)**2.0)*(vv-
     &            (4.0*k-1.0)**2.0)/(k*(2.0*k-1.0)*x2)
30            px=px+rp
	   qx=1.0d0
	   rq=1.0d0
	   do 35 k=1,k0
	      rq=-0.78125d-2*rq*(vv-(4.0*k-1.0)**2.0)*(vv-
     &            (4.0*k+1.0)**2.0)/(k*(2.0*k+1.0)*x2)
35            qx=qx+rq
	   qx=0.125d0*(vv-1.0d0)*qx/x
	   xk=x-(0.5d0*(j+v0)+0.25d0)*pi
	   a0=dsqrt(rp2/x)
	   ck=dcos(xk)
	   sk=dsin(xk)
	   if (j.eq.0) bjv0=a0*(px*ck-qx*sk)
	   if (j.eq.1) bjv1=a0*(px*ck-qx*sk)
40      continue
	if (v0.eq.0.0d0) then
	   ga=1.0d0
	else
	   call gam0(v0,ga)
	   ga=v0*ga
	endif
	fac=(2.0d0/x)**v0*ga
	vl(0)=bjv0
	dl(0)=-bjv1+v0/x*bjv0
	vl(1)=bjv1
	dl(1)=bjv0-(1.0d0+v0)/x*bjv1
	r0=2.0d0*(1.0d0+v0)/x
	if (n.le.1) then
	   vl(0)=fac*vl(0)
	   dl(0)=fac*dl(0)-v0/x*vl(0)
	   vl(1)=fac*r0*vl(1)
	   dl(1)=fac*r0*dl(1)-(1.0d0+v0)/x*vl(1)
	   return
	endif
	if (n.ge.2.and.n.le.int(0.9*x)) then
	   f0=bjv0
	   f1=bjv1
	   do 45 k=2,n
	      f=2.0d0*(k+v0-1.0d0)/x*f1-f0
	      f0=f1
	      f1=f
45            vl(k)=f
	else if (n.ge.2) then
	   m=msta1(x,200)
	   if (m.lt.n) then
	      n=m
	   else
	      m=msta2(x,n,15)
	   endif
	   f2=0.0d0
	   f1=1.0d-100
	   do 50 k=m,0,-1
	      f=2.0d0*(v0+k+1.0d0)/x*f1-f2
	      if (k.le.n) vl(k)=f
	      f2=f1
50            f1=f
	   if (dabs(bjv0).gt.dabs(bjv1)) cs=bjv0/f
	   if (dabs(bjv0).le.dabs(bjv1)) cs=bjv1/f2
	   do 55 k=0,n
55            vl(k)=cs*vl(k)
	endif
	vl(0)=fac*vl(0)
	do 65 j=1,n
	   rc=fac*r0
	   vl(j)=rc*vl(j)
	   dl(j-1)=-0.5d0*x/(j+v0)*vl(j)
65         r0=2.0d0*(j+v0+1)/x*r0
	dl(n)=2.0d0*(v0+n)*(vl(n-1)-vl(n))/x
	vm=n+v0
	return
	end


	subroutine gam0 (x,ga)
c
c       ================================================
c       purpose: compute gamma function â(x)
c       input :  x  --- argument of â(x)  ( |x| ó 1 )
c       output:  ga --- â(x)
c       ================================================
c
	implicit double precision (a-h,o-z)
	dimension g(25)
	data g/1.0d0,0.5772156649015329d0,
     &       -0.6558780715202538d0, -0.420026350340952d-1,
     &        0.1665386113822915d0, -.421977345555443d-1,
     &        -.96219715278770d-2, .72189432466630d-2,
     &        -.11651675918591d-2, -.2152416741149d-3,
     &         .1280502823882d-3, -.201348547807d-4,
     &        -.12504934821d-5, .11330272320d-5,
     &        -.2056338417d-6, .61160950d-8,
     &         .50020075d-8, -.11812746d-8,
     &         .1043427d-9, .77823d-11,
     &        -.36968d-11, .51d-12,
     &        -.206d-13, -.54d-14, .14d-14/
	gr=(25)
	do 20 k=24,1,-1
20         gr=gr*x+g(k)
	ga=1.0d0/(gr*x)
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
