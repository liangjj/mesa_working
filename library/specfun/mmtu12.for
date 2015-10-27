	program mmtu12
c
c       ===============================================================
c       purpose: this program computes the modified mathieu functions 
c                of the first and second kinds, mcm(1)(2)(x,q) and 
c                msm(1)(2)(x,q), and their derivatives using 
c                subroutine mtu12
c       input:   kf --- function code
c                       kf=1 for computing mcm(x,q)
c                       kf=2 for computing msm(x,q)
c                kc --- function code
c                       kc=1 for computing mcm(1)(x,q) and mcm(1)'(x,q)
c                            or msm(1)(x,q) and msm(1)'(x,q)
c                       kc=2 for computing mcm(2)(x,q) and mcm(2)'(x,q)
c                            or msm(2)(x,q) and msm(2)'(x,q)
c                       kc=3 for both modified mathieu functions of the
c                            first and second kinds, and their
c                            derivatives
c                m  --- order of mathieu functions
c                q  --- parameter of mathieu functions
c                x  --- argument of mathieu functions
c       output:  f1r --- mcm(1)(x,q) or msm(1)(x,q)
c                d1r --- derivative of mcm(1)(x,q) or msm(1)(x,q)
c                f2r --- mcm(2)(x,q) or msm(2)(x,q)
c                d2r --- derivative of mcm(2)(x,q) or msm(2)(x,q)
c       ===============================================================
c
	implicit double precision (a-h,o-z)
	write(*,*)'please enter kf, m, q and x '
	read(*,*)kf,m,q,x
	write(*,10)kf,m,q,x
	kc=3
	call mtu12(kf,kc,m,q,x,f1r,d1r,f2r,d2r)
	write(*,*)
	if (kf.eq.1) then
	   write(*,*)'   x      mcm(1)(x,q)    mcm(1)''(x,q)',
     &               '    mcm(2)(x,q)     mcm(2)''(x,q)'
	else
	   write(*,*)'   x      msm(1)(x,q)    msm(1)''(x,q)',
     &               '    msm(2)(x,q)     msm(2)''(x,q)'
	endif
	write(*,*)' --------------------------------------',
     &            '-------------------------------'
	write(*,20)x,f1r,d1r,f2r,d2r
	write(*,*)
	write(*,30)f1r*d2r-f2r*d1r,.63661977236758d0
	write(*,40)
10      format(1x,4hkf =,i2,',  ',3hm =,i3,',  ',
     &         3hq =,f5.1,',  ',3hx =,f5.1)
20      format(1x,f5.1,4d16.8)
30      format(1x,'wronskian=',e16.8,3x,'should equal   2/pi=',e16.8)
40      format(1x,/1x,'caution: this check is not accurate if it ',
     &        'involves',/1x,'         the subtraction of two ',
     &        'similar numbers')
	end


	subroutine mtu12(kf,kc,m,q,x,f1r,d1r,f2r,d2r)
c
c       ==============================================================
c       purpose: compute modified mathieu functions of the first and
c                second kinds, mcm(1)(2)(x,q) and msm(1)(2)(x,q),
c                and their derivatives
c       input:   kf --- function code
c                       kf=1 for computing mcm(x,q)
c                       kf=2 for computing msm(x,q)
c                kc --- function code
c                       kc=1 for computing the first kind
c                       kc=2 for computing the second kind
c                            or msm(2)(x,q) and msm(2)'(x,q)
c                       kc=3 for computing both the first
c                            and second kinds
c                m  --- order of mathieu functions
c                q  --- parameter of mathieu functions ( q ò 0 )
c                x  --- argument of mathieu functions
c       output:  f1r --- mcm(1)(x,q) or msm(1)(x,q)
c                d1r --- derivative of mcm(1)(x,q) or msm(1)(x,q)
c                f2r --- mcm(2)(x,q) or msm(2)(x,q)
c                d2r --- derivative of mcm(2)(x,q) or msm(2)(x,q)
c       routines called:
c            (1) cva2 for computing the characteristic values
c            (2) fcoef for computing expansion coefficients
c            (3) jynb for computing jn(x), yn(x) and their
c                derivatives
c       ==============================================================
c
	implicit double precision (a-h,o-z)
	dimension fg(251),bj1(0:251),dj1(0:251),bj2(0:251),dj2(0:251),
     &            by1(0:251),dy1(0:251),by2(0:251),dy2(0:251)
	eps=1.0d-14
	if (kf.eq.1.and.m.eq.2*int(m/2)) kd=1
	if (kf.eq.1.and.m.ne.2*int(m/2)) kd=2
	if (kf.eq.2.and.m.ne.2*int(m/2)) kd=3
	if (kf.eq.2.and.m.eq.2*int(m/2)) kd=4
	call cva2(kd,m,q,a)
	if (q.le.1.0d0) then
	   qm=7.5+56.1*sqrt(q)-134.7*q+90.7*sqrt(q)*q
	else
	   qm=17.0+3.1*sqrt(q)-.126*q+.0037*sqrt(q)*q
	endif
	km=int(qm+0.5*m)              
	call fcoef(kd,m,q,a,fg)
	ic=int(m/2)+1
	if (kd.eq.4) ic=m/2
	c1=dexp(-x)
	c2=dexp(x)
	u1=dsqrt(q)*c1
	u2=dsqrt(q)*c2
	call jynb(km,u1,nm,bj1,dj1,by1,dy1)
	call jynb(km,u2,nm,bj2,dj2,by2,dy2)
	if (kc.eq.2) go to 50
	f1r=0.0d0
	do 30 k=1,km
	   if (kd.eq.1) then
	      f1r=f1r+(-1)**(ic+k)*fg(k)*bj1(k-1)*bj2(k-1)
	   else if (kd.eq.2.or.kd.eq.3) then
	      f1r=f1r+(-1)**(ic+k)*fg(k)*(bj1(k-1)*bj2(k)
     &            +(-1)**kd*bj1(k)*bj2(k-1))
	   else
	      f1r=f1r+(-1)**(ic+k)*fg(k)*(bj1(k-1)*bj2(k+1)
     &            -bj1(k+1)*bj2(k-1))
	   endif
	   if (k.ge.5.and.dabs(f1r-w1).lt.dabs(f1r)*eps) go to 35
30         w1=f1r
35      f1r=f1r/fg(1)
	d1r=0.0d0
	do 40 k=1,km
	   if (kd.eq.1) then
	      d1r=d1r+(-1)**(ic+k)*fg(k)*(c2*bj1(k-1)*dj2(k-1)
     &            -c1*dj1(k-1)*bj2(k-1))
	   else if (kd.eq.2.or.kd.eq.3) then
	      d1r=d1r+(-1)**(ic+k)*fg(k)*(c2*(bj1(k-1)*dj2(k)
     &            +(-1)**kd*bj1(k)*dj2(k-1))-c1*(dj1(k-1)*bj2(k)
     &            +(-1)**kd*dj1(k)*bj2(k-1)))
	   else
	      d1r=d1r+(-1)**(ic+k)*fg(k)*(c2*(bj1(k-1)*dj2(k+1)
     &            -bj1(k+1)*dj2(k-1))-c1*(dj1(k-1)*bj2(k+1)
     &            -dj1(k+1)*bj2(k-1)))
	   endif
	   if (k.ge.5.and.dabs(d1r-w2).lt.dabs(d1r)*eps) go to 45
40         w2=d1r
45      d1r=d1r*dsqrt(q)/fg(1)
	if (kc.eq.1) return
50      f2r=0.0d0
	do 55 k=1,km
	   if (kd.eq.1) then
	      f2r=f2r+(-1)**(ic+k)*fg(k)*bj1(k-1)*by2(k-1)
	   else if (kd.eq.2.or.kd.eq.3) then
	      f2r=f2r+(-1)**(ic+k)*fg(k)*(bj1(k-1)*by2(k)
     &            +(-1)**kd*bj1(k)*by2(k-1))
	   else
	      f2r=f2r+(-1)**(ic+k)*fg(k)*(bj1(k-1)*by2(k+1)
     &            -bj1(k+1)*by2(k-1))
	   endif
	   if (k.ge.5.and.dabs(f2r-w1).lt.dabs(f2r)*eps) go to 60
55         w1=f2r
60      f2r=f2r/fg(1)
	d2r=0.0d0
	do 65 k=1,km
	   if (kd.eq.1) then
	      d2r=d2r+(-1)**(ic+k)*fg(k)*(c2*bj1(k-1)*dy2(k-1)
     &            -c1*dj1(k-1)*by2(k-1))
	   else if (kd.eq.2.or.kd.eq.3) then
	      d2r=d2r+(-1)**(ic+k)*fg(k)*(c2*(bj1(k-1)*dy2(k)
     &            +(-1)**kd*bj1(k)*dy2(k-1))-c1*(dj1(k-1)*by2(k)
     &            +(-1)**kd*dj1(k)*by2(k-1)))
	   else
	      d2r=d2r+(-1)**(ic+k)*fg(k)*(c2*(bj1(k-1)*dy2(k+1)
     &            -bj1(k+1)*dy2(k-1))-c1*(dj1(k-1)*by2(k+1)
     &            -dj1(k+1)*by2(k-1)))
	   endif
	   if (k.ge.5.and.dabs(d2r-w2).lt.dabs(d2r)*eps) go to 70
65         w2=d2r
70         d2r=d2r*dsqrt(q)/fg(1)
	return
	end


	subroutine fcoef(kd,m,q,a,fc)
c
c       =====================================================
c       purpose: compute expansion coefficients for mathieu
c                functions and modified mathieu functions
c       input :  m  --- order of mathieu functions
c                q  --- parameter of mathieu functions
c                kd --- case code
c                       kd=1 for cem(x,q)  ( m = 0,2,4,...)
c                       kd=2 for cem(x,q)  ( m = 1,3,5,...)
c                       kd=3 for sem(x,q)  ( m = 1,3,5,...)
c                       kd=4 for sem(x,q)  ( m = 2,4,6,...)
c                a  --- characteristic value of mathieu
c                       functions for given m and q
c       output:  fc(k) --- expansion coefficients of mathieu
c                       functions ( k= 1,2,...,km )
c                       fc(1),fc(2),fc(3),... correspond to
c                       a0,a2,a4,... for kd=1 case, a1,a3,
c                       a5,... for kd=2 case, b1,b3,b5,...
c                       for kd=3 case and b2,b4,b6,... for
c                       kd=4 case
c       =====================================================
c
	implicit double precision (a-h,o-z)
	dimension fc(251)
	if (q.le.1.0d0) then
	   qm=7.5+56.1*sqrt(q)-134.7*q+90.7*sqrt(q)*q
	else
	   qm=17.0+3.1*sqrt(q)-.126*q+.0037*sqrt(q)*q
	endif
	km=int(qm+0.5*m)                   
	if (q.eq.0.0d0) then
	   do 10 k=1,km
10            fc(k)=0.0d0
	   if (kd.eq.1) then
	      fc((m+2)/2)=1.0d0
	      if (m.eq.0) fc(1)=1.0d0/dsqrt(2.0d0)
	   else if (kd.eq.4) then
	      fc(m/2)=1.0d0
	   else
	      fc((m+1)/2)=1.0d0
	   endif
	   return
	endif
	kb=0
	s=0.0d0
	f=1.0d-100
	u=0.0d0
	fc(km)=0.0d0
	if (kd.eq.1) then
	   do 25 k=km,3,-1
	      v=u
	      u=f
	      f=(a-4.0d0*k*k)*u/q-v
	      if (dabs(f).lt.dabs(fc(k+1))) then
		 kb=k
		 fc(1)=1.0d-100
		 sp=0.0d0
		 f3=fc(k+1)
		 fc(2)=a/q*fc(1)
		 fc(3)=(a-4.0d0)*fc(2)/q-2.0d0*fc(1)
		 u=fc(2)
		 f1=fc(3)
		 do 15 i=3,kb
		    v=u
		    u=f1
		    f1=(a-4.0d0*(i-1.0d0)**2)*u/q-v
		    fc(i+1)=f1
		    if (i.eq.kb) f2=f1
		    if (i.ne.kb) sp=sp+f1*f1
15               continue
		 sp=sp+2.0d0*fc(1)**2+fc(2)**2+fc(3)**2
		 ss=s+sp*(f3/f2)**2
		 s0=dsqrt(1.0d0/ss)
		 do 20 j=1,km
		    if (j.le.kb+1) then
		       fc(j)=s0*fc(j)*f3/f2
		    else
		       fc(j)=s0*fc(j)
		    endif
20               continue
		 go to 85
	      else
		 fc(k)=f
		 s=s+f*f
	      endif
25         continue
	   fc(2)=q*fc(3)/(a-4.0d0-2.0d0*q*q/a)
	   fc(1)=q/a*fc(2)
	   s=s+2.0d0*fc(1)**2+fc(2)**2
	   s0=dsqrt(1.0d0/s)
	   do 30 k=1,km
30            fc(k)=s0*fc(k)
	else if (kd.eq.2.or.kd.eq.3) then
	   do 35 k=km,3,-1
	      v=u
	      u=f
	      f=(a-(2.0d0*k-1)**2)*u/q-v
	      if (dabs(f).ge.dabs(fc(k))) then
		 fc(k-1)=f
		 s=s+f*f
	      else
		 kb=k
		 f3=fc(k)
		 go to 45
	      endif
35         continue
	   fc(1)=q/(a-1.0d0-(-1)**kd*q)*fc(2)
	   s=s+fc(1)*fc(1)
	   s0=dsqrt(1.0d0/s)
	   do 40 k=1,km
40            fc(k)=s0*fc(k)
	   go to 85
45         fc(1)=1.0d-100
	   fc(2)=(a-1.0d0-(-1)**kd*q)/q*fc(1)
	   sp=0.0d0
	   u=fc(1)
	   f1=fc(2)
	   do 50 i=2,kb-1
	      v=u
	      u=f1
	      f1=(a-(2.0d0*i-1.0d0)**2)*u/q-v
	      if (i.ne.kb-1) then
		 fc(i+1)=f1
		 sp=sp+f1*f1
	      else
		 f2=f1
	      endif
50         continue
	   sp=sp+fc(1)**2+fc(2)**2
	   ss=s+sp*(f3/f2)**2
	   s0=1.0d0/dsqrt(ss)
	   do 55 j=1,km
	      if (j.lt.kb) fc(j)=s0*fc(j)*f3/f2
	      if (j.ge.kb) fc(j)=s0*fc(j)
55         continue
	else if (kd.eq.4) then
	   do 60 k=km,3,-1
	      v=u
	      u=f
	      f=(a-4.0d0*k*k)*u/q-v
	      if (dabs(f).ge.dabs(fc(k))) then
		 fc(k-1)=f
		 s=s+f*f
	      else
		 kb=k
		 f3=fc(k)
		 go to 70
	      endif
60         continue
	   fc(1)=q/(a-4.0d0)*fc(2)
	   s=s+fc(1)*fc(1)
	   s0=dsqrt(1.0d0/s)
	   do 65 k=1,km
65            fc(k)=s0*fc(k)
	   go to 85
70         fc(1)=1.0d-100
	   fc(2)=(a-4.0d0)/q*fc(1)
	   sp=0.0d0
	   u=fc(1)
	   f1=fc(2)
	   do 75 i=2,kb-1
	      v=u
	      u=f1
	      f1=(a-4.0d0*i*i)*u/q-v
	      if (i.ne.kb-1) then
		 fc(i+1)=f1
		 sp=sp+f1*f1
	      else
		 f2=f1
	      endif
75         continue
	   sp=sp+fc(1)**2+fc(2)**2
	   ss=s+sp*(f3/f2)**2
	   s0=1.0d0/dsqrt(ss)
	   do 80 j=1,km
	      if (j.lt.kb) fc(j)=s0*fc(j)*f3/f2
	      if (j.ge.kb) fc(j)=s0*fc(j)
80         continue
	endif
85      if (fc(1).lt.0.0d0) then
	   do 90 j=1,km
90            fc(j)=-fc(j)
	endif
	return
	end


	subroutine cva2(kd,m,q,a)
c
c       ======================================================
c       purpose: calculate a specific characteristic value of
c                mathieu functions
c       input :  m  --- order of mathieu functions
c                q  --- parameter of mathieu functions
c                kd --- case code
c                       kd=1 for cem(x,q)  ( m = 0,2,4,...)
c                       kd=2 for cem(x,q)  ( m = 1,3,5,...)
c                       kd=3 for sem(x,q)  ( m = 1,3,5,...)
c                       kd=4 for sem(x,q)  ( m = 2,4,6,...)
c       output:  a  --- characteristic value
c       routines called:
c             (1) refine for finding accurate characteristic
c                 values using an iteration method
c             (2) cv0 for finding initial characteristic
c                 values using polynomial approximation
c             (3) cvqm for computing initial characteristic
c                 values for q ó 3*m
c             (3) cvql for computing initial characteristic
c                 values for q ò m*m
c       ======================================================
c
	implicit double precision (a-h,o-z)
	if (m.le.12.or.q.le.3.0*m.or.q.gt.m*m) then
	    call cv0(kd,m,q,a)
	    if (q.ne.0.0d0) call refine(kd,m,q,a,1)
	else
	   ndiv=10
	   delta=(m-3.0)*m/ndiv
	   if ((q-3.0*m).le.(m*m-q)) then
5             nn=int((q-3.0*m)/delta)+1
	      delta=(q-3.0*m)/nn
	      q1=2.0*m
	      call cvqm(m,q1,a1)
	      q2=3.0*m
	      call cvqm(m,q2,a2)
	      qq=3.0*m
	      do 10 i=1,nn
		 qq=qq+delta
		 a=(a1*q2-a2*q1+(a2-a1)*qq)/(q2-q1)
		 iflag=1
		 if (i.eq.nn) iflag=-1
		 call refine(kd,m,qq,a,iflag)
		 q1=q2
		 q2=qq
		 a1=a2
		 a2=a
10            continue
	      if (iflag.eq.-10) then
		 ndiv=ndiv*2
		 delta=(m-3.0)*m/ndiv
		 go to 5
	      endif
	   else
15            nn=int((m*m-q)/delta)+1
	      delta=(m*m-q)/nn
	      q1=m*(m-1.0)
	      call cvql(kd,m,q1,a1)
	      q2=m*m
	      call cvql(kd,m,q2,a2)
	      qq=m*m
	      do 20 i=1,nn
		 qq=qq-delta
		 a=(a1*q2-a2*q1+(a2-a1)*qq)/(q2-q1)
		 iflag=1
		 if (i.eq.nn) iflag=-1
		 call refine(kd,m,qq,a,iflag)
		 q1=q2
		 q2=qq
		 a1=a2
		 a2=a
20            continue
	      if (iflag.eq.-10) then
		 ndiv=ndiv*2
		 delta=(m-3.0)*m/ndiv
		 go to 15
	      endif
	   endif
	endif
	return
	end


	subroutine refine(kd,m,q,a,iflag)
c
c       =====================================================
c       purpose: calculate the accurate characteristic value
c                by the secant method
c       input :  m --- order of mathieu functions
c                q --- parameter of mathieu functions
c                a --- initial characteristic value
c       output:  a --- refineed characteristic value
c       routine called:  cvf for computing the value of f for
c                        characteristic equation
c       ========================================================
c
	implicit double precision (a-h,o-z)
	eps=1.0d-14
	mj=10+m
	ca=a
	delta=0.0d0
	x0=a
	call cvf(kd,m,q,x0,mj,f0)
	x1=1.002*a
	call cvf(kd,m,q,x1,mj,f1)
5       do 10 it=1,100
	   mj=mj+1
	   x=x1-(x1-x0)/(1.0d0-f0/f1)
	   call cvf(kd,m,q,x,mj,f)
	   if (abs(1.0-x1/x).lt.eps.or.f.eq.0.0) go to 15
	   x0=x1
	   f0=f1
	   x1=x
10         f1=f
15      a=x
	if (delta.gt.0.05) then
	   a=ca
	   if (iflag.lt.0) then
	      iflag=-10
	   endif
	   return
	endif
	if (abs((a-ca)/ca).gt.0.05) then
	   x0=ca
	   delta=delta+0.005d0
	   call cvf(kd,m,q,x0,mj,f0)
	   x1=(1.0d0+delta)*ca
	   call cvf(kd,m,q,x1,mj,f1)
	   go to 5
	endif
	return
	end


	subroutine cvf(kd,m,q,a,mj,f)
c
c       ======================================================
c       purpose: compute the value of f for characteristic
c                equation of mathieu functions
c       input :  m --- order of mathieu functions
c                q --- parameter of mathieu functions
c                a --- characteristic value
c       output:  f --- value of f for characteristic equation
c       ======================================================
c
	implicit double precision (a-h,o-z)
	b=a
	ic=int(m/2)
	l=0
	l0=0
	j0=2
	jf=ic
	if (kd.eq.1) l0=2
	if (kd.eq.1) j0=3
	if (kd.eq.2.or.kd.eq.3) l=1
	if (kd.eq.4) jf=ic-1
	t1=0.0d0
	do 10 j=mj,ic+1,-1
10         t1=-q*q/((2.0d0*j+l)**2-b+t1)
	if (m.le.2) then
	   t2=0.0d0
	   if (kd.eq.1.and.m.eq.0) t1=t1+t1
	   if (kd.eq.1.and.m.eq.2) t1=-2.0*q*q/(4.0-b+t1)-4.0
	   if (kd.eq.2.and.m.eq.1) t1=t1+q
	   if (kd.eq.3.and.m.eq.1) t1=t1-q
	else
	   if (kd.eq.1) t0=4.0d0-b+2.0d0*q*q/b
	   if (kd.eq.2) t0=1.0d0-b+q
	   if (kd.eq.3) t0=1.0d0-b-q
	   if (kd.eq.4) t0=4.0d0-b
	   t2=-q*q/t0
	   do 15 j=j0,jf
15            t2=-q*q/((2.0d0*j-l-l0)**2-b+t2)
	endif
	f=(2.0d0*ic+l)**2+t1+t2-b
	return
	end


	subroutine cv0(kd,m,q,a0)
c
c       =====================================================
c       purpose: compute the initial characteristic value of
c                mathieu functions for m ó 12  or q ó 300 or
c                q ò m*m
c       input :  m  --- order of mathieu functions
c                q  --- parameter of mathieu functions
c       output:  a0 --- characteristic value
c       routines called:
c             (1) cvqm for computing initial characteristic
c                 value for q ó 3*m
c             (2) cvql for computing initial characteristic
c                 value for q ò m*m
c       ====================================================
c
	implicit double precision (a-h,o-z)
	q2=q*q
	if (m.eq.0) then
	   if (q.le.1.0) then
	      a0=(((.0036392*q2-.0125868)*q2+.0546875)*q2-.5)*q2
	   else if (q.le.10.0) then
	      a0=((3.999267d-3*q-9.638957d-2)*q-.88297)*q
     &           +.5542818
	   else
	      call cvql(kd,m,q,a0)
	   endif
	else if (m.eq.1) then
	   if (q.le.1.0.and.kd.eq.2) then
	      a0=(((-6.51e-4*q-.015625)*q-.125)*q+1.0)*q+1.0
	   else if (q.le.1.0.and.kd.eq.3) then
	      a0=(((-6.51e-4*q+.015625)*q-.125)*q-1.0)*q+1.0
	   else if (q.le.10.0.and. kd.eq.2) then
	      a0=(((-4.94603d-4*q+1.92917d-2)*q-.3089229)
     &           *q+1.33372)*q+.811752
	   else if (q.le.10.0.and.kd.eq.3) then
	      a0=((1.971096d-3*q-5.482465d-2)*q-1.152218)
     &           *q+1.10427
	   else
	      call cvql(kd,m,q,a0)
	   endif
	else if (m.eq.2) then
	   if (q.le.1.0.and.kd.eq.1) then
	      a0=(((-.0036391*q2+.0125888)*q2-.0551939)*q2
     &           +.416667)*q2+4.0
	   else if (q.le.1.0.and.kd.eq.4) then
	      a0=(.0003617*q2-.0833333)*q2+4.0
	   else if (q.le.15.and.kd.eq.1) then
	      a0=(((3.200972d-4*q-8.667445d-3)*q
     &           -1.829032d-4)*q+.9919999)*q+3.3290504
	   else if (q.le.10.0.and.kd.eq.4) then
	      a0=((2.38446d-3*q-.08725329)*q-4.732542d-3)
     &           *q+4.00909
	   else
	      call cvql(kd,m,q,a0)
	   endif
	else if (m.eq.3) then
	   if (q.le.1.0.and.kd.eq.2) then
	      a0=((6.348e-4*q+.015625)*q+.0625)*q2+9.0
	   else if (q.le.1.0.and.kd.eq.3) then
	      a0=((6.348e-4*q-.015625)*q+.0625)*q2+9.0
	   else if (q.le.20.0.and.kd.eq.2) then
	      a0=(((3.035731d-4*q-1.453021d-2)*q
     &           +.19069602)*q-.1039356)*q+8.9449274
	   else if (q.le.15.0.and.kd.eq.3) then
	      a0=((9.369364d-5*q-.03569325)*q+.2689874)*q
     &           +8.771735
	   else
	      call cvql(kd,m,q,a0)
	   endif
	else if (m.eq.4) then
	   if (q.le.1.0.and.kd.eq.1) then
	      a0=((-2.1e-6*q2+5.012e-4)*q2+.0333333)*q2+16.0
	   else if (q.le.1.0.and.kd.eq.4) then
	      a0=((3.7e-6*q2-3.669e-4)*q2+.0333333)*q2+16.0
	   else if (q.le.25.0.and.kd.eq.1) then
	      a0=(((1.076676d-4*q-7.9684875d-3)*q
     &           +.17344854)*q-.5924058)*q+16.620847
	   else if (q.le.20.0.and.kd.eq.4) then
	      a0=((-7.08719d-4*q+3.8216144d-3)*q
     &           +.1907493)*q+15.744
	   else
	      call cvql(kd,m,q,a0)
	   endif
	else if (m.eq.5) then
	   if (q.le.1.0.and.kd.eq.2) then
	      a0=((6.8e-6*q+1.42e-5)*q2+.0208333)*q2+25.0
	   else if (q.le.1.0.and.kd.eq.3) then
	      a0=((-6.8e-6*q+1.42e-5)*q2+.0208333)*q2+25.0
	   else if (q.le.35.0.and.kd.eq.2) then
	      a0=(((2.238231d-5*q-2.983416d-3)*q
     &           +.10706975)*q-.600205)*q+25.93515
	   else if (q.le.25.0.and.kd.eq.3) then
	      a0=((-7.425364d-4*q+2.18225d-2)*q
     &           +4.16399d-2)*q+24.897
	   else
	      call cvql(kd,m,q,a0)
	   endif
	else if (m.eq.6) then
	   if (q.le.1.0) then
	      a0=(.4d-6*q2+.0142857)*q2+36.0
	   else if (q.le.40.0.and.kd.eq.1) then
	      a0=(((-1.66846d-5*q+4.80263d-4)*q
     &           +2.53998d-2)*q-.181233)*q+36.423
	   else if (q.le.35.0.and.kd.eq.4) then
	      a0=((-4.57146d-4*q+2.16609d-2)*q-2.349616d-2)*q
     &           +35.99251
	   else
	      call cvql(kd,m,q,a0)
	   endif
	else if (m.eq.7) then
	   if (q.le.10.0) then
	      call cvqm(m,q,a0)
	   else if (q.le.50.0.and.kd.eq.2) then
	      a0=(((-1.411114d-5*q+9.730514d-4)*q
     &           -3.097887d-3)*q+3.533597d-2)*q+49.0547
	   else if (q.le.40.0.and.kd.eq.3) then
	      a0=((-3.043872d-4*q+2.05511d-2)*q
     &           -9.16292d-2)*q+49.19035
	   else
	      call cvql(kd,m,q,a0)
	   endif
	else if (m.ge.8) then
	   if (q.le.3.*m) then
	      call cvqm(m,q,a0)
	   else if (q.gt.m*m) then
	      call cvql(kd,m,q,a0)
	   else
	      if (m.eq.8.and.kd.eq.1) then
		 a0=(((8.634308d-6*q-2.100289d-3)*q+.169072)*q
     &              -4.64336)*q+109.4211
	      else if (m.eq.8.and.kd.eq.4) then
		 a0=((-6.7842d-5*q+2.2057d-3)*q+.48296)*q+56.59
	      else if (m.eq.9.and.kd.eq.2) then
		 a0=(((2.906435d-6*q-1.019893d-3)*q+.1101965)*q
     &              -3.821851)*q+127.6098
	      else if (m.eq.9.and.kd.eq.3) then
		 a0=((-9.577289d-5*q+.01043839)*q+.06588934)*q
     &              +78.0198
	      else if (m.eq.10.and.kd.eq.1) then
		 a0=(((5.44927d-7*q-3.926119d-4)*q+.0612099)*q
     &              -2.600805)*q+138.1923
	      else if (m.eq.10.and.kd.eq.4) then
		 a0=((-7.660143d-5*q+.01132506)*q-.09746023)*q
     &              +99.29494
	      else if (m.eq.11.and.kd.eq.2) then
		 a0=(((-5.67615d-7*q+7.152722d-6)*q+.01920291)*q
     &              -1.081583)*q+140.88
	      else if (m.eq.11.and.kd.eq.3) then
		 a0=((-6.310551d-5*q+.0119247)*q-.2681195)*q
     &              +123.667
	      else if (m.eq.12.and.kd.eq.1) then
		 a0=(((-2.38351d-7*q-2.90139d-5)*q+.02023088)*q
     &              -1.289)*q+171.2723
	      else if (m.eq.12.and.kd.eq.4) then
		 a0=(((3.08902d-7*q-1.577869d-4)*q+.0247911)*q
     &              -1.05454)*q+161.471
	      endif
	   endif
	endif
	return
	end


	subroutine cvql(kd,m,q,a0)
c
c       ========================================================
c       purpose: compute the characteristic value of mathieu
c                functions  for q ò 3m
c       input :  m  --- order of mathieu functions
c                q  --- parameter of mathieu functions
c       output:  a0 --- initial characteristic value
c       ========================================================
c
	implicit double precision (a-h,o-z)
	if (kd.eq.1.or.kd.eq.2) w=2.0d0*m+1.0d0
	if (kd.eq.3.or.kd.eq.4) w=2.0d0*m-1.0d0
	w2=w*w
	w3=w*w2
	w4=w2*w2
	w6=w2*w4
	d1=5.0+34.0/w2+9.0/w4
	d2=(33.0+410.0/w2+405.0/w4)/w
	d3=(63.0+1260.0/w2+2943.0/w4+486.0/w6)/w2
	d4=(527.0+15617.0/w2+69001.0/w4+41607.0/w6)/w3
	c1=128.0
	p2=q/w4
	p1=dsqrt(p2)
	cv1=-2.0*q+2.0*w*dsqrt(q)-(w2+1.0)/8.0
	cv2=(w+3.0/w)+d1/(32.0*p1)+d2/(8.0*c1*p2)
	cv2=cv2+d3/(64.0*c1*p1*p2)+d4/(16.0*c1*c1*p2*p2)
	a0=cv1-cv2/(c1*p1)
	return
	end


	subroutine cvqm(m,q,a0)
c
c       =====================================================
c       purpose: compute the characteristic value of mathieu
c                functions for q ó m*m
c       input :  m  --- order of mathieu functions
c                q  --- parameter of mathieu functions
c       output:  a0 --- initial characteristic value
c       =====================================================
c
	implicit double precision (a-h,o-z)
	hm1=.5*q/(m*m-1.0)
	hm3=.25*hm1**3/(m*m-4.0)
	hm5=hm1*hm3*q/((m*m-1.0)*(m*m-9.0))
	a0=m*m+q*(hm1+(5.0*m*m+7.0)*hm3
     &     +(9.0*m**4+58.0*m*m+29.0)*hm5)
	return
	end


	subroutine jynb(n,x,nm,bj,dj,by,dy)
c
c       =====================================================
c       purpose: compute bessel functions jn(x), yn(x) and
c                their derivatives
c       input :  x --- argument of jn(x) and yn(x) ( x ò 0 )
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
