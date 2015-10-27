	program maswfb
c
c       ============================================================
c       purpose: this program computes the prolate and oblate 
c                spheroidal angular functions of the first kind
c                and their derivatives using subroutine aswfb
c       input :  m  --- mode parameter,  m = 0,1,2,...
c                n  --- mode parameter,  n = m,m+1,...
c                c  --- spheroidal parameter
c                x  --- argument of angular function, |x| ó 1.0
c                kd --- function code
c                       kd=1 for prolate;  kd=-1 for oblate
c                cv --- characteristic value
c       output:  s1f --- angular function of the first kind
c                s1d --- derivative of the angular function of
c                        the first kind
c       examples:
c               kd = 1, m = 2, n = 3, c = 3.0 and cv = 14.8277782138
c                  x         smn(c,x)            smn'(c,x)
c                --------------------------------------------
c                 0.2      .28261309d+01       .12418631d+02
c                 0.5      .49938554d+01       .92761604d+00
c                 0.8      .31693975d+01      -.12646552d+02
c
c               kd =-1, m = 2, n = 3, c = 3.0 and cv = 8.8093939208
c                  x         smn(-ic,x)         smn'(-ic,x)
c                --------------------------------------------
c                 0.2      .29417848d+01       .14106305d+02
c                 0.5      .64138827d+01       .76007194d+01
c                 0.8      .60069873d+01      -.14387479d+02
c       ============================================================
c
	implicit double precision (a-h,o-z)
	dimension eg(200)
	write(*,*)'please enter kd, m, n and c '
	read(*,*)kd,m,n,c
	write(*,10)kd,m,n,c
	call segv(m,n,c,kd,cv,eg)
	write(*,20)cv
	write(*,*)
	if (kd.eq.1 ) then
	   write(*,*)'    x         smn(c,x)            smn''(c,x)'
	else if (kd.eq.-1) then
	   write(*,*)'    x         smn(-ic,x)         smn''(-ic,x)'
	endif
	write(*,*)'  --------------------------------------------'
	do 5 i=0,20
	   x=-1.0d0+0.1d0*i
	   call aswfb(m,n,c,x,kd,cv,s1f,s1d)  
5          write(*,30)x,s1f,s1d
10      format(1x,'kd ='i2,', ','m =',i2,', ','n =',i2,', ','c =',f5.1)
20      format(1x,' cv =',f18.10)
30      format(1x,f5.1,2d20.8)
	end


	subroutine aswfb(m,n,c,x,kd,cv,s1f,s1d)
c
c       ===========================================================
c       purpose: compute the prolate and oblate spheroidal angular
c                functions of the first kind and their derivatives
c       input :  m  --- mode parameter,  m = 0,1,2,...
c                n  --- mode parameter,  n = m,m+1,...
c                c  --- spheroidal parameter
c                x  --- argument of angular function, |x| < 1.0
c                kd --- function code
c                       kd=1 for prolate;  kd=-1 for oblate
c                cv --- characteristic value
c       output:  s1f --- angular function of the first kind
c                s1d --- derivative of the angular function of
c                        the first kind
c       routines called:
c            (1) sdmn for computing expansion coefficients dk
c            (2) lpmns for computing associated legendre function
c                of the first kind pmn(x)
c       ===========================================================
c
	implicit double precision (a-h,o-z)
	dimension df(200),pm(0:251),pd(0:251)
	eps=1.0d-14
	ip=1
	if (n-m.eq.2*int((n-m)/2)) ip=0
	nm=25+int((n-m)/2+c)
	nm2=2*nm+m
	call sdmn(m,n,c,cv,kd,df)
	call lpmns(m,nm2,x,pm,pd)
	su1=0.0d0
	do 10 k=1,nm
	   mk=m+2*(k-1)+ip
	   su1=su1+df(k)*pm(mk)
	   if (dabs(sw-su1).lt.dabs(su1)*eps) goto 15
10         sw=su1
15      s1f=(-1)**m*su1
	su1=0.0d0
	do 20 k=1,nm
	   mk=m+2*(k-1)+ip
	   su1=su1+df(k)*pd(mk)
	   if (dabs(sw-su1).lt.dabs(su1)*eps) goto 25
20         sw=su1
25      s1d=(-1)**m*su1
	return
	end


	subroutine sdmn(m,n,c,cv,kd,df)
c
c       =====================================================
c       purpose: compute the expansion coefficients of the
c                prolate and oblate spheroidal functions, dk
c       input :  m  --- mode parameter
c                n  --- mode parameter
c                c  --- spheroidal parameter
c                cv --- characteristic value
c                kd --- function code
c                       kd=1 for prolate; kd=-1 for oblate
c       output:  df(k) --- expansion coefficients dk;
c                          df(1), df(2), ... correspond to
c                          d0, d2, ... for even n-m and d1,
c                          d3, ... for odd n-m
c       =====================================================
c
	implicit double precision (a-h,o-z)
	dimension a(200),d(200),g(200),df(200)
	nm=25+int(0.5*(n-m)+c)
	if (c.lt.1.0d-10) then
	   do 5 i=1,nm
5             df(i)=0d0
	   df((n-m)/2+1)=1.0d0
	   return
	endif   
	cs=c*c*kd
	ip=1
	if (n-m.eq.2*int((n-m)/2)) ip=0
	do 10 i=1,nm+2
	   if (ip.eq.0) k=2*(i-1)
	   if (ip.eq.1) k=2*i-1
	   dk0=m+k
	   dk1=m+k+1
	   dk2=2*(m+k)
	   d2k=2*m+k
	   a(i)=(d2k+2.0)*(d2k+1.0)/((dk2+3.0)*(dk2+5.0))*cs
	   d(i)=dk0*dk1+(2.0*dk0*dk1-2.0*m*m-1.0)/((dk2-1.0)
     &          *(dk2+3.0))*cs
	   g(i)=k*(k-1.0)/((dk2-3.0)*(dk2-1.0))*cs
10      continue
	fs=1.0d0
	f1=0.0d0
	f0=1.0d-100
	kb=0
	df(nm+1)=0.0d0
	do 30 k=nm,1,-1
	   f=-((d(k+1)-cv)*f0+a(k+1)*f1)/g(k+1)
	   if (dabs(f).gt.dabs(df(k+1))) then
	      df(k)=f
	      f1=f0
	      f0=f
	      if (dabs(f).gt.1.0d+100) then
		 do 12 k1=k,nm
12                  df(k1)=df(k1)*1.0d-100
		 f1=f1*1.0d-100
		 f0=f0*1.0d-100
	      endif  
	   else
	      kb=k
	      fl=df(k+1)
	      f1=1.0d-100
	      f2=-(d(1)-cv)/a(1)*f1
	      df(1)=f1
	      if (kb.eq.1) then
		 fs=f2
	      else if (kb.eq.2) then
		 df(2)=f2
		 fs=-((d(2)-cv)*f2+g(2)*f1)/a(2)
	      else 
		 df(2)=f2
		 do 20 j=3,kb+1
		    f=-((d(j-1)-cv)*f2+g(j-1)*f1)/a(j-1)
		    if (j.le.kb) df(j)=f
		    if (dabs(f).gt.1.0d+100) then
		       do 15 k1=1,j
15                        df(k1)=df(k1)*1.0d-100
		       f=f*1.0d-100
		       f2=f2*1.0d-100
		    endif  
		    f1=f2
20                  f2=f
		 fs=f
	      endif
	      go to 35
	   endif
30      continue
35      su1=0.0d0
	r1=1.0d0
	do 40 j=m+ip+1,2*(m+ip)
40         r1=r1*j
	su1=df(1)*r1
	do 45 k=2,kb
	   r1=-r1*(k+m+ip-1.5d0)/(k-1.0d0)
45           su1=su1+r1*df(k)
	su2=0.0d0
	do 50 k=kb+1,nm
	   if (k.ne.1) r1=-r1*(k+m+ip-1.5d0)/(k-1.0d0)
	   su2=su2+r1*df(k)
	   if (dabs(sw-su2).lt.dabs(su2)*1.0d-14) goto 55
50         sw=su2
55      r3=1.0d0
	do 60 j=1,(m+n+ip)/2
60         r3=r3*(j+0.5d0*(n+m+ip))
	r4=1.0d0
	do 65 j=1,(n-m-ip)/2
65         r4=-4.0d0*r4*j
	s0=r3/(fl*(su1/fs)+su2)/r4
	do 70 k=1,kb
70         df(k)=fl/fs*s0*df(k)
	do 75 k=kb+1,nm
75         df(k)=s0*df(k)
	return
	end


	subroutine segv(m,n,c,kd,cv,eg)
c
c       =========================================================
c       purpose: compute the characteristic values of spheroidal
c                wave functions
c       input :  m  --- mode parameter
c                n  --- mode parameter
c                c  --- spheroidal parameter
c                kd --- function code
c                       kd=1 for prolate; kd=-1 for oblate
c       output:  cv --- characteristic value for given m, n and c
c                eg(l) --- characteristic value for mode m and n'
c                          ( l = n' - m + 1 )
c       =========================================================
c
	implicit double precision (a-h,o-z)
	dimension b(100),h(100),d(300),e(300),f(300),cv0(100),
     &            a(300),g(300),eg(200)
	if (c.lt.1.0d-10) then
	   do 5 i=1,n
5             eg(i)=(i+m)*(i+m-1.0d0)
	   go to 70
	endif                                           
	icm=(n-m+2)/2
	nm=10+int(0.5*(n-m)+c)
	cs=c*c*kd
	do 60 l=0,1
	   do 10 i=1,nm
	      if (l.eq.0) k=2*(i-1)
	      if (l.eq.1) k=2*i-1
	      dk0=m+k
	      dk1=m+k+1
	      dk2=2*(m+k)
	      d2k=2*m+k
	      a(i)=(d2k+2.0)*(d2k+1.0)/((dk2+3.0)*(dk2+5.0))*cs
	      d(i)=dk0*dk1+(2.0*dk0*dk1-2.0*m*m-1.0)/((dk2-1.0)
     &             *(dk2+3.0))*cs
10            g(i)=k*(k-1.0)/((dk2-3.0)*(dk2-1.0))*cs
	   do 15 k=2,nm
	      e(k)=dsqrt(a(k-1)*g(k))
15            f(k)=e(k)*e(k)
	   f(1)=0.0d0
	   e(1)=0.0d0
	   xa=d(nm)+dabs(e(nm))
	   xb=d(nm)-dabs(e(nm))
	   nm1=nm-1
	   do 20 i=1,nm1
	      t=dabs(e(i))+dabs(e(i+1))
	      t1=d(i)+t
	      if (xa.lt.t1) xa=t1
	      t1=d(i)-t
	      if (t1.lt.xb) xb=t1
20         continue
	   do 25 i=1,icm
	      b(i)=xa
25            h(i)=xb
	   do 55 k=1,icm
	      do 30 k1=k,icm
		 if (b(k1).lt.b(k)) then
		    b(k)=b(k1)
		    go to 35
		 endif
30            continue
35            if (k.ne.1.and.h(k).lt.h(k-1)) h(k)=h(k-1)
40            x1=(b(k)+h(k))/2.0d0
	      cv0(k)=x1
	      if (dabs((b(k)-h(k))/x1).lt.1.0d-14) go to 50
	      j=0
	      s=1.0d0
	      do 45 i=1,nm
		 if (s.eq.0.0d0) s=s+1.0d-30
		 t=f(i)/s
		 s=d(i)-t-x1
		 if (s.lt.0.0d0) j=j+1
45            continue
	      if (j.lt.k) then
		 h(k)=x1
	      else
		 b(k)=x1
		 if (j.ge.icm) then
		    b(icm)=x1
		 else
		    if (h(j+1).lt.x1) h(j+1)=x1
		    if (x1.lt.b(j)) b(j)=x1
		 endif
	      endif
	      go to 40
50            cv0(k)=x1
	      if (l.eq.0) eg(2*k-1)=cv0(k)
	      if (l.eq.1) eg(2*k)=cv0(k)
55         continue
60      continue
70      cv=eg(n-m+1)
	return
	end


	subroutine lpmns(m,n,x,pm,pd)
c
c       ========================================================
c       purpose: compute associated legendre functions pmn(x)
c                and pmn'(x) for a given order
c       input :  x --- argument of pmn(x)
c                m --- order of pmn(x),  m = 0,1,2,...,n
c                n --- degree of pmn(x), n = 0,1,2,...,n
c       output:  pm(n) --- pmn(x)
c                pd(n) --- pmn'(x)
c       ========================================================
c
	implicit double precision (a-h,o-z)
	dimension pm(0:n),pd(0:n)
	do 10 k=0,n
	   pm(k)=0.0d0
10         pd(k)=0.0d0
	if (dabs(x).eq.1.0d0) then
	   do 15 k=0,n
	      if (m.eq.0) then
		 pm(k)=1.0d0
		 pd(k)=0.5d0*k*(k+1.0)
		 if (x.lt.0.0) then
		    pm(k)=(-1)**k*pm(k)
		    pd(k)=(-1)**(k+1)*pd(k)
		 endif
	      else if (m.eq.1) then
		 pd(k)=1.0d+300
	      else if (m.eq.2) then
		 pd(k)=-0.25d0*(k+2.0)*(k+1.0)*k*(k-1.0)
		 if (x.lt.0.0) pd(k)=(-1)**(k+1)*pd(k)
	      endif
15         continue
	   return
	endif
	x0=dabs(1.0d0-x*x)
	pm0=1.0d0
	pmk=pm0
	do 20 k=1,m
	   pmk=(2.0d0*k-1.0d0)*dsqrt(x0)*pm0
20         pm0=pmk
	pm1=(2.0d0*m+1.0d0)*x*pm0
	pm(m)=pmk
	pm(m+1)=pm1
	do 25 k=m+2,n
	   pm2=((2.0d0*k-1.0d0)*x*pm1-(k+m-1.0d0)*pmk)/(k-m)
	   pm(k)=pm2
	   pmk=pm1
25         pm1=pm2
	pd(0)=((1.0d0-m)*pm(1)-x*pm(0))/(x*x-1.0)  
	do 30 k=1,n
30         pd(k)=(k*x*pm(k)-(k+m)*pm(k-1))/(x*x-1.0d0)
	return
	end
