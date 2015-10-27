	program msdmn
c
c       ===========================================================
c       purpose: this program computes the expansion coefficients 
c                of the prolate and oblate spheroidal functions, 
c                dk, using subroutine sdmn
c       input :  m  --- mode parameter
c                n  --- mode parameter
c                c  --- spheroidal parameter
c                cv --- characteristic value
c                kd --- function code
c                       kd=1 for prolate; kd=-1 for oblate
c       output:  df(k) --- expansion coefficients dk;
c                          df(1), df(2),... correspond to
c                          d0, d2,... for even n-m and d1,
c                          d3,... for odd n-m
c       example: compute the first 12 expansion coefficients for
c                kd= 1, m=2, n=2, c=3.0 and cv=7.1511005241; and
c                kd=-1, m=2, n=2, c=3.0 and cv=4.5264604622
c
c                coefficients of prolate and oblate functions
c
c                  r          dr(c)             dr(-ic)
c                -------------------------------------------
c                  0     .9237882817d+00    .1115434000d+01
c                  2    -.2901607696d-01    .4888489020d-01
c                  4     .8142246173d-03    .1600845667d-02
c                  6    -.1632270292d-04    .3509183384d-04
c                  8     .2376699010d-06    .5416293446d-06
c                 10    -.2601391701d-08    .6176624069d-08
c                 12     .2209142844d-10    .5407431236d-10
c                 14    -.1494812074d-12    .3745889118d-12
c                 16     .8239302207d-15    .2103624480d-14
c                 18    -.3768260778d-17    .9768323113d-17
c                 20     .1452384658d-19    .3812753620d-19
c                 22    -.4780280430d-22    .1268321726d-21
c       ===========================================================
c
	implicit double precision (a-h,o-z)
	dimension df(200),eg(200)
	write(*,*)'please enter kd, m, n and c '
	read(*,*)kd,m,n,c
	call segv(m,n,c,kd,cv,eg)
	write(*,30)kd,m,n,c,cv
	call sdmn(m,n,c,cv,kd,df)
	write(*,*)
	if (kd.eq.1) then
	   write(*,*)'coefficients of prolate function'
	   write(*,*)
	   write(*,*)'   r             dr(c)'
	else
	   write(*,*)'coefficients of oblate function'
	   write(*,*)
	   write(*,*)'   r            dr(-ic)'
	endif
	write(*,*)'----------------------------'
	nm=25+int(0.5*(n-m)+c)         
	do 10 k=1,nm 
	   if (n-m.eq.2*int((n-m)/2)) then
	      j=2*(k-1)
	   else
	      j=2*k-1
	   endif
10         write(*,20)j,df(k)
20      format(2x,i3,4x,d18.10)
30      format(1x,3hkd=,i3,',  ',2hm=,i3,',  ',2hn=,i3,',  ',2hc=,
     &         f5.1,',  ',4hcv =,f18.10)
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
c                          df(1), df(2),... correspond to
c                          d0, d2,... for even n-m and d1,
c                          d3,... for odd n-m
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
