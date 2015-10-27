	program msckb
c
c       ============================================================
c       purpose: this program computes the expansion coefficients 
c                of the prolate and oblate spheroidal functions, 
c                c2k, using subroutine sckb
c       input :  m  --- mode parameter
c                n  --- mode parameter
c                c  --- spheroidal parameter
c                cv --- characteristic value
c                kd --- function code
c                       kd=1 for prolate; kd=-1 for oblate
c       output:  ck(k) --- expansion coefficients ck;
c                          ck(1), ck(2), ... correspond to
c                          c0, c2, ...
c       example: compute the first 13 expansion coefficients c2k for
c                kd= 1, m=2, n=3, c=3.0 and cv=14.8277782138; and
c                kd=-1, m=2, n=3, c=3.0 and cv=8.80939392077
c
c                coefficients of prolate and oblate functions
c
c                  k         c2k(c)              c2k(-ic)
c                ---------------------------------------------
c                  0     .9173213327d+01     .2489664942d+02
c                  1     .4718258929d+01    -.1205287032d+02
c                  2     .9841212916d+00     .2410564082d+01
c                  3     .1151870224d+00    -.2735821590d+00
c                  4     .8733916403d-02     .2026057157d-01
c                  5     .4663888254d-03    -.1061946315d-02
c                  6     .1853910398d-04     .4158091152d-04
c                  7     .5708084895d-06    -.1264400411d-05
c                  8     .1402786472d-07     .3074963448d-07
c                  9     .2817194508d-09    -.6120579463d-09
c                 10     .4712094447d-11     .1015900041d-10
c                 11     .6667838485d-13    -.1427953361d-12
c                 12     .8087995432d-15     .1721924955d-14
c       ============================================================
c
	implicit double precision (a-h,o-z)
	dimension ck(200),df(200),eg(200)
	write(*,*)'please kd, m, n and c '
	read(*,*)kd,m,n,c
	call segv(m,n,c,kd,cv,eg)
	write(*,30)kd,m,n,c,cv
	call sdmn(m,n,c,cv,kd,df)
	call sckb(m,n,c,df,ck)
	write(*,*)
	if (kd.eq.1) then
	   write(*,*)'coefficients of prolate function'
	   write(*,*)
	   write(*,*)'   k            c2k(c)'
	else
	   write(*,*)'coefficients of oblate function'
	   write(*,*)
	   write(*,*)'   k           c2k(-ic)'
	endif
	write(*,*)'----------------------------'
	nm=25+int((n-m)/2+c)
	do 10 k=1,nm
10         write(*,20)k-1,ck(k)
20      format(2x,i3,4x,d18.10)
30      format(1x,3hkd=,i3,',  ',2hm=,i3,',  ',2hn=,i3,
     &         ',  ',2hc=,f5.1,',  ',4hcv =,f18.10)
	end


	subroutine sckb(m,n,c,df,ck)
c
c       ======================================================
c       purpose: compute the expansion coefficients of the
c                prolate and oblate spheroidal functions, c2k
c       input :  m  --- mode parameter
c                n  --- mode parameter
c                c  --- spheroidal parameter
c                df(k) --- expansion coefficients dk
c       output:  ck(k) --- expansion coefficients ck;
c                          ck(1), ck(2), ... correspond to
c                          c0, c2, ...
c       ======================================================
c
	implicit double precision (a-h,o-z)
	dimension df(200),ck(200)
	if (c.le.1.0d-10) c=1.0d-10
	nm=25+int(0.5*(n-m)+c)
	ip=1
	if (n-m.eq.2*int((n-m)/2)) ip=0
	fac=-0.5d0**m
	reg=1.0d0
	if (m+nm.gt.80) reg=1.0d-200
	do 35 k=0,nm-1
	   fac=-fac
	   i1=2*k+ip+1
	   r=reg
	   do 10 i=i1,i1+2*m-1
10            r=r*i
	   i2=k+m+ip
	   do 15 i=i2,i2+k-1
15            r=r*(i+0.5d0)
	   sum=r*df(k+1)
	   do 20 i=k+1,nm
	      d1=2.0d0*i+ip
	      d2=2.0d0*m+d1
	      d3=i+m+ip-0.5d0
	      r=r*d2*(d2-1.0d0)*i*(d3+k)/(d1*(d1-1.0d0)*(i-k)*d3)
	      sum=sum+r*df(i+1)
	      if (dabs(sw-sum).lt.dabs(sum)*1.0d-14) go to 25
20            sw=sum
25         r1=reg
	   do 30 i=2,m+k
30            r1=r1*i
35         ck(k+1)=fac*sum/r1
	return
	end 


	subroutine sdmn(m,n,c,cv,kd,df)
c
c       =====================================================
c       purpose: compute the expansion coefficients of the
c                prolate and oblate spheroidal functions
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
