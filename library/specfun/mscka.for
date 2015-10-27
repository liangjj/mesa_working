	program mscka
c
c       ============================================================
c       purpose: this program computes the expansion coefficients 
c                of the prolate and oblate spheroidal functions, 
c                c2k, using subroutine scka
c       input :  m  --- mode parameter
c                n  --- mode parameter
c                c  --- spheroidal parameter
c                cv --- characteristic value
c                kd --- function code
c                       kd=1 for prolate; kd=-1 for oblate
c       output:  ck(k) --- expansion coefficients ck;
c                          ck(1), ck(2),... correspond to
c                          c0, c2,...
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
	dimension ck(200),eg(200)
	write(*,*)'please kd, m, n and c '
	read(*,*)kd,m,n,c
	call segv(m,n,c,kd,cv,eg)
	write(*,30)kd,m,n,c,cv
	call scka(m,n,c,cv,kd,ck)
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
30      format(1x,3hkd=,i3,',  ',2hm=,i3,',  ',2hn=,i3,',  ',2hc=,
     &         f5.1,',  ',4hcv =,f18.10)
	end


	subroutine scka(m,n,c,cv,kd,ck)
c
c       ======================================================
c       purpose: compute the expansion coefficients of the
c                prolate and oblate spheroidal functions, c2k
c       input :  m  --- mode parameter
c                n  --- mode parameter
c                c  --- spheroidal parameter
c                cv --- characteristic value
c                kd --- function code
c                       kd=1 for prolate; kd=-1 for oblate
c       output:  ck(k) --- expansion coefficients ck;
c                          ck(1), ck(2),... correspond to
c                          c0, c2,...
c       ======================================================
c
	implicit double precision (a-h,o-z)
	dimension ck(200)
	if (c.le.1.0d-10) c=1.0d-10
	nm=25+int((n-m)/2+c)
	cs=c*c*kd
	ip=1
	if (n-m.eq.2*int((n-m)/2)) ip=0
	fs=1.0d0
	f1=0.0d0
	f0=1.0d-100
	kb=0
	ck(nm+1)=0.0d0
	do 15 k=nm,1,-1
	   f=(((2.0d0*k+m+ip)*(2.0d0*k+m+1.0d0+ip)-cv+cs)*f0
     &       -4.0d0*(k+1.0d0)*(k+m+1.0d0)*f1)/cs
	   if (dabs(f).gt.dabs(ck(k+1))) then
	      ck(k)=f
	      f1=f0
	      f0=f
	      if (dabs(f).gt.1.0d+100) then
		 do 5 k1=nm,k,-1
5                   ck(k1)=ck(k1)*1.0d-100
		 f1=f1*1.0d-100
		 f0=f0*1.0d-100
	      endif
	   else
	      kb=k
	      fl=ck(k+1)
	      f1=1.0d0
	      f2=0.25d0*((m+ip)*(m+ip+1.0)-cv+cs)/(m+1.0)*f1
	      ck(1)=f1
	      if (kb.eq.1) then
		 fs=f2
	      else if (kb.eq.2) then
		 ck(2)=f2
		 fs=0.125d0*(((m+ip+2.0)*(m+ip+3.0)-cv+cs)*f2
     &              -cs*f1)/(m+2.0)
	      else
		 ck(2)=f2
		 do 10 j=3,kb+1
		    f=0.25d0*(((2.0*j+m+ip-4.0)*(2.0*j+m+ip-
     &                3.0)-cv+cs)*f2-cs*f1)/((j-1.0)*(j+m-1.0))
		    if (j.le.kb) ck(j)=f
		    f1=f2
10                  f2=f
		 fs=f
	      endif
	      go to 20
	   endif
15      continue
20      su1=0.0d0
	do 25 k=1,kb
25         su1=su1+ck(k)
	su2=0.0d0
	do 30 k=kb+1,nm
30         su2=su2+ck(k)
	r1=1.0d0
	do 35 j=1,(n+m+ip)/2
35         r1=r1*(j+0.5d0*(n+m+ip))
	r2=1.0d0
	do 40 j=1,(n-m-ip)/2
40         r2=-r2*j
	if (kb.eq.0) then
	    s0=r1/(2.0d0**n*r2*su2)
	else
	    s0=r1/(2.0d0**n*r2*(fl/fs*su1+su2))
	endif
	do 45 k=1,kb
45         ck(k)=fl/fs*s0*ck(k)
	do 50 k=kb+1,nm
50         ck(k)=s0*ck(k)
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
