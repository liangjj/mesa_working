	program mrswfo
c
c       =============================================================
c       purpose: this program computes the radial oblate spheriodal 
c                functions of the first and second kinds, and their 
c                derivatives using subroutine rswfo
c       input :  m  --- mode parameter, m = 0,1,2,...
c                n  --- mode parameter, n = m,m+1,m+2,...
c                c  --- spheroidal parameter
c                x  --- argument (x ע 0)
c                cv --- characteristic value
c                kf --- function code
c                       kf=1 for the first kind
c                       kf=2 for the second kind
c                       kf=3 for both the first and second kinds
c       output:  r1f --- radial function of the first kind
c                r1d --- derivative of the radial function of
c                        the first kind
c                r2f --- radial function of the second kind
c                r2d --- derivative of the radial function of
c                        the second kind
c       example:
c                kd=-1, m = 2, n = 3, c = 5.0 and cv = 2.10980581604
c
c    x  r23(1)(-ic,ix) r23(1)'(-ic,ix) r23(2)(-ic,ix) r23(2)'(-ic,ix)
c  ------------------------------------------------------------------
c   0.0      0.0000000 (-1) 4.9911346  (-1)-4.0071049  (-1) 3.3065682
c   0.5 (-1) 2.0215966 (-1) 1.9559661  (-1)-1.5154835  (-1) 6.4482526
c   1.0 (-1) 1.0526737 (-1)-5.4010624  (-1) 1.3065731  (-1) 2.7958491
c   1.5 (-1)-1.1611246 (-2)-8.9414690  (-2) 3.7405936  (-1)-5.0118500
c   5.0 (-2) 3.6463251 (-2)-8.3283065  (-2) 1.5549615  (-1) 1.7544481
c       ==============================================================
c
	implicit double precision (a-h,o-z)
	dimension eg(200)
	write(*,*)'please enter kf'
	read(*,*)kf
	write(*,10)kf
	write(*,*)'please enter m, n, c and x ( x ע 0 )'
	read(*,*)m,n,c,x
	call segv(m,n,c,-1,cv,eg)
	write(*,20)m,n,c,cv,x
	write(*,*)
	call rswfo(m,n,c,x,cv,kf,r1f,r1d,r2f,r2d)
	write(*,*)'   x    rmn(1)(-ic,ix)  rmn''(1)(-ic,ix)  ',
     &            ' rmn(2)(-ic,ix)  rmn''(2)(-ic,ix)'
	write(*,*)'  ---------------------------------------------',
     &            '---------------------------'
	if (kf.eq.1) then
	   write(*,30)x,r1f,r1d
	else if (kf.eq.2) then
	   write(*,40)x,r2f,r2d
	else if (kf.eq.3) then
	   write(*,30)x,r1f,r1d,r2f,r2d
	endif
	if (kf.eq.3) then
	   write(*,50)r1f*r2d-r2f*r1d,1.0d0/(c*(x*x+1.0d0))
	   write(*,60)
	endif
10      format(1x,3hkf=,i3)
20      format(1x,2hm=,i2,',   ',2hn=,i2,',   ',2hc=,f5.1,
     &         ',    ',4hcv =,f18.10,',  ',2hx=,f5.2)
30      format(1x,f5.1,4d17.8)
40      format(1x,f5.1,34x,4d17.8)
50      format(1x,/1x,'wronskian check:',/1x,'computed value =',
     &         d17.8,5x,'exact value =',d17.8)
60      format(1x,/1x,'caution: this check is not accurate if it ',
     &        'involves',/1x,'         the subtraction of two ',
     &        'similar numbers')
	end


	subroutine rswfo(m,n,c,x,cv,kf,r1f,r1d,r2f,r2d)
c
c       ==========================================================
c       purpose: compute oblate radial functions of the first
c                and second kinds, and their derivatives
c       input :  m  --- mode parameter,  m = 0,1,2,...
c                n  --- mode parameter,  n = m,m+1,m+2,...
c                c  --- spheroidal parameter
c                x  --- argument (x ע 0)
c                cv --- characteristic value
c                kf --- function code
c                       kf=1 for the first kind
c                       kf=2 for the second kind
c                       kf=3 for both the first and second kinds
c       output:  r1f --- radial function of the first kind
c                r1d --- derivative of the radial function of
c                        the first kind
c                r2f --- radial function of the second kind
c                r2d --- derivative of the radial function of
c                        the second kind
c       routines called:
c            (1) sdmn for computing expansion coefficients dk
c            (2) rmn1 for computing prolate or oblate radial
c                function of the first kind
c            (3) rmn2l for computing prolate or oblate radial
c                function of the second kind for a large argument
c            (4) rmn2so for computing oblate radial functions of
c                the second kind for a small argument
c       ==========================================================
c
	implicit double precision (a-h,o-z)
	dimension df(200)
	kd=-1
	call sdmn(m,n,c,cv,kd,df)
	if (kf.ne.2) then
	   call rmn1(m,n,c,x,df,kd,r1f,r1d)
	endif
	if (kf.gt.1) then
	   id=10
	   if (x.gt.1.0d-8) then
	      call rmn2l(m,n,c,x,df,kd,r2f,r2d,id)
	   endif
	   if (id.gt.-1) then
	      call rmn2so(m,n,c,x,cv,df,kd,r2f,r2d)
	   endif
	endif
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


	subroutine rmn1(m,n,c,x,df,kd,r1f,r1d)
c
c       =======================================================
c       purpose: compute prolate and oblate spheroidal radial
c                functions of the first kind for given m, n,
c                c and x
c       routines called:
c            (1) sckb for computing expansion coefficients c2k
c            (2) sphj for computing the spherical bessel
c                functions of the first kind     
c       =======================================================
c
	implicit double precision (a-h,o-z)
	dimension ck(200),df(200),sj(0:251),dj(0:251)
	eps=1.0d-14
	ip=1
	nm1=int((n-m)/2)
	if (n-m.eq.2*nm1) ip=0
	nm=25+nm1+int(c)
	reg=1.0d0
	if (m+nm.gt.80) reg=1.0d-200
	r0=reg
	do 10 j=1,2*m+ip
10         r0=r0*j
	r=r0    
	suc=r*df(1)
	do 15 k=2,nm
	   r=r*(m+k-1.0)*(m+k+ip-1.5d0)/(k-1.0d0)/(k+ip-1.5d0)
	   suc=suc+r*df(k)
	   if (k.gt.nm1.and.dabs(suc-sw).lt.dabs(suc)*eps) go to 20
15         sw=suc
20      continue
	if (x.eq.0.0) then
	   call sckb(m,n,c,df,ck)
	   sum=0.0d0
	   do 25 j=1,nm
	      sum=sum+ck(j)
	      if (dabs(sum-sw1).lt.dabs(sum)*eps) go to 30
25            sw1=sum
30         r1=1.0d0
	   do 35 j=1,(n+m+ip)/2
35            r1=r1*(j+0.5d0*(n+m+ip))
	   r2=1.0d0
	   do 40 j=1,m
40            r2=2.0d0*c*r2*j
	   r3=1.0d0
	   do 45 j=1,(n-m-ip)/2
45            r3=r3*j
	   sa0=(2.0*(m+ip)+1.0)*r1/(2.0**n*c**ip*r2*r3)
	   if (ip.eq.0) then
	      r1f=sum/(sa0*suc)*df(1)*reg
	      r1d=0.0d0
	   else if (ip.eq.1) then
	      r1f=0.0d0
	      r1d=sum/(sa0*suc)*df(1)*reg
	   endif
	   return
	endif
	cx=c*x
	nm2=2*nm+m
	call sphj(nm2,cx,nm2,sj,dj)
	a0=(1.0d0-kd/(x*x))**(0.5d0*m)/suc  
	r1f=0.0d0
	do 50 k=1,nm
	   l=2*k+m-n-2+ip
	   if (l.eq.4*int(l/4)) lg=1
	   if (l.ne.4*int(l/4)) lg=-1
	   if (k.eq.1) then
	      r=r0
	   else
	      r=r*(m+k-1.0)*(m+k+ip-1.5d0)/(k-1.0d0)/(k+ip-1.5d0)
	   endif
	   np=m+2*k-2+ip
	   r1f=r1f+lg*r*df(k)*sj(np)
	   if (k.gt.nm1.and.dabs(r1f-sw).lt.dabs(r1f)*eps) go to 55
50         sw=r1f
55      r1f=r1f*a0
	b0=kd*m/x**3.0d0/(1.0-kd/(x*x))*r1f    
	sud=0.0d0
	do 60 k=1,nm
	   l=2*k+m-n-2+ip
	   if (l.eq.4*int(l/4)) lg=1
	   if (l.ne.4*int(l/4)) lg=-1
	   if (k.eq.1) then
	      r=r0
	   else
	      r=r*(m+k-1.0)*(m+k+ip-1.5d0)/(k-1.0d0)/(k+ip-1.5d0)
	   endif
	   np=m+2*k-2+ip
	   sud=sud+lg*r*df(k)*dj(np)
	   if (k.gt.nm1.and.dabs(sud-sw).lt.dabs(sud)*eps) go to 65
60         sw=sud
65      r1d=b0+a0*c*sud
	return
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
	reg=1.0d0
	if (m+nm.gt.80) reg=1.0d-200
	fac=-0.5d0**m
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
	      if (dabs(sw-sum).lt.dabs(sum)*1.0d-14) goto 25
20            sw=sum
25         r1=reg
	   do 30 i=2,m+k
30            r1=r1*i
35         ck(k+1)=fac*sum/r1
	return
	end 


	subroutine sphj(n,x,nm,sj,dj)
c
c       =======================================================
c       purpose: compute spherical bessel functions jn(x) and
c                their derivatives
c       input :  x --- argument of jn(x)
c                n --- order of jn(x)  ( n = 0,1,תתת )
c       output:  sj(n) --- jn(x)
c                dj(n) --- jn'(x)
c                nm --- highest order computed
c       routines called:
c                msta1 and msta2 for computing the starting
c                point for backward recurrence
c       =======================================================
c
	implicit double precision (a-h,o-z)
	dimension sj(0:n),dj(0:n)
	nm=n
	if (dabs(x).eq.1.0d-100) then
	   do 10 k=0,n
	      sj(k)=0.0d0
10            dj(k)=0.0d0
	   sj(0)=1.0d0
	   dj(1)=.3333333333333333d0
	   return
	endif
	sj(0)=dsin(x)/x
	sj(1)=(sj(0)-dcos(x))/x
	if (n.ge.2) then
	   sa=sj(0)
	   sb=sj(1)
	   m=msta1(x,200)
	   if (m.lt.n) then
	      nm=m
	   else
	      m=msta2(x,n,15)
	   endif
	   f0=0.0d0
	   f1=1.0d0-100
	   do 15 k=m,0,-1
	      f=(2.0d0*k+3.0d0)*f1/x-f0
	      if (k.le.nm) sj(k)=f
	      f0=f1
15            f1=f
	   if (dabs(sa).gt.dabs(sb)) cs=sa/f
	   if (dabs(sa).le.dabs(sb)) cs=sb/f0
	   do 20 k=0,nm
20            sj(k)=cs*sj(k)
	endif      
	dj(0)=(dcos(x)-dsin(x)/x)/x
	do 25 k=1,nm
25         dj(k)=sj(k-1)-(k+1.0d0)*sj(k)/x
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


	subroutine rmn2l(m,n,c,x,df,kd,r2f,r2d,id)
c
c       ========================================================
c       purpose: compute prolate and oblate spheroidal radial
c                functions of the second kind for given m, n,
c                c and a large cx
c       routine called:
c                sphy for computing the spherical bessel
c                functions of the second kind      
c       ========================================================
c
	implicit double precision (a-h,o-z)
	dimension df(200),sy(0:251),dy(0:251)
	eps=1.0d-14
	ip=1
	nm1=int((n-m)/2)
	if (n-m.eq.2*nm1) ip=0
	nm=25+nm1+int(c)
	reg=1.0d0
	if (m+nm.gt.80) reg=1.0d-200
	nm2=2*nm+m
	cx=c*x
	call sphy(nm2,cx,nm2,sy,dy)
	r0=reg
	do 10 j=1,2*m+ip
10         r0=r0*j
	r=r0    
	suc=r*df(1)
	do 15 k=2,nm
	   r=r*(m+k-1.0)*(m+k+ip-1.5d0)/(k-1.0d0)/(k+ip-1.5d0)
	   suc=suc+r*df(k)
	   if (k.gt.nm1.and.dabs(suc-sw).lt.dabs(suc)*eps) go to 20
15         sw=suc
20      a0=(1.0d0-kd/(x*x))**(0.5d0*m)/suc
	r2f=0.0
	do 50 k=1,nm
	   l=2*k+m-n-2+ip
	   if (l.eq.4*int(l/4)) lg=1
	   if (l.ne.4*int(l/4)) lg=-1
	   if (k.eq.1) then
	      r=r0
	   else
	      r=r*(m+k-1.0)*(m+k+ip-1.5d0)/(k-1.0d0)/(k+ip-1.5d0)
	   endif
	   np=m+2*k-2+ip
	   r2f=r2f+lg*r*(df(k)*sy(np))
	   eps1=dabs(r2f-sw)
	   if (k.gt.nm1.and.eps1.lt.dabs(r2f)*eps) go to 55
50         sw=r2f
55      id1=int(log10(eps1/dabs(r2f)+eps))
	r2f=r2f*a0
	if (np.ge.nm2) then
	   id=10
	   return
	endif
	b0=kd*m/x**3.0d0/(1.0-kd/(x*x))*r2f                
	sud=0.0d0
	do 60 k=1,nm
	   l=2*k+m-n-2+ip
	   if (l.eq.4*int(l/4)) lg=1
	   if (l.ne.4*int(l/4)) lg=-1
	   if (k.eq.1) then
	      r=r0
	   else
	      r=r*(m+k-1.0)*(m+k+ip-1.5d0)/(k-1.0d0)/(k+ip-1.5d0)
	   endif
	   np=m+2*k-2+ip
	   sud=sud+lg*r*(df(k)*dy(np))
	   eps2=dabs(sud-sw)
	   if (k.gt.nm1.and.eps2.lt.dabs(sud)*eps) go to 65
60         sw=sud
65      r2d=b0+a0*c*sud       
	id2=int(log10(eps2/dabs(sud)+eps))
	id=max(id1,id2)
	return
	end


	subroutine sphy(n,x,nm,sy,dy)
c
c       ======================================================
c       purpose: compute spherical bessel functions yn(x) and
c                their derivatives
c       input :  x --- argument of yn(x) ( x ע 0 )
c                n --- order of yn(x) ( n = 0,1,תתת )
c       output:  sy(n) --- yn(x)
c                dy(n) --- yn'(x)
c                nm --- highest order computed
c       ======================================================
c
	implicit double precision (a-h,o-z)
	dimension sy(0:n),dy(0:n)
	nm=n
	if (x.lt.1.0d-60) then
	   do 10 k=0,n
	      sy(k)=-1.0d+300
10            dy(k)=1.0d+300
	   return
	endif
	sy(0)=-dcos(x)/x
	sy(1)=(sy(0)-dsin(x))/x
	f0=sy(0)
	f1=sy(1)
	do 15 k=2,n
	   f=(2.0d0*k-1.0d0)*f1/x-f0
	   sy(k)=f
	   if (dabs(f).ge.1.0d+300) go to 20              
	   f0=f1
15         f1=f
20      nm=k-1
	dy(0)=(dsin(x)+dcos(x)/x)/x
	do 25 k=1,nm
25         dy(k)=sy(k-1)-(k+1.0d0)*sy(k)/x
	return
	end


	subroutine rmn2so(m,n,c,x,cv,df,kd,r2f,r2d)
c
c       =============================================================
c       purpose: compute oblate radial functions of the second kind
c                with a small argument, rmn(-ic,ix) & rmn'(-ic,ix)
c       routines called:
c            (1) sckb for computing the expansion coefficients c2k
c            (2) kmn for computing the joining factors
c            (3) qstar for computing the factor defined in (15.7.3)
c            (4) cbk for computing the the expansion coefficient
c                defined in (15.7.6)
c            (5) gmn for computing the function defined in (15.7.4)
c            (6) rmn1 for computing the radial function of the first
c                kind
c       =============================================================
c
	implicit double precision (a-h,o-z)
	dimension bk(200),ck(200),df(200),dn(200)
	if (dabs(df(1)).le.1.0d-280) then
	   r2f=1.0d+300
	   r2d=1.0d+300
	   return
	endif
	eps=1.0d-14
	pi=3.141592653589793d0
	nm=25+int((n-m)/2+c)
	ip=1
	if (n-m.eq.2*int((n-m)/2)) ip=0
	call sckb(m,n,c,df,ck)
	call kmn(m,n,c,cv,kd,df,dn,ck1,ck2)
	call qstar(m,n,c,ck,ck1,qs,qt)
	call cbk(m,n,c,cv,qt,ck,bk)
	if (x.eq.0.0d0) then
	   sum=0.0d0
	   do 10 j=1,nm
	      sum=sum+ck(j)
	      if (dabs(sum-sw).lt.dabs(sum)*eps) go to 15
10            sw=sum
15         if (ip.eq.0) then
	      r1f=sum/ck1
	      r2f=-0.5d0*pi*qs*r1f
	      r2d=qs*r1f+bk(1)
	   else if (ip.eq.1) then
	      r1d=sum/ck1
	      r2f=bk(1)
	      r2d=-0.5d0*pi*qs*r1d
	   endif
	   return
	else
	   call gmn(m,n,c,x,bk,gf,gd)
	   call rmn1(m,n,c,x,df,kd,r1f,r1d)
	   h0=datan(x)-0.5d0*pi
	   r2f=qs*r1f*h0+gf
	   r2d=qs*(r1d*h0+r1f/(1.0d0+x*x))+gd
	endif
	return
	end


	subroutine qstar(m,n,c,ck,ck1,qs,qt)
c
c       =========================================================
c       purpose: compute q*mn(-ic) for oblate radial functions
c                with a small argument
c       =========================================================
c
	implicit double precision (a-h,o-z)
	dimension ap(200),ck(200)
	ip=1
	if (n-m.eq.2*int((n-m)/2)) ip=0
	r=1.0d0/ck(1)**2
	ap(1)=r
	do 20 i=1,m
	   s=0.0d0
	   do 15 l=1,i
	      sk=0.0d0
	      do 10 k=0,l
10               sk=sk+ck(k+1)*ck(l-k+1)
15            s=s+sk*ap(i-l+1)
20      ap(i+1)=-r*s
	qs0=ap(m+1)     
	do 30 l=1,m
	   r=1.0d0
	   do 25 k=1,l
25            r=r*(2.0d0*k+ip)*(2.0d0*k-1.0d0+ip)/(2.0d0*k)**2
30         qs0=qs0+ap(m-l+1)*r
	qs=(-1)**ip*ck1*(ck1*qs0)/c
	qt=-2.0d0/ck1*qs
	return
	end


	subroutine cbk(m,n,c,cv,qt,ck,bk)
c
c       =====================================================
c       purpose: compute coefficient bk's for oblate radial
c                functions with a small argument
c       =====================================================
c
	implicit double precision (a-h,o-z)
	dimension bk(200),ck(200),u(200),v(200),w(200)
	eps=1.0d-14
	ip=1
	if (n-m.eq.2*int((n-m)/2)) ip=0
	nm=25+int(0.5*(n-m)+c)
	u(1)=0.0d0
	n2=nm-2
	do 10 j=2,n2
10         u(j)=c*c
	do 15 j=1,n2
15         v(j)=(2.0*j-1.0-ip)*(2.0*(j-m)-ip)+m*(m-1.0)-cv
	do 20 j=1,nm-1
20         w(j)=(2.0*j-ip)*(2.0*j+1.0-ip)
	if (ip.eq.0) then
	   do 40 k=0,n2-1
	      s1=0.0d0
	      i1=k-m+1
	      do 30 i=i1,nm
		 if (i.lt.0) go to 30
		 r1=1.0d0
		 do 25 j=1,k
25                  r1=r1*(i+m-j)/j
		 s1=s1+ck(i+1)*(2.0*i+m)*r1
		 if (dabs(s1-sw).lt.dabs(s1)*eps) go to 35
		 sw=s1
30            continue
35            bk(k+1)=qt*s1
40         continue
	else if (ip.eq.1) then
	   do 60 k=0,n2-1
	      s1=0.0d0
	      i1=k-m+1
	      do 50 i=i1,nm
		 if (i.lt.0) go to 50
		 r1=1.0d0
		 do 45 j=1,k
45                  r1=r1*(i+m-j)/j
		 if (i.gt.0) s1=s1+ck(i)*(2.0*i+m-1)*r1
		 s1=s1-ck(i+1)*(2.0*i+m)*r1
		 if (dabs(s1-sw).lt.dabs(s1)*eps) go to 55
		 sw=s1
50            continue
55            bk(k+1)=qt*s1
60         continue
	endif
	w(1)=w(1)/v(1)
	bk(1)=bk(1)/v(1)
	do 65 k=2,n2
	   t=v(k)-w(k-1)*u(k)
	   w(k)=w(k)/t
65         bk(k)=(bk(k)-bk(k-1)*u(k))/t
	do 70 k=n2-1,1,-1
70         bk(k)=bk(k)-w(k)*bk(k+1)
	return
	end


	subroutine gmn(m,n,c,x,bk,gf,gd)
c
c       ===========================================================
c       purpose: compute gmn(-ic,ix) and its derivative for oblate
c                radial functions with a small argument
c       ===========================================================
c
	implicit double precision (a-h,o-z)
	dimension bk(200)
	eps=1.0d-14
	ip=1
	if (n-m.eq.2*int((n-m)/2)) ip=0
	nm=25+int(0.5*(n-m)+c)
	xm=(1.0d0+x*x)**(-0.5d0*m)
	gf0=0.0d0
	do 10 k=1,nm
	   gf0=gf0+bk(k)*x**(2.0*k-2.0)
	   if (dabs((gf0-gw)/gf0).lt.eps.and.k.ge.10) go to 15
10         gw=gf0
15      gf=xm*gf0*x**(1-ip)
	gd1=-m*x/(1.0d0+x*x)*gf
	gd0=0.0d0
	do 20 k=1,nm
	   if (ip.eq.0) then
	      gd0=gd0+(2.0d0*k-1.0)*bk(k)*x**(2.0*k-2.0)
	   else
	      gd0=gd0+2.0d0*k*bk(k+1)*x**(2.0*k-1.0)
	   endif
	   if (dabs((gd0-gw)/gd0).lt.eps.and.k.ge.10) go to 25
20         gw=gd0
25      gd=gd1+xm*gd0
	return
	end


	subroutine kmn(m,n,c,cv,kd,df,dn,ck1,ck2)
c
c       ===================================================
c       purpose: compute the expansion coefficients of the
c                prolate and oblate spheroidal functions
c                and joining factors
c       ===================================================
c
	implicit double precision (a-h,o-z)
	dimension u(200),v(200),w(200),df(200),dn(200),
     &            tp(200),rk(200)
	nm=25+int(0.5*(n-m)+c)
	nn=nm+m
	cs=c*c*kd
	ip=1
	if (n-m.eq.2*int((n-m)/2)) ip=0
	do 10 i=1,nn+3     
	   if (ip.eq.0) k=-2*(i-1)
	   if (ip.eq.1) k=-(2*i-3)
	   gk0=2.0d0*m+k
	   gk1=(m+k)*(m+k+1.0d0)
	   gk2=2.0d0*(m+k)-1.0d0
	   gk3=2.0d0*(m+k)+3.0d0
	   u(i)=gk0*(gk0-1.0d0)*cs/(gk2*(gk2+2.0d0))
	   v(i)=gk1-cv+(2.0d0*(gk1-m*m)-1.0d0)*cs/(gk2*gk3)
10         w(i)=(k+1.0d0)*(k+2.0d0)*cs/((gk2+2.0d0)*gk3)
	do 20 k=1,m
	   t=v(m+1)
	   do 15 l=0,m-k-1
15            t=v(m-l)-w(m-l+1)*u(m-l)/t
20         rk(k)=-u(k)/t
	r=1.0d0
	do 25 k=1,m
	   r=r*rk(k)
25         dn(k)=df(1)*r
	tp(nn)=v(nn+1)
	do 30 k=nn-1,m+1,-1
	   tp(k)=v(k+1)-w(k+2)*u(k+1)/tp(k+1)
	   if (k.gt.m+1) rk(k)=-u(k)/tp(k)
30      continue
	if (m.eq.0) dnp=df(1)
	if (m.ne.0) dnp=dn(m)
	dn(m+1)=(-1)**ip*dnp*cs/((2.0*m-1.0)*(2.0*m+1.0-4.0*ip)
     &          *tp(m+1))
	do 35 k=m+2,nn
35         dn(k)=rk(k)*dn(k-1)
	r1=1.0d0
	do 40 j=1,(n+m+ip)/2
40         r1=r1*(j+0.5d0*(n+m+ip))
	nm1=(n-m)/2
	r=1.0d0
	do 45 j=1,2*m+ip
45         r=r*j
	su0=r*df(1)
	do 50 k=2,nm
	   r=r*(m+k-1.0)*(m+k+ip-1.5d0)/(k-1.0d0)/(k+ip-1.5d0)
	   su0=su0+r*df(k)
	   if (k.gt.nm1.and.dabs((su0-sw)/su0).lt.1.0d-14) go to 55
50         sw=su0
55      if (kd.eq.1) goto 70
	r2=1.0d0
	do 60 j=1,m
60         r2=2.0d0*c*r2*j
	r3=1.0d0
	do 65 j=1,(n-m-ip)/2
65         r3=r3*j
	sa0=(2.0*(m+ip)+1.0)*r1/(2.0**n*c**ip*r2*r3*df(1))
	ck1=sa0*su0
	if (kd.eq.-1) return
70      r4=1.0d0
	do 75 j=1,(n-m-ip)/2
75         r4=4.0d0*r4*j
	r5=1.0d0
	do 80 j=1,m
80         r5=r5*(j+m)/c
	g0=dn(m)
	if (m.eq.0) g0=df(1)
	sb0=(ip+1.0)*c**(ip+1)/(2.0*ip*(m-2.0)+1.0)/(2.0*m-1.0)
	ck2=(-1)**ip*sb0*r4*r5*g0/r1*su0
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
	      dk1=m+k+1.0d0
	      dk2=2.0d0*(m+k)
	      d2k=2.0d0*m+k
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

