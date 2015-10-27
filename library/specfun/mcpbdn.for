	program mcpbdn
c
c       =============================================================
c       purpose: this program computes parabolic cylinder functions 
c                dn(z) for an integer order and a complex argument
c                using subroutine cpbdn
c       input :  x --- real part of z
c                y --- imaginary part of z
c                n --- order of dn(z)
c       output:  cpb(|n|) --- dn(z)
c                cpd(|n|) --- dn'(z)
c       example:
c                z = 5.0+ 5.0 i
c
c     n     re[dn(z)]      im[dn(z)]      re[dn'(z)]     im[dn'(z)]
c   -----------------------------------------------------------------
c     0   .99779828d+00  .66321897d-01 -.23286910d+01 -.26603004d+01
c     1   .46573819d+01  .53206009d+01  .26558457d+01 -.24878635d+02
c     2  -.43138931d+01  .49823592d+02  .14465848d+03 -.10313305d+03
c     3  -.28000219d+03  .21690729d+03  .12293320d+04  .30720802d+03
c     4  -.24716057d+04 -.46494526d+03  .38966424d+04  .82090067d+04
c    -1   .10813809d+00 -.90921592d-01 -.50014908d+00 -.23280660d-01
c    -2   .24998820d-02 -.19760577d-01 -.52486940d-01  .47769856d-01
c    -3  -.15821033d-02 -.23090595d-02 -.68249161d-03  .10032670d-01
c    -4  -.37829961d-03 -.10158757d-03  .89032322d-03  .11093416d-02
c       =============================================================
c
	implicit double precision (x,y)
	implicit complex*16 (c,z)
	dimension cpb(0:100),cpd(0:100)
	write(*,*)'please enter n, x and y '
	read(*,*)n,x,y
	write(*,20)n,x,y
	z=cmplx(x,y)
	n0=abs(n)
	call cpbdn(n,z,cpb,cpd)
	write(*,*)
	if (n.ge.0) then
	   write(*,*)'  n     re[dn(z)]       im[dn(z)]       ',
     &               're[dn''(z)]      im[dn''(z)]'
	else
	   write(*,*)' -n     re[dn(z)]       im[dn(z)]       ',
     &               're[dn''(z)]      im[dn''(z)]'
	endif
	write(*,*)'-------------------------------------------',
     &            '-------------------------'
	do 5 i=0,n0
5          write(*,10)i,cpb(i),cpd(i)
10      format(1x,i3,4d16.8)
20      format(1x,'n =',i3,',   z =x+iy :',f6.2,'+',f6.2,' i')
	end


	subroutine cpbdn(n,z,cpb,cpd)
c
c       ==================================================
c       purpose: compute the parabolic cylinder functions 
c                 dn(z) and dn'(z) for a complex argument
c       input:   z --- complex argument of dn(z)
c                n --- order of dn(z)  ( n=0,ס1,ס2,תתת )
c       output:  cpb(|n|) --- dn(z)
c                cpd(|n|) --- dn'(z)
c       routines called:
c            (1) cpdsa for computing dn(z) for a small |z|
c            (2) cpdla for computing dn(z) for a large |z|
c       ==================================================
c
	implicit double precision (a-b,d-h,o-y)
	implicit complex*16 (c,z)
	dimension cpb(0:*),cpd(0:*)
	pi=3.141592653589793d0
	x=real(z)
	a0=cdabs(z)
	c0=(0.0d0,0.0d0)
	ca0=cdexp(-0.25d0*z*z)
	if (n.ge.0) then
	   cf0=ca0
	   cf1=z*ca0
	   cpb(0)=cf0
	   cpb(1)=cf1
	   do 10 k=2,n
	      cf=z*cf1-(k-1.0d0)*cf0
	      cpb(k)=cf
	      cf0=cf1
10            cf1=cf
	else
	   n0=-n
	   if (x.le.0.0.or.cdabs(z).eq.0.0) then
	      cf0=ca0
	      cpb(0)=cf0
	      z1=-z
	      if (a0.le.7.0) then
		 call cpdsa(-1,z1,cf1)
	      else
		 call cpdla(-1,z1,cf1)
	      endif
	      cf1=dsqrt(2.0d0*pi)/ca0-cf1
	      cpb(1)=cf1
	      do 15 k=2,n0
		 cf=(-z*cf1+cf0)/(k-1.0d0)
		 cpb(k)=cf
		 cf0=cf1
15               cf1=cf
	   else
	      if (a0.le.3.0) then
		 call cpdsa(-n0,z,cfa)
		 cpb(n0)=cfa
		 n1=n0+1
		 call cpdsa(-n1,z,cfb)
		 cpb(n1)=cfb
		 nm1=n0-1
		 do 20 k=nm1,0,-1
		    cf=z*cfa+(k+1.0d0)*cfb
		    cpb(k)=cf
		    cfb=cfa
20                  cfa=cf
	      else
		 m=100+abs(n)
		 cfa=c0
		 cfb=(1.0d-30,0.0d0)
		 do 25 k=m,0,-1
		    cf=z*cfb+(k+1.0d0)*cfa
		    if (k.le.n0) cpb(k)=cf
		    cfa=cfb
25                  cfb=cf
		 cs0=ca0/cf
		 do 30 k=0,n0
30                  cpb(k)=cs0*cpb(k)
	      endif
	   endif
	endif
	cpd(0)=-0.5d0*z*cpb(0)
	if (n.ge.0) then
	   do 35 k=1,n
35            cpd(k)=-0.5d0*z*cpb(k)+k*cpb(k-1)
	else
	   do 40 k=1,n0
40            cpd(k)=0.5d0*z*cpb(k)-cpb(k-1)
	endif
	return
	end


	subroutine cpdsa(n,z,cdn)
c
c       ===========================================================
c       purpose: compute complex parabolic cylinder function dn(z)
c                for small argument
c       input:   z   --- complex argument of d(z)
c                n   --- order of d(z) (n = 0,-1,-2,תתת)
c       output:  cdn --- dn(z)
c       routine called: gaih for computing ג(x), x=n/2 (n=1,2,...)
c       ===========================================================
c
	implicit double precision (a-b,d-h,o-y)
	implicit complex*16 (c,z)
	eps=1.0d-15
	pi=3.141592653589793d0
	sq2=dsqrt(2.0d0)
	ca0=cdexp(-.25d0*z*z)
	va0=0.5d0*(1.0d0-n)
	if (n.eq.0.0) then
	   cdn=ca0
	else
	   if (cdabs(z).eq.0.0) then
	      if (va0.le.0.0.and.va0.eq.int(va0)) then
		 cdn=0.0d0
	      else
		 call gaih(va0,ga0)
		 pd=dsqrt(pi)/(2.0d0**(-.5d0*n)*ga0)
		 cdn=cmplx(pd,0.0d0)
	      endif
	   else
	      xn=-n
	      call gaih(xn,g1)
	      cb0=2.0d0**(-0.5d0*n-1.0d0)*ca0/g1
	      vt=-.5d0*n
	      call gaih(vt,g0)
	      cdn=cmplx(g0,0.0d0)
	      cr=(1.0d0,0.0d0)
	      do 10 m=1,250
		 vm=.5d0*(m-n)
		 call gaih(vm,gm)
		 cr=-cr*sq2*z/m
		 cdw=gm*cr
		 cdn=cdn+cdw
		 if (cdabs(cdw).lt.cdabs(cdn)*eps) go to 20
10            continue
20            cdn=cb0*cdn
	   endif
	endif
	return
	end


	subroutine cpdla(n,z,cdn)
c
c       ===========================================================
c       purpose: compute complex parabolic cylinder function dn(z)
c                for large argument
c       input:   z   --- complex argument of dn(z)
c                n   --- order of dn(z) (n = 0,ס1,ס2,תתת)
c       output:  cdn --- dn(z)
c       ===========================================================
c
	implicit double precision (a-b,d-h,o-y)
	implicit complex*16 (c,z)
	cb0=z**n*cdexp(-.25d0*z*z)
	cr=(1.0d0,0.0d0)
	cdn=(1.0d0,0.0d0)
	do 10 k=1,16
	   cr=-0.5d0*cr*(2.0*k-n-1.0)*(2.0*k-n-2.0)/(k*z*z)
	   cdn=cdn+cr
	   if (cdabs(cr).lt.cdabs(cdn)*1.0d-12) go to 15
10      continue
15      cdn=cb0*cdn
	return
	end


	subroutine gaih(x,ga)
c
c       =====================================================
c       purpose: compute gamma function ג(x)
c       input :  x  --- argument of ג(x), x = n/2, n=1,2,תתת
c       output:  ga --- ג(x)
c       =====================================================
c
	implicit double precision (a-h,o-z)
	pi=3.141592653589793d0
	if (x.eq.int(x).and.x.gt.0.0) then
	   ga=1.0d0
	   m1=int(x-1.0)
	   do 10 k=2,m1
10            ga=ga*k
	else if (x+.5d0.eq.int(x+.5d0).and.x.gt.0.0) then
	   m=int(x)
	   ga=dsqrt(pi)
	   do 15 k=1,m
15            ga=0.5d0*ga*(2.0d0*k-1.0d0)
	endif
	return
	end
