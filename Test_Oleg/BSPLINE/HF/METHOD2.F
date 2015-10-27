*     ------------------------------------------------------------------
*               m e t h o d
*     ------------------------------------------------------------------
*
*       uses m1, m2, or m3 to solve the radial equation. if the input
*   data indicated meth(i) = 3, then this  solution  is  returned  to
*   de.  otherwise,  the routine searches for an acceptable  solution
*   which  is  both  positive  near  the  origin and has the required
*   number  of nodes.
*
      SUBROUTINE method(i,first,ed,pde,ms)
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      PARAMETER(nt=100,kx=8,nwd=10)
*
      COMMON /coeff/coef(100),ijptr(5,5)
      character config*50,el*3,atom*6,term*6
      COMMON /label/config,el(nwd),atom,term
*
      COMMON /ORBITS/p(nt,nwd), n(nwd), l(nwd), max(nwd)
      COMMON /ORBSCF/ sum(nwd), s(nwd), dpm(nwd), acc(nwd), meth(nwd),
     :             iord(nwd), e(nwd,nwd)
*
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
*
      COMMON /SPLINE/gr(nt,kx),grm(nt,kx),grw(nt,kx),
     :        bsp(nt,kx,kx),bspd(nt,kx,kx,2),db1(nt,kx),db2(nt,kx),
     :        sb(nt,kx),r1(nt,kx),rm1(nt,kx),rm2(nt,kx),
     :        h0(nt,kx), rkm(nt,nt,kx,kx), lh, krk
*
      LOGICAL fail,omit,rel,all,trace,print,strong
      COMMON /TEST/fail,omit,rel,all,trace
      DOUBLE PRECISION pde(nt),h1(nt,kx),bpv(nt),w(nt),yrm(nt,kx),
     :          xv(nt),a(nt,nt),b(nt,nt),eigval(nt),w1(nt),w2(nt)
      INTEGER ipvtd(nt)
      LOGICAL first
*
*  *****  'first' must be 'true' the first time solve is called for
*  *****  potential and exchange to be computed
*  *****  'eu' is the upper bound of the energy parameter
*  *****  'em' is the minimum value of the energy parameter
*
      fail = .false.
      mit = 0
      if (first) then
        first = .false.
        em = 0.d0
        eu = ((z - dmin1(0.5d0*s(i),d2*s(i)))/n(i))**2
        fu = eu
	node = n(i)-l(i) -1
 	call potl(i,rel,yrm)
	call xch(i,rel,xv)
*	print *, 'xv after xch'
*	print '(6f12.8)',(xv(ii),ii=1,ns)
	call hlm(z,l(i),nt,kx,ks,ns,db2,rm1,rm2,h1)
	do 10 m = 1,ks
	  do 10 ii = 1,ns
	    h1(ii,m) = h1(ii,m) - yrm(ii,m)
10      continue
	call bxv(nt,ks,ns,h1,p(1,i),bpv)
	do 12 ii = 1,ns
	  bpv(ii) = bpv(ii) - xv(ii)
12      continue
	ed = ddot(ns,p(1,i),1,bpv,1)
	print *, 'ed =',ed
	do 14 j = 1,nwf
	  eij = e(i,j)
	  print *, 'eij',eij
	  if( eij .ne. 0.d0 .and. eij .ne. 1.d-10 .and. i.ne.j) then
	    call bxvpw(nt,ks,ns,eij,sb,p(1,j),xv,w)
	  end if
14	continue
	print *, 'xv'
	print '(6f12.8)',(xv(ii),ii=1,ns)
      end if
17    call solve(i,ed,h1,xv,pde,ms)
      print *, 'After solve', ms
      print '(6f12.8)', (pde(ii),ii=1,ns)
*
*  *****  if mk equals 1,2 check the nodes
* 
      if (mk .le. 2) then
*
*  *****  count the number of nodes
*
        mn = ms
        nc = nodec(pde,ns,mn)
        if (trace) write (6,99) el(i),nc,mn,pde(mn),ed,eu,em
99      format(2x,a3,' nc =',i3,' mn =',i3,' pde(mn) =',
     :         d10.2,' ed =',d10.2,' eu =',d10.2,' em =',d10.2)
*
*  *****  if node count is off or pn too small set up matrix
*
        if (iabs(nc-node).eq.0 .and. abs(pde(mn)).gt. .05d0/sqrt(z))then
*
*         .. we have the correct number of nodes
          if (pde(mn) .lt. 0.d0 ) then
            do 52 j = 1,ns
              pde(j) = - pde(j)
52          continue
	  end if
	  e(i,i) = ed
          return
	else  
*
	  if (mit .eq. 0) then
	    mit = mit+1
	    call coulom(nt,kx,ks,ns,sb,h1,a,b)
*           .. apply boundary condition
            do 500 ii = 2,ks
        	    a(1,ii) = 0.d0
        	    a(ii,1) = 0.d0
        	    a(1,1) = -1.d10
        	    b(1,ii) = 0.d0
        	    b(ii,1) = 0.d0
        	    b(1,1) = 1.0d0
500         continue
	    
	    call rsg(nt,ns-1,a,b,eigval,0,eigvec,w1,w2,ier)
	    if (ier .ne. 0) then
	      print *, 'Error in Eigenvalues routine'
	      stop
	    else
	      print *, ' Method called RSG: eigenvalues'
	      print '(4(I3,1PD12.4))', (ii, eigval(ii),ii=ns-1,1,-1)
	    end if
	    ed = eigval(ns-node-1)
	    print *, ' RSG eig = ',ed
*           .. regenerate the matrix
            call coulom(nt,kx,ks,ns,sb,h1,a,b)
*           .. from  (A - e B)
	    do 502 jj = 1,ns
	      do 502 ii = 1,ns
		a(ii,jj) = a(ii,jj) - ed*b(ii,jj)
502         continue
*           .. apply boundary condition
            do 510 ii = 2,ks
        	    a(1,ii) = 0.d0
        	    a(ii,1) = 0.d0
        	    a(1,1) = -1.d10
        	    b(1,ii) = 0.d0
        	    b(ii,1) = 0.d0
        	    b(1,1) = 1.0d0
510         continue
 	    call ginvit(a,nt,ns-1,pde,ipvtd,w1,ier)
	    if (ier .ne. 0) stop
	    pde(ns) = 0.d0
	    print *, 'Pde after ginvit'
	    print '(6f12.8)', (pde(ii),ii=1,ns)
	    pn = 1.d0/sqrt(quad(ns,ks,pde,pde))
	    if (pde(4) .lt. 0.d0) pn = -pn
	    do 520 ii = 1,ns
	      pde(ii) = pn*pde(ii)
520         continue
	    print *, ' Normalized  after ginvit'
	    print '(6f12.8)', (pde(ii),ii=1,ns)
	    ed = ed - ddot(ns,p(1,i),1,xv,1)
	    print *, 'New ed',ed
            go to 17
	  end if
	  fail = .true.
	end if
      else
	e(i,i) = ed
      end if
      end
