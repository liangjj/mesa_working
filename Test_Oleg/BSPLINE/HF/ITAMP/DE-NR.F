*     ------------------------------------------------------------------
*    		D E 
*     ------------------------------------------------------------------
*
*       This routine controls the improvement of the radial functions
*   with the same orbital quantum number.  
*
*   The Newton Raphson iteration is used to solve the
*   non-linear systems of equations.
*
      SUBROUTINE de(iel)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER (nt=100,kx=8,ktx=3*kx-2,nwd=10)
*
      CHARACTER config*50,el*3,atom*6,term*6
      COMMON /LABEL/config,el(nwd),atom,term
*
      COMMON /ORBITS/p(nt,nwd), n(nwd), l(nwd), max(nwd)
      COMMON /ORBSCF/sum(nwd), s(nwd), dpm(nwd), acc(nwd), meth(nwd),
     :               iord(nwd), e(nwd,nwd)
*
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
*
      COMMON /SPLINE/gr(nt,kx),grm(nt,kx),grw(nt,kx),
     :        bsp(nt,kx,kx),bspd(nt,kx,kx,2),db1(nt,kx),db2(nt,kx),
     :        sb(nt,kx),r1(nt,kx),rm1(nt,kx),rm2(nt,kx),
     :        h0(nt,kx), rkm(nt,nt,kx,kx), lh, krk
*
       DOUBLE PRECISION bandde(ktx,nt),h1(nt,kx),bpv(nt),w(nt),
     :        yrm(nt,kx), xv(nt), ipvtd(nt)
*--------------------------------------------------------------------
*
	print *, 'DE:  nwf',nwf, 'iel',iel
 	call potl(iel,rel,yrm)
	call xch(iel,rel,xv)
	call hlm(z,l(iel),nt,kx,ks,ns,db2,rm1,rm2,h1)
	do 10 m = 1,ks
	  do 10 ii = 1,ns
	    h1(ii,m) = h1(ii,m) - yrm(ii,m)
10      continue
	call bxv(nt,ks,ns,h1,p(1,iel),bpv)
	do 12 ii = 1,ns
	  bpv(ii) = bpv(ii) - xv(ii)
12      continue
	eii = ddot(ns,p(1,iel),1,bpv,1)
	print *,'eii', eii
	do 14 m = 1,ks
	  do 14 ii = 1,nt
	     h1(ii,m) = h1(ii,m) - eii*sb(ii,m)
14      continue
	do 16 j = 1,nwf
	  eij = e(iel,j)
	  if( eij .ne. 0 .and. iej .ne. 1.d-10 .and. iel.ne.j) then
	    call bxvpw(nt,ks,nc,eij,sb,p(1,j),xv,w)
	  end if
16	continue
*
*  clear the bandde array
*
  	do 20 i=1,ks-1
	  do 20 j=1,ns
            bandde(i,j)=0.d0
 20	continue 
*
*  set up bandde
*
*       ... lower portion
	do 30 j=1,ks
   	  do 30 i=ks-j+1,ns
	    bandde(3*ks-1-j,i-ks+j)=h1(i,j)
 30	continue            
*      ... upper portion
	do 40 j=2,ks
	  do 40 i=1,ns-j+1
	    bandde(2*ks-j,i+j-1)= h1(i+j-1,ks-j+1)
 40	continue
*
* apply the zero boundary condition at the origin and at rmax
*
	do 50 i=1,ks
*           .. diagonals
      	    j=2*ks-i
	    bandde(j,i)=0.d0
	    j=3*ks-1-i
	    bandde(j,ns-ks+i)=0.d0
*           .. columns
	    bandde(2*ks -2 +i,1) = 0.d0
	    bandde(ks-1+i,ns) = 0.d0
 50	continue
 	bandde(2*ks-1,1)=1.d0 
        bandde(2*ks-1,ns) = 1.d0

*	print *, ' Band DE'
*	do 99 m = 1,3*ks-2
*	  print '(6f12.8)', (bandde(m,ii),ii=1,ns)
*99      continue
*	print *, ' h1'
*	do 199 m = 1,ks
*	  print '(6f12.8)', (h1(ii,m),ii=1,ns)
*199     continue
*
* factorize bandd
*
	call dgbfa(bandde,ktx,ns,ks-1,ks-1,ipvtd,ier)
	if (ier .ne. 0) then
          print *, 'bandDE is non_invertible'
	  print *, ' Ier = ',ier
	  stop
	end if
	call dgbsl(bandDE,ktx,ns,ks-1,ks-1,ipvtd,xv,0)
* 
	return
	end	
