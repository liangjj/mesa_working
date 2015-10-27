*
*     ------------------------------------------------------------------
*               s o l v e
*     ------------------------------------------------------------------
*
*       when first is .true., solve computes the potential and exchange
*   function and initializes variables for the i'th radial  equation.
*   the vector p1 is the solution of the radial equation and p2 the
*   variation of the solution with respect to the energy parameter
*   e(i,i).
*
*
      SUBROUTINE solve(iel,first,ed,h1,xv,pde,ms)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER(nt=100,kx=8,ktx=3*kx-2,nwd=10)
      character config*50,el*3,atom*6,term*6
      COMMON /label/config,el(nwd),atom,term
      COMMON /ORBITS/p(nt,nwd), n(nwd), l(nwd), max(nwd)
      COMMON /ORBSCF/ sum(nwd), s(nwd), dpm(nwd), acc(nwd), meth(nwd),
     :             iord(nwd), e(nwd,nwd)
*
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
      LOGICAL fail,omit,rel,all,trace,print,strong
      COMMON /TEST/fail,omit,rel,all,trace
      LOGICAL first
*
      COMMON /SPLINE/gr(nt,kx),grm(nt,kx),grw(nt,kx),
     :        bsp(nt,kx,kx),bspd(nt,kx,kx,2),db1(nt,kx),db2(nt,kx),
     :        sb(nt,kx),r1(nt,kx),rm1(nt,kx),rm2(nt,kx),
     :        h0(nt,kx), rkm(nt,nt,kx,kx), lh, krk
*
       DOUBLE PRECISION bandde(ktx,nt),h1(nt,kx),xv(nt),w(nt),pde(nt)
       INTEGER ipvtd(nt)
*--------------------------------------------------------------------
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
	    bandde(3*ks-1-j,i-ks+j)=h1(i,j) -ed*sb(i,j)
 30	continue            
*      ... upper portion
	do 40 j=2,ks
	  do 40 i=1,ns-j+1
	    bandde(2*ks-j,i+j-1)= h1(i+j-1,ks-j+1)-ed*sb(i+j-1,ks-j+1)
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
*       
*       copy xv to pde
	do 60 i = 1,ns
	  pde(i) = xv(i)
60      continue
	pde(1) = 0.d0
	pde(ns) = 0.d0
	call dgbsl(bandDE,ktx,ns,ks-1,ks-1,ipvtd,pde,0)
* 
	if (meth(i) .eq. 2) then
	  do 70 i = 1,ns
	    w(i) = p(i,iel)
70        continue
	  call dgbsl(bandDE,ktx,ns,ks-1,ks-1,ipvtd,w,0)
	  a00 = quad(ns,ks,pde,pde)
	  a01 = quad(ns,ks,pde,w)
	  a11 = quad(ns,ks,w,w)
	  print *, 'a00,a01,a11',a00,a01,a11
	  b = 2.d0*a01/a11
          c = (a00 - 1.d0)/a11
	  if (b .gt. 0.d0) then
	    den = b + sqrt(b*b-4*c)
	  else
	    den = b - sqrt(b*b-4*c)
          end if
	  deltae = - 2.d0*c/den
	  print *, 'b,c,deltae',b,c,deltae
	  do 80 i = 1,ns
	    pde(i) = pde(i) + deltae*w(i)
80        continue
        end if
	ms = ns
	end	
