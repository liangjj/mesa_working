*=====================================================================
	subroutine ykf(j1,j2,k,rel)
*=====================================================================
*
*     Computes the spline solution of the differential equation
*           (d^2-k(k+1)/r^2)yk(r)= -(2k+1)(1/r) P_a(r)P_b(r)
*
*     SUBROUTINES called:
*       sdyk
*       yvalc
*       mintm
*
*     on exit
*     -------
*     yk    spline expansion of the solution
*     fyk   values at the gaussian points * (1/r)
*
*---------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER (nt=100,kx=8,ktx=3*kx-2,nwd=10)
*
      COMMON /ORBITS/p(nt,nwd), n(nwd), l(nwd), max(nwd)
*
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
*
      COMMON /SPLINE/gr(nt,kx),grm(nt,kx),grw(nt,kx),
     :        bsp(nt,kx,kx),bspd(nt,kx,kx,2),db1(nt,kx),db2(nt,kx),
     :        sb(nt,kx),r1(nt,kx),rm1(nt,kx),rm2(nt,kx),
     :        h0(nt,kx), rkm(nt,nt,kx,kx),lh,krk

      COMMON /SLATER/fc1(nt,kx),fc2(nt,kx),yk(nt),fyk(nt,kx),
     :        dyk(ktx,nt),ipvtd(nt),fc(nt,kx),ic1,ic2,iy1,iy2,ky
*
      LOGICAL same, REL
*---------------------------------------------------------------------
*       
        nint = ns-ks+1
        if ( k.eq.0 .and. j1 .eq. j2) then
  	  ykn   = 1.d0
        else
	  ykn   = 0.0d0
        end if
        if (j1 .eq. j2) then
	   same = .true.
        else
	   same = .false.
        end if
*
	nint=ns-ks+1
*
*  set up and factor the differential operator
*
	print *, 'YKF calling facdyk'
        call facdyk(nt,kx,ktx,ks,k,nint,t,bsp,bspd,db2,rm2,ipvtd,dyk)
*
* clear the array yk
*
	do 10 i=1,ns
	  yk(i)=0.d0
 10 	continue
*
*       .. create the numerical values of the orbitals on the grid
*   
	if (.not. same) then
	  if (ic1 .ne. j1) then
	    call yval(0,0,nt,kx,ks,nint,gr,bsp,p(1,j1),fc1)
	    ic1 = j1
	  end if
	  if (ic2 .ne. j2) then
	    call yval(0,0,nt,kx,ks,nint,gr,bsp,p(1,j2),fc2)
	    ic2 = j2
          end if
*
* set up the array yk (spline approx. of (1/r)fc1(r)fc2(r)
*
	  do 100 m=1,ks
	    do 100 ith=1,ks
	      do 100 ii = 1,nint
   	      i=ii+ith-1
	        yk(i)=yk(i)+grw(ii,m)*grm(ii,m)*bsp(ii,m,ith)*
     :  	     fc1(ii,m)*fc2(ii,m)
 100	  continue
	else
	  if (ic1 .ne. j1) then
	    call yval(0,0,nt,kx,ks,nint,gr,bsp,p(1,j1),fc1)
	    ic1 = j1
	  end if
*
* set up the array yk (spline approx. of (1/r)fc1(r)fc2(r)
*
	  do 110 m=1,ks
	    do 110 ith=1,ks
	      do 110 ii = 1,nint
   	      i=ii+ith-1
	        yk(i)=yk(i)+grw(ii,m)*grm(ii,m)*bsp(ii,m,ith)*
     :  	     fc1(ii,m)*fc1(ii,m)
 110	continue
	end if
*
	const=-(2*k+1.d0)
	do 200 i=1,ns
	  yk(i)=const*yk(i)
 200	continue
	yk(1)=0.d0
	yk(ns) = ykn
*
* solve the matrix equation
*
	call dgbsl(dyk,ktx,ns,ks-1,ks-1,ipvtd,yk,0)
*
*  evaluates the function yk{isp,ith} at all the gaussian points
*
        call yval(0,1,nt,kx,ks,nint,grm,bsp,yk,fyk)
*
*	print *, 'Values of r and yk at first gaussian point'
*	do 210 i = 1,nint
*	rr = 1/grm(i,1)
*	print *, rr, fyk(i,1)*rr, exact(rr)
*  210   continue
*	stop

*
*  add a relativistic correction if rel=TRUE
*
	if (rel) then
	  c = (2*k+1)*fine
	  if (.not. same) then
            do 300 m = 1,ks
	       do 300 i = 1,nint
	          fyk(i,m)=fyk(i,m)+ c*grm(i,m)*fc1(i,m)*fc2(i,m)
  300       continue
	  else
            do 310 m = 1,ks
	       do 310 i = 1,nint
	          fyk(i,m)=fyk(i,m)+ c*grm(i,m)*fc1(i,m)*fc1(i,m)
  310       continue
	  end if
        end if
	end
