*     ------------------------------------------------------------------
*               x c h
*     ------------------------------------------------------------------
*
*       this routine computes functions associated with the exchange
*   function for the i'th radial equation,  including   contributions
*   from   the  interactions. 
*
      SUBROUTINE xch(i,rel,xv)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER (nt=100,kx=8,ktx=3*kx-2,nwd=10)
      DOUBLE PRECISION xv(nt)
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
     :        h0(nt,kx), rkm(nt,nt,kx,kx),lh,krk

      COMMON /SLATER/fc1(nt,kx),fc2(nt,kx),yk(nt),fyk(nt,kx),
     :        dyk(ktx,nt),ipvtd(nt),fc(nt,kx),ic1,ic2,iy1,iy2,ky
*
      LOGICAL kfirst, REL
      DOUBLE PRECISION xr(nt,kx)
*---------------------------------------------------------------------
*       
        nint = ns-ks+1
*       .. clear xr
	do 15 m = 1,ks
	  do 15 ii = 1,ns
	    xr(ii,m) = 0.d0
   15   continue
        if (ic1 .ne. i) then
           call yval(0,0,nt,kx,ks,nint,gr,bsp,p(1,i),fc1)
	   ic1 = i
        end if
*
*       .. compute contribution from each k
	do 20  k = 0,2*lmax,2
	  ykn   = 0.0d0
	  kfirst = .true.
	  do 30 j = 1,nwf
	    if (j .ne. i) then
	      if (iabs(l(i)-l(j)) .le. k .and. k.le. l(i)+l(j)) then
		if (mod(k,2) .eq. mod(l(i)+l(j),2)) then
                  cc = 2.d0*b(i,j,k)
		  print *, 'i,j,k,b(i,j,k)',i,j,k,cc
	          if (abs(cc) .gt. 1.d-10) then
		    if ( kfirst ) then
		      kfirst = .false.
*                     .. clear the array yk
	              do 45 ii=1,ns
	                yk(ii)=0.d0
 45 	              continue
*
*                     set up and factor the differential operator
*
                      call facdyk(nt,kx,ktx,ks,k,nint,t,bsp,bspd,
     :                            db2,rm2,ipvtd,dyk)
*   
		    end if
	            if (ic2 .ne. j) then
	              call yval(0,0,nt,kx,ks,nint,gr,bsp,p(1,j),fc2)
	              ic2 = j
	            end if
		    do 50 m = 1,ks
		      do 50 ii = 1,nint
		        fyk(ii,m) = cc*fc1(ii,m)*fc2(ii,m)
 50                 continue
	          end if
*
*                 add a relativistic correction if rel=TRUE
*
	          if (rel) then
	            c = (2*k+1)*fine
                    do 60 m = 1,ks
	              do 60 ii = 1,nint
	                xr(ii,m)= xr(ii,m)+ c*fyk(ii,m)
  60                continue
	          end if
*
* set up the array xk (spline approx. of SUM ()(1/r)fc1(r)fc2(r)
*
  	          do 70 m=1,ks
	            do 70 ith=1,ks
	              do 70 ii = 1,nint
         	        irow=ii+ith-1
	                yk(irow)=yk(irow)+
     :                      grw(ii,m)*grm(ii,m)*bsp(ii,m,ith)*fyk(ii,m)
 70	          continue
*
  	          const=-(2*k+1.d0)
	          do 80 ii=1,ns
	            yk(ii)=const*yk(ii)
 80	          continue
	          yk(1)=0.d0
	          yk(ns) = ykn
*
* solve the matrix equation
*
		  print *, 'yk before dgbsl'
		  print '(6f12.8)', (yk(ii),ii=1,ns)
	          call dgbsl(dyk,ktx,ns,ks-1,ks-1,ipvtd,yk,0)
                  call yval(0,0,nt,kx,ks,nint,grm,bsp,yk,fyk)
		  do 90 m = 1,ks
		    do 90 ii = 1,nint
		      xr(ii,m) = xr(ii,m) + fyk(ii,m)*fc2(ii,m)
 90               continue
		end if
	      end if
	    end if
 30       continue
*
 20    continue
*
*  evaluates the function x(r) at all the gaussian points
*
	do 100 m = 1,ks
	  do 100 ii = 1,nint
	    xr(ii,m) = xr(ii,m)*grw(ii,m)*grm(ii,m)
 100    continue
	call vinty(nt,kx,ks,nint,bsp,xr,xv)
*
	end
