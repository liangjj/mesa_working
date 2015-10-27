*     ------------------------------------------------------------------
*               p o t l
*     ------------------------------------------------------------------
*
*       computes and stores the potential function
*                              2(k-1)
*              yr = sum  a    y      (j,j;r)
*                   j,k   ijk
*
      SUBROUTINE potl(i,rel,yrm)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER(nt=100,kx=8,ktx=3*kx-2,nwd=10)
      DOUBLE PRECISION yrm(nt,kx)
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
      DOUBLE PRECISION yr(nt,kx), yv(nt)
*---------------------------------------------------------------------
        print *, 'Entering Potl for orbital ',i
        nint = ns-ks+1
*       .. clear yv
        do 10 ii = 1,ns
	  yv(ii) = 0.d0
   10   continue
*       .. clear yr
	do 15 m = 1,ks
	  do 15 ii = 1,ns
	    yr(ii,m) = 0.d0
   15   continue
*
*       .. compute contribution from each k
	do 20  k = 0,2*lmax,2
	  ykn   = 0.0d0
	  kfirst = .true.
	  do 30 j = 1,nwf
	    if (k .le. 2*min(l(i),l(j))) then
              cc = 2.d0*a(i,j,k)
	      print *, 'i,j,k,a(i,j,k)',i,j,k,cc
	      if (abs(cc) .gt. 1.d-10) then
		if (k .eq. 0) ykn = ykn + cc
		if ( kfirst ) then
		   kfirst = .false.
*                  .. clear array for fyk
		   do 40 m = 1,ks
		     do 40 jj = 1,nint
		       fyk(jj,m) = 0.d0
 40		   continue
*                  .. clear the array yk
	           do 45 ii=1,ns
	             yk(ii)=0.d0
 45 	           continue
*
*                  set up and factor the differential operator
*
		   print *, ' Potl: calling facdyk'
                   call facdyk(nt,kx,ktx,ks,k,nint,t,bsp,bspd,
     :                         db2,rm2,ipvtd,dyk)
		end if
*   
	        if (ic1 .ne. j) then
	          call yval(0,0,nt,kx,ks,nint,gr,bsp,p(1,j),fc1)
	          ic1 = j
	        end if
		do 50 m = 1,ks
		  do 50 ii = 1,nint
		    fyk(ii,m) = fyk(ii,m) + cc*fc1(ii,m)*fc1(ii,m)
 50             continue
	      end if
	    end if
 30       continue
*
*         add a relativistic correction if rel=TRUE
*
	  if (rel) then
	    c = (2*k+1)*fine
            do 60 m = 1,ks
	      do 60 ii = 1,nint
	        yr(ii,m)= yr(ii,m)+ c*fyk(ii,m)
  60       continue
	  end if
*
* set up the array yk (spline approx. of SUM ()(1/r)fc1(r)fc2(r)
*
	  do 70 m=1,ks
	    do 70 ith=1,ks
	      do 70 ii = 1,nint
   	        irow=ii+ith-1
	        yk(irow)=yk(irow)+grw(ii,m)*grm(ii,m)*bsp(ii,m,ith)*
     :  	      fyk(ii,m)
 70	  continue
*
  	  const=-(2*k+1.d0)
	  do 80 ii=1,ns
	    yk(ii)=const*yk(ii)
 80	  continue
	  yk(1)=0.d0
	  yk(ns) = ykn
          print *, ' Potl: yk'
	  print '(6f12.8)', (yk(ii),ii=1,ns)
*
* solve the matrix equation
*
	  call dgbsl(dyk,ktx,ns,ks-1,ks-1,ipvtd,yk,0)
       print *, 'Potl: after dgbsl', i
       print '(6f12.8)', (yk(ii),ii=1,ns)
*
*          .. add to yv
	   do 90 ii = 1,ns
	     yv(ii) = yv(ii) + yk(ii)
 90        continue
 20    continue 
       print *, 'Potl: yv', i
       print '(6f12.8)', (yv(ii),ii=1,ns)
*
*  evaluates the function Y(r) at all the gaussian points
*
        call yval(1,0,nt,kx,ks,nint,grm,bsp,yv,yr)
	print *, ' Y(r)'
*	do 99 m = 1,ks
	  m = ks/2
	  print '(6f12.8)', (yr(ii,m),ii=1,nint)
99      continue
	do 100 m = 1,ks
	  do 100 ii = 1,nint
 	    yr(ii,m) = yr(ii,m)*grw(ii,m)*grm(ii,m)
 100    continue
	call minty(0,nt,kx,ks,nint,bsp,yr,yrm)
        print *, ' Test yrm'
	call bxv(nt,ks,ns,yrm,p(1,i),yv)
        ans1 = ddot(ns,yv,1,p(1,i),1)
	print *, '<p1,potl,p1>', ans1, quadr(1,1,-1)
*
	end
