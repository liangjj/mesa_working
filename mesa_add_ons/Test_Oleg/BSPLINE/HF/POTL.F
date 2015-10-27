!=======================================================================
  SUBROUTINE potl 
!=======================================================================
!
!       computes and stores the potential function
!                              2(k-1)
!              yr = sum  a    y      (j,j;r)
!                   j,k   ijk
!----------------------------------------------------------------------
!
    Use Spline_Param
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(nv) :: yv
    REAL(KIND=8), DIMENTION(nv,ks) :: yr 

    print *, 'Entering Potl for orbital ',i
    
    !   .. clear yv
    yv = 0
    !  .. clear yr
    yr = 0

    !  .. compute contribution from each k
    do k = 0,2*lmax,2
      ykn   = 0.0d0
      kfirst = .true.
      do j = 1,nwf
	if (k <= 2*min(l(i),l(j))) then
          cc = 2.d0*a(i,j,k)
	  print *, 'i,j,k,a(i,j,k)',i,j,k,cc
	  if (abs(cc) > 1.d-10) then
            if (k == 0) ykn = ykn + cc
	    if ( kfirst ) then
	      kfirst = .false.
              ! .. clear array for fyk
              fyk = 0
              ! .. clear the array yk
              yk = 0
		do 50 m = 1,ks
		  do 50 ii = 1,nint
		    fyk(ii,m) = fyk(ii,m) + cc*fc1(ii,m)*fc1(ii,m)
 50             continue
	      end if
	    end if
 30       continue
 
          .. if there is a contribution, add
 
          if (.not. kfirst) then
 
          add a relativistic correction if rel=TRUE
 
	    if (rel) then
	      c = (2*k+1)*fine
              do 60 m = 1,ks
	        do 60 ii = 1,nint
	          yr(ii,m)= yr(ii,m)+ c*fyk(ii,m)*grm(ii,m)
  60          continue
	    end if
 
    set up the array yk (spline approx. of SUM ()(1/r)fc1(r)fc2(r)
 
	    do 70 m=1,ks
	      do 70 ith=1,ks
	        do 70 ii = 1,nint
     	        irow=ii+ith-1
	          yk(irow)=yk(irow)+grw(ii,m)*grm(ii,m)*bsp(ii,m,ith)*
     :  	        fyk(ii,m)
   70	    continue
 
    	    const=-(2*k+1.d0)
	    do 80 ii=1,ns
	      yk(ii)=const*yk(ii)
   80	    continue
	    yk(1)=0.d0
	    yk(ns) = ykn
            print *, ' Potl: yk'
            print '(6f12.8)', (yk(ii),ii=1,ns)
 
    solve the matrix equation
 
	    call dgbsl(dyk,ktx,ns,ks-1,ks-1,ipvtd,yk,0)
            print *, 'Potl: after dgbsl', i
            print '(6f12.8)', (yk(ii),ii=1,ns)
 
             .. add to yv
	    do 90 ii = 1,ns
	      yv(ii) = yv(ii) + yk(ii)
 90         continue
	  end if
 20     continue 
       print *, 'Potl: yv', i
       print '(6f12.8)', (yv(ii),ii=1,ns)
 
   evaluates the function Y(r) at all the gaussian points
 
        call yval(1,0,nt,kx,ks,nint,grm,bsp,yv,yr)
	print *, ' Y(r)'
	do 100 m = 1,ks
	  do 100 ii = 1,nint
 	    yr(ii,m) = yr(ii,m)*grw(ii,m)*grm(ii,m)
 100    continue
	call minty(0,nt,kx,ks,nint,bsp,yr,yrm)
        print *, ' Test yrm'
	call bxv(nt,ks,ns,yrm,p(1,i),yv)
        ans1 = ddot(ns,yv,1,p(1,i),1)
	print *, '<p1,potl,p1>', ans1, quadr(1,1,-1)
 
	end
