*===================================================================
	subroutine minty(icase,nt,kx,ks,nint,bsp,ygr,ym)
*===================================================================
*
*    	Computes the array elements
*
*            <B_i,  y(r) B_j>.
*
* 	SUBROUTINES called
*	     none
*
* 	This version dated 11/3/1989
*
*-------------------------------------------------------------------
*
	implicit double precision (a-h,o-z)
	dimension ygr(nt,kx), bsp(nt,kx,kx), ym(nt,kx)
*
*       on entry
*       --------
*    icase 0 if i >=j (subdiagonals); 1 if i ,=j (superdiagonals)
*    nt    the leading dimension of array, must be greater than the
*          number of elements in the knot sequence
*    kx    maximum order of the spline
*    ks    the order of the spline (<=kx)
*    nint  the number of intervals starting at t(k)
*    bsp   array of B-spline values at the gaussian points
*    ygr   array of values of a specific function  y(r) at the
*          gaussian points of each interval
*
*    on exit
*    -------
*    ym    <B_i, y(r) B_j> in symmetric upper storage mode
*
*--------------------------------------------------------------------
	ns = nint+ks-1
*
*    clear the ym array
*
      do 10 m = 1,ks
	do 10 i = 1,ns
	  ym(i,m) = 0.d0
 10   continue
*     
      if (icase .ne. 0) then
        do 20 m = 1,ks
	  do 20 ith = 1,ks
	    do 20 jth = ith,ks
	      jcol = jth-ith+1
	      do 20 i = 1,nint
	        irow = i+ith-1
	        ym(irow,jcol) = ym(irow,jcol)+
     :                             ygr(i,m)*bsp(i,m,ith)*bsp(i,m,jth)
 20     continue 
      else
        do 30 m = 1,ks
	  do 30 ith = 1,ks
	    do 30 jth = 1,ith
	      jcol = jth-ith+ks
	      do 30 i = 1,nint
	        irow = i+ith-1
	        ym(irow,jcol) = ym(irow,jcol)+
     :                             ygr(i,m)*bsp(i,m,ith)*bsp(i,m,jth)
 30     continue 
      end if
*
       end 
