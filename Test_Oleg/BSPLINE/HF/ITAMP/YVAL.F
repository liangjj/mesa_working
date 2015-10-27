*=======================================================================
        subroutine yval(iadd,mm,nt,kx,ks,nint,gr,bsp,yv,ygr)
*=======================================================================
*
*   This routine computes the values of r^mm f(r)
*   at the gaussian points of each interval,  where f(r) is defined by
*   the spline expansion vector yv. 
*
*   VBSP. this version dated 8/10/89
*------------------------------------------------------------------------
	implicit double precision (a-h,o-z)
	dimension gr(nt,kx),bsp(nt,kx,kx),yv(nt),ygr(nt,kx)
*
*   on entry
*   --------
*     iadd    if iadd > 0, results are added to previous
*     mm      integer defining the power of r
*     nt      the leading dimension of arrays, must be greater than
*             the number of elements in the knot sequence
*     kx      maximum order of the spline
*     ks      the order of the spline (<= kx)
*     nint    the number of intervals starting at t(k)
*     gr      the values of r at the gaussian points
*     bsp     array of B-spline values at the gaussian points
*     yv      the spline expansion vector for the funciton f(r)
*
*
*   on exit
*   -------
*     ygr     array of values of y(x) = Sum_i c_i B_i(x) at the
*             gaussian points of each interval
*
*  Note:  Negative values of mm can be avoided by replacing the gr
*     array in the argument by the array of reciprocals, grm
*-----------------------------------------------------------------------

*
*  ...  clear the ygr array
*   
*     print *, 'Yval: iadd,mm,ks,nint',iadd,mm,ks,nint
*     print *, 'yv'
*     print '(6f12.8)', (yv(i),i=1,ns)
      if (iadd .le. 0) then
        do 10 m=1,ks
          do 10 i=1,nint
            ygr(i,m)=0.d0
 10     continue
      end if
*
* ...  form the values
*
      if (mm .gt. 1 .or. mm .lt. 0) then
	do 20 m = 1,ks
	   do 20 ith = 1,ks
	      do 20 i = 1,nint
		ygr(i,m)=ygr(i,m)+
     :                   yv(i+ith-1)*bsp(i,m,ith)*gr(i,m)**mm
 20     continue
*    
      else if (mm .eq. 1) then
	do 30 m = 1,ks
	   do 30 ith = 1,ks
	      do 30 i = 1,nint
		ygr(i,m)=ygr(i,m)+
     :                   yv(i+ith-1)*bsp(i,m,ith)*gr(i,m)
 30     continue
      else if (mm .eq. 0) then
	do 40 m = 1,ks
	   do 40 ith = 1,ks
	      do 40 i = 1,nint
		ygr(i,m)=ygr(i,m)+yv(i+ith-1)*bsp(i,m,ith)
 40     continue
      end if
      return
      end
