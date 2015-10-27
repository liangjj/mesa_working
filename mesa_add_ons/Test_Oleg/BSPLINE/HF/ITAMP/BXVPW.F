*=======================================================================
        subroutine bxvpw(nt,k,n,c,b,v,y,w)
*=======================================================================
*   Computes y = c* b * v + y   where b is a symmetric, banded matrix
*   and v, w are vectors
*
*   Written by C. F. Fischer
*------------------------------------------------------------------------
*
	implicit double precision (a-h,o-z)
	dimension b(nt,*),v(*),y(*),w(*)
*
*   on entry
*   --------
*     nt      the leading dimension of arrays.
*     k       the number of diagonals
*     n       the order of the matrix
*     c,ic    coefficients
*     b       the symmetric, banded matrix in column storage mode
*     v       vector
*     w       working array
*
*   on exit
*   -------
*     y       y = c*B*v +y
*-----------------------------------------------------------------------

*
* ...   initialize the w array
*
       do 1 j = 1,n
	 w(j) = 0.d0
 1     continue
*
*      .. contribution from sub-diagonals
*
       do 10 i=1,k
         do 10 j=k-i+1,n
	   w(j) = w(j)  +b(j,i)*v(j-k+i)
10     continue
*
*      .. contribution from super-diagonals
*
       do 20 i=1,k-1
         do 20 j=1,n-k+i
	   w(j) = w(j)  +b(j+k-i,i)*v(j+k-i)
20     continue
*
*      print *, ' ... y before call to bxvpw'
*      print *, (y(i),i=1,n)
       print *, ' C was',c
       if ( c .ne. 1.d0) then
         do 50 i = 1,n
           y(i) = y(i) + c*w(i)
   50    Continue
       else
         Do 70 i = 1,n
           y(i) = y(i) + w(i)
   70    Continue
       end if
*      print *, ' ... y after call to bxvpw'
*      print *, (y(i),i=1,n)
       END
