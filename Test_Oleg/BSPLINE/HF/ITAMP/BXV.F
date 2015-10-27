*=======================================================================
        subroutine bxv(nt,k,n,b,v,y)
*=======================================================================
*   Computes y = b * v    where b is a symmetric, banded matrix,
*   and v, y are vectors
*
*   Written by C. F. Fischer
*------------------------------------------------------------------------
*
	implicit double precision (a-h,o-z)
	dimension b(nt,*),v(*),y(*)
*
*   on entry
*   --------
*     nt      the leading dimension of arrays.
*     k       the number of diagonals
*     n       the order of the matrix
*     b       the symmetric, banded matrix in column storage mode
*     v       vector
*
*   on exit
*   -------
*     y       y = B*v 
*-----------------------------------------------------------------------
*
* ...   initialize the y array
*
*	  print *, 'nt,n,k',nt,n,k
*	  print *, 'Matrix'
*	  do 99 m = 1,k
*	    print '(6f12.8)', (b(ii,m),ii=1,n)
*99        continue
*	  print *, 'Vector'
*	  print '(6f12.8)', (v(ii),ii=1,n)

          do 1 j = 1,n
  	    y(j) = 0.d0
 1        continue
*
*      .. contribution from sub-diagonals
*
       do 10 i=1,k
         do 10 j=k-i+1,n
	   y(j) = y(j)  +b(j,i)*v(j-k+i)
10     continue
*
*      .. contribution from super-diagonals
*
       do 20 i=1,k-1
         do 20 j=1,n-k+i
	   y(j) = y(j)  +b(j+k-i,i)*v(j+k-i)
20     continue
*
	  print *, 'BXV Result'
	  print '(6f12.8)', (y(ii),ii=1,n)
       END
