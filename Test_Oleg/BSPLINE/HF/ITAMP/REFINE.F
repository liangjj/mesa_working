*=======================================================================
        subroutine refine(nt,n,a,b,c,eigval,eigvec,ipvt,x,ierr)
*=======================================================================
*   Given an estimate of the eigenvalue, 
*           (a - lambda b)v = c v = 0
*   computes and eigenvector and updates the eigenvalue
*
*   Written by C. F. Fischer
*------------------------------------------------------------------------
*
	implicit double precision (a-h,o-z)
	dimension a(nt,nt), b(nt,nt), c(nt,nt),
     :            eigvec(nt),ipvt(nt),x(nt)
*
*   on entry
*   --------
*     nt     dimensions of a, b
*     n       size of the matrix
*     a,b     matrices
*     eigval  estimate of lambda
*     ipvtd   integer work array
*     x       double precision work array
*
*   on exit
*   -------
*     eigval  improved estimate of the eigenvalue
*     eigvec  the eigenvector corresponding to eigval
*
*-----------------------------------------------------------------------
*
* 	..Form the c matrix
	s = eigval
	do 10 j = 1,n
	   do 10 i = 1,n
	     c(i,j) = a(i,j) -s*b(i,j)
   10   continue
*
*       .. form an initial guess of the eigenvector
	do 20 i = 1,n
	  eigvec(i) = 0.d0
   20   continue
	eigvec(2) = .5d0
	eigvec(3) = .5d0
	eigvec(4) = .5d0
	eigvec(5) = 1.d0

   30   call ginvit(c,nt,n,eigvec,ipvt,x,ierr)
*
	eiv = 0.d0
	den = 0.d0
	do 50 i = 1,n
	  eiv = eiv + a(i,i)*eigvec(i)**2
	  den = den + b(i,i)*eigvec(i)**2
   50   continue
*
*       .. the off-diagonal part
	eiv2 = 0.d0
	den2 = 0.d0
	do 60 j = 1,n-1
	  t = 0.d0
	  d = 0.d0
	  do 70 i = j+1,n
             t= t + a(i,j)*eigvec(i)
	     d= d + b(i,j)*eigvec(i)
   70     continue
	  eiv2 = eiv2 + t*eigvec(j)
	  den2 = den2 + d*eigvec(j)
   60   continue
 
	eivnew = (eiv + 2*eiv2)/(den +2*den2)
	eigval = eivnew
*
*       .. normalize the eigenvector
	dn = 1.d0/sqrt(den+2*den2)
	do 80 i = 1,n
	  eigvec(i) = eigvec(i)*dn
   80   continue

	return
	end
