*deck bksol
      subroutine bksol (n, a, x)
c***begin prologue  bksol
c***subsidiary
c***purpose  subsidiary to bvsup
c***library   slatec
c***type      single precision (bksol-s, dbksol-d)
c***author  watts, h. a., (snla)
c***description
c
c **********************************************************************
c     solution of an upper triangular linear system by
c     back-substitution
c
c     the matrix a is assumed to be stored in a linear
c     array proceeding in a row-wise manner. the
c     vector x contains the given constant vector on input
c     and contains the solution on return.
c     the actual diagonal of a is unity while a diagonal
c     scaling matrix is stored there.
c **********************************************************************
c
c***see also  bvsup
c***routines called  sdot
c***revision history  (yymmdd)
c   750601  date written
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910722  updated author section.  (als)
c***end prologue  bksol
c
      dimension a(*),x(*)
c
c***first executable statement  bksol
      m=(n*(n+1))/2
      x(n)=x(n)*a(m)
      if (n .eq. 1) go to 20
      nm1=n-1
      do 10 k=1,nm1
      j=n-k
      m=m-k-1
   10 x(j)=x(j)*a(m) - sdot(k,a(m+1),1,x(j+1),1)
c
   20 return
      end
