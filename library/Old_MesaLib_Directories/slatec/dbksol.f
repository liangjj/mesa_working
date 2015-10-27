*deck dbksol
      subroutine dbksol (n, a, x)
c***begin prologue  dbksol
c***subsidiary
c***purpose  subsidiary to dbvsup
c***library   slatec
c***type      double precision (bksol-s, dbksol-d)
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
c***see also  dbvsup
c***routines called  ddot
c***revision history  (yymmdd)
c   750601  date written
c   890831  modified array declarations.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910722  updated author section.  (als)
c***end prologue  dbksol
c
      double precision ddot
      integer j, k, m, n, nm1
      double precision a(*), x(*)
c
c***first executable statement  dbksol
      m = (n*(n + 1))/2
      x(n) = x(n)*a(m)
      nm1 = n - 1
      if (nm1 .lt. 1) go to 20
      do 10 k = 1, nm1
         j = n - k
         m = m - k - 1
         x(j) = x(j)*a(m) - ddot(k,a(m+1),1,x(j+1),1)
   10 continue
   20 continue
c
      return
      end
