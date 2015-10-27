*deck trisp
      subroutine trisp (n, a, b, c, d, u, z)
c***begin prologue  trisp
c***subsidiary
c***purpose  subsidiary to sepeli
c***library   slatec
c***type      single precision (trisp-s)
c***author  (unknown)
c***description
c
c     this subroutine solves for a non-zero eigenvector corresponding
c     to the zero eigenvalue of the transpose of the rank
c     deficient one matrix with subdiagonal a, diagonal b, and
c     superdiagonal c , with a(1) in the (1,n) position, with
c     c(n) in the (n,1) position, and all other elements zero.
c
c***see also  sepeli
c***routines called  (none)
c***revision history  (yymmdd)
c   801001  date written
c   890831  modified array declarations.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c***end prologue  trisp
c
      dimension       a(*)       ,b(*)       ,c(*)       ,d(*)       ,
     1                u(*)       ,z(*)
c***first executable statement  trisp
      bn = b(n)
      d(1) = a(2)/b(1)
      v = a(1)
      u(1) = c(n)/b(1)
      nm2 = n-2
      do  10 j=2,nm2
         den = b(j)-c(j-1)*d(j-1)
         d(j) = a(j+1)/den
         u(j) = -c(j-1)*u(j-1)/den
         bn = bn-v*u(j-1)
         v = -v*d(j-1)
   10 continue
      den = b(n-1)-c(n-2)*d(n-2)
      d(n-1) = (a(n)-c(n-2)*u(n-2))/den
      an = c(n-1)-v*d(n-2)
      bn = bn-v*u(n-2)
      den = bn-an*d(n-1)
c
c     set last component equal to one
c
      z(n) = 1.0
      z(n-1) = -d(n-1)
      nm1 = n-1
      do  20 j=2,nm1
         k = n-j
         z(k) = -d(k)*z(k+1)-u(k)*z(n)
   20 continue
      return
      end
