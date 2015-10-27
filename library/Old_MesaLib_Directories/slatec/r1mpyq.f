*deck r1mpyq
      subroutine r1mpyq (m, n, a, lda, v, w)
c***begin prologue  r1mpyq
c***subsidiary
c***purpose  subsidiary to snsq and snsqe
c***library   slatec
c***type      single precision (r1mpyq-s, d1mpyq-d)
c***author  (unknown)
c***description
c
c     given an m by n matrix a, this subroutine computes a*q where
c     q is the product of 2*(n - 1) transformations
c
c           gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)
c
c     and gv(i), gw(i) are givens rotations in the (i,n) plane which
c     eliminate elements in the i-th and n-th planes, respectively.
c     q itself is not given, rather the information to recover the
c     gv, gw rotations is supplied.
c
c     the subroutine statement is
c
c       subroutine r1mpyq(m,n,a,lda,v,w)
c
c     where
c
c       m is a positive integer input variable set to the number
c         of rows of a.
c
c       n is a positive integer input variable set to the number
c         of columns of a.
c
c       a is an m by n array. on input a must contain the matrix
c         to be postmultiplied by the orthogonal matrix q
c         described above. on output a*q has replaced a.
c
c       lda is a positive integer input variable not less than m
c         which specifies the leading dimension of the array a.
c
c       v is an input array of length n. v(i) must contain the
c         information necessary to recover the givens rotation gv(i)
c         described above.
c
c       w is an input array of length n. w(i) must contain the
c         information necessary to recover the givens rotation gw(i)
c         described above.
c
c***see also  snsq, snsqe
c***routines called  (none)
c***revision history  (yymmdd)
c   800301  date written
c   890831  modified array declarations.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900326  removed duplicate information from description section.
c           (wrb)
c   900328  added type section.  (wrb)
c***end prologue  r1mpyq
      integer m,n,lda
      real a(lda,*),v(*),w(*)
      integer i,j,nmj,nm1
      real cos,one,sin,temp
      save one
      data one /1.0e0/
c***first executable statement  r1mpyq
      nm1 = n - 1
      if (nm1 .lt. 1) go to 50
      do 20 nmj = 1, nm1
         j = n - nmj
         if (abs(v(j)) .gt. one) cos = one/v(j)
         if (abs(v(j)) .gt. one) sin = sqrt(one-cos**2)
         if (abs(v(j)) .le. one) sin = v(j)
         if (abs(v(j)) .le. one) cos = sqrt(one-sin**2)
         do 10 i = 1, m
            temp = cos*a(i,j) - sin*a(i,n)
            a(i,n) = sin*a(i,j) + cos*a(i,n)
            a(i,j) = temp
   10       continue
   20    continue
c
c     apply the second set of givens rotations to a.
c
      do 40 j = 1, nm1
         if (abs(w(j)) .gt. one) cos = one/w(j)
         if (abs(w(j)) .gt. one) sin = sqrt(one-cos**2)
         if (abs(w(j)) .le. one) sin = w(j)
         if (abs(w(j)) .le. one) cos = sqrt(one-sin**2)
         do 30 i = 1, m
            temp = cos*a(i,j) + sin*a(i,n)
            a(i,n) = -sin*a(i,j) + cos*a(i,n)
            a(i,j) = temp
   30       continue
   40    continue
   50 continue
      return
c
c     last card of subroutine r1mpyq.
c
      end
