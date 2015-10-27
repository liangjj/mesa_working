*deck dwupdt
      subroutine dwupdt (n, r, ldr, w, b, alpha, cos, sin)
c***begin prologue  dwupdt
c***subsidiary
c***purpose  subsidiary to dnls1 and dnls1e
c***library   slatec
c***type      double precision (rwupdt-s, dwupdt-d)
c***author  (unknown)
c***description
c
c     given an n by n upper triangular matrix r, this subroutine
c     computes the qr decomposition of the matrix formed when a row
c     is added to r. if the row is specified by the vector w, then
c     dwupdt determines an orthogonal matrix q such that when the
c     n+1 by n matrix composed of r augmented by w is premultiplied
c     by (q transpose), the resulting matrix is upper trapezoidal.
c     the orthogonal matrix q is the product of n transformations
c
c           g(1)*g(2)* ... *g(n)
c
c     where g(i) is a givens rotation in the (i,n+1) plane which
c     eliminates elements in the i-th plane. dwupdt also
c     computes the product (q transpose)*c where c is the
c     (n+1)-vector (b,alpha). q itself is not accumulated, rather
c     the information to recover the g rotations is supplied.
c
c     the subroutine statement is
c
c       subroutine dwupdt(n,r,ldr,w,b,alpha,cos,sin)
c
c     where
c
c       n is a positive integer input variable set to the order of r.
c
c       r is an n by n array. on input the upper triangular part of
c         r must contain the matrix to be updated. on output r
c         contains the updated triangular matrix.
c
c       ldr is a positive integer input variable not less than n
c         which specifies the leading dimension of the array r.
c
c       w is an input array of length n which must contain the row
c         vector to be added to r.
c
c       b is an array of length n. on input b must contain the
c         first n elements of the vector c. on output b contains
c         the first n elements of the vector (q transpose)*c.
c
c       alpha is a variable. on input alpha must contain the
c         (n+1)-st element of the vector c. on output alpha contains
c         the (n+1)-st element of the vector (q transpose)*c.
c
c       cos is an output array of length n which contains the
c         cosines of the transforming givens rotations.
c
c       sin is an output array of length n which contains the
c         sines of the transforming givens rotations.
c
c     **********
c
c***see also  dnls1, dnls1e
c***routines called  (none)
c***revision history  (yymmdd)
c   800301  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900326  removed duplicate information from description section.
c           (wrb)
c   900328  added type section.  (wrb)
c***end prologue  dwupdt
      integer n,ldr
      double precision alpha
      double precision r(ldr,*),w(*),b(*),cos(*),sin(*)
      integer i,j,jm1
      double precision cotan,one,p5,p25,rowj,tan,temp,zero
      save one, p5, p25, zero
      data one,p5,p25,zero /1.0d0,5.0d-1,2.5d-1,0.0d0/
c***first executable statement  dwupdt
      do 60 j = 1, n
         rowj = w(j)
         jm1 = j - 1
c
c        apply the previous transformations to
c        r(i,j), i=1,2,...,j-1, and to w(j).
c
         if (jm1 .lt. 1) go to 20
         do 10 i = 1, jm1
            temp = cos(i)*r(i,j) + sin(i)*rowj
            rowj = -sin(i)*r(i,j) + cos(i)*rowj
            r(i,j) = temp
   10       continue
   20    continue
c
c        determine a givens rotation which eliminates w(j).
c
         cos(j) = one
         sin(j) = zero
         if (rowj .eq. zero) go to 50
         if (abs(r(j,j)) .ge. abs(rowj)) go to 30
            cotan = r(j,j)/rowj
            sin(j) = p5/sqrt(p25+p25*cotan**2)
            cos(j) = sin(j)*cotan
            go to 40
   30    continue
            tan = rowj/r(j,j)
            cos(j) = p5/sqrt(p25+p25*tan**2)
            sin(j) = cos(j)*tan
   40    continue
c
c        apply the current transformation to r(j,j), b(j), and alpha.
c
         r(j,j) = cos(j)*r(j,j) + sin(j)*rowj
         temp = cos(j)*b(j) + sin(j)*alpha
         alpha = -sin(j)*b(j) + cos(j)*alpha
         b(j) = temp
   50    continue
   60    continue
      return
c
c     last card of subroutine dwupdt.
c
      end
