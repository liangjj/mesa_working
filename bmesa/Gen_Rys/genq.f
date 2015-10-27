*deck @(#)genq.f	1.1 9/9/91
      subroutine genq( t, w, b, muzero, kpts, endpts, n)
c
c        b        real scratch array of length n
c
c     output parameters (both arrays of length n)
c
c        t        will contain the desired nodes x(1),,,x(n)
c        w        will contain the desired weights c(1),,,c(n)
c
      implicit real*8 (a-h,o-z)
      real*8  muzero
      dimension  b(n), t(n), w(n), endpts(2)
      common/io/inp,iout
c
c     the matrix of coefficients is assumed to be symmetric.
c     the array t contains the diagonal elements, the array
c     b the off-diagonal elements.
c     make appropriate changes in the lower right 2 by 2
c     submatrix if lobatto or radau.
c
      if(kpts.eq.1) then
c
c        if kpts=1, only t(n) must be changed
c
         t(n) =gbslve(endpts(1), n, t, b)*b(n-1)**2 + endpts(1)
      elseif(kpts.eq.2) then
c
c        if kpts=2, t(n) and b(n-1) must be recomputed
c
         gam =gbslve(endpts(1), n, t, b)
         t1 = ((endpts(1) - endpts(2))/(gbslve(endpts(2), n, t, b)
     1                    - gam))
         b(n-1) =  sqrt(t1)
         t(n) = endpts(1) + gam*t1
      endif
c
c     note that the indices of the elements of b run from 1 to n-1
c     and thus the value of b(n) is arbitrary.
c     now compute the eigenvalues of the symmetric tridiagonal
c     matrix, which has been modified as necessary.
c     the method used is a ql-type method with origin shifting
c
      call rzero(w,n)
      w(1) = 1.0d0
c
      call gbtql2 (n, t, b, w, ierr)
      do 10 i = 1, n
         w(i) = muzero * w(i) * w(i)
 10   continue   
c
      return
      end
