*deck sptsl
      subroutine sptsl (n, d, e, b)
c***begin prologue  sptsl
c***purpose  solve a positive definite tridiagonal linear system.
c***library   slatec (linpack)
c***category  d2b2a
c***type      single precision (sptsl-s, dptsl-d, cptsl-c)
c***keywords  linear algebra, linpack, matrix, positive definite, solve,
c             tridiagonal
c***author  dongarra, j., (anl)
c***description
c
c     sptsl given a positive definite tridiagonal matrix and a right
c     hand side will find the solution.
c
c     on entry
c
c        n        integer
c                 is the order of the tridiagonal matrix.
c
c        d        real(n)
c                 is the diagonal of the tridiagonal matrix.
c                 on output, d is destroyed.
c
c        e        real(n)
c                 is the offdiagonal of the tridiagonal matrix.
c                 e(1) through e(n-1) should contain the
c                 offdiagonal.
c
c        b        real(n)
c                 is the right hand side vector.
c
c     on return
c
c        b        contains the solution.
c
c***references  j. j. dongarra, j. r. bunch, c. b. moler, and g. w.
c                 stewart, linpack users' guide, siam, 1979.
c***routines called  (none)
c***revision history  (yymmdd)
c   780814  date written
c   890505  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900326  removed duplicate information from description section.
c           (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  sptsl
      integer n
      real d(*),e(*),b(*)
c
      integer k,kbm1,ke,kf,kp1,nm1,nm1d2
      real t1,t2
c
c     check for 1 x 1 case
c
c***first executable statement  sptsl
      if (n .ne. 1) go to 10
         b(1) = b(1)/d(1)
      go to 70
   10 continue
         nm1 = n - 1
         nm1d2 = nm1/2
         if (n .eq. 2) go to 30
            kbm1 = n - 1
c
c           zero top half of subdiagonal and bottom half of
c           superdiagonal
c
            do 20 k = 1, nm1d2
               t1 = e(k)/d(k)
               d(k+1) = d(k+1) - t1*e(k)
               b(k+1) = b(k+1) - t1*b(k)
               t2 = e(kbm1)/d(kbm1+1)
               d(kbm1) = d(kbm1) - t2*e(kbm1)
               b(kbm1) = b(kbm1) - t2*b(kbm1+1)
               kbm1 = kbm1 - 1
   20       continue
   30    continue
         kp1 = nm1d2 + 1
c
c        clean up for possible 2 x 2 block at center
c
         if (mod(n,2) .ne. 0) go to 40
            t1 = e(kp1)/d(kp1)
            d(kp1+1) = d(kp1+1) - t1*e(kp1)
            b(kp1+1) = b(kp1+1) - t1*b(kp1)
            kp1 = kp1 + 1
   40    continue
c
c        back solve starting at the center, going towards the top
c        and bottom
c
         b(kp1) = b(kp1)/d(kp1)
         if (n .eq. 2) go to 60
            k = kp1 - 1
            ke = kp1 + nm1d2 - 1
            do 50 kf = kp1, ke
               b(k) = (b(k) - e(k)*b(k+1))/d(k)
               b(kf+1) = (b(kf+1) - e(kf)*b(kf))/d(kf+1)
               k = k - 1
   50       continue
   60    continue
         if (mod(n,2) .eq. 0) b(1) = (b(1) - e(1)*b(2))/d(1)
   70 continue
      return
      end
