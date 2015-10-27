*deck cppdi
      subroutine cppdi (ap, n, det, job)
c***begin prologue  cppdi
c***purpose  compute the determinant and inverse of a complex hermitian
c            positive definite matrix using factors from cppco or cppfa.
c***library   slatec (linpack)
c***category  d2d1b, d3d1b
c***type      complex (sppdi-s, dppdi-d, cppdi-c)
c***keywords  determinant, inverse, linear algebra, linpack, matrix,
c             packed, positive definite
c***author  moler, c. b., (u. of new mexico)
c***description
c
c     cppdi computes the determinant and inverse
c     of a complex hermitian positive definite matrix
c     using the factors computed by cppco or cppfa .
c
c     on entry
c
c        ap      complex (n*(n+1)/2)
c                the output from cppco or cppfa.
c
c        n       integer
c                the order of the matrix  a .
c
c        job     integer
c                = 11   both determinant and inverse.
c                = 01   inverse only.
c                = 10   determinant only.
c
c     on return
c
c        ap      the upper triangular half of the inverse .
c                the strict lower triangle is unaltered.
c
c        det     real(2)
c                determinant of original matrix if requested.
c                otherwise not referenced.
c                determinant = det(1) * 10.0**det(2)
c                with  1.0 .le. det(1) .lt. 10.0
c                or  det(1) .eq. 0.0 .
c
c     error condition
c
c        a division by zero will occur if the input factor contains
c        a zero on the diagonal and the inverse is requested.
c        it will not occur if the subroutines are called correctly
c        and if cpoco or cpofa has set info .eq. 0 .
c
c***references  j. j. dongarra, j. r. bunch, c. b. moler, and g. w.
c                 stewart, linpack users' guide, siam, 1979.
c***routines called  caxpy, cscal
c***revision history  (yymmdd)
c   780814  date written
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900326  removed duplicate information from description section.
c           (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  cppdi
      integer n,job
      complex ap(*)
      real det(2)
c
      complex t
      real s
      integer i,ii,j,jj,jm1,j1,k,kj,kk,kp1,k1
c***first executable statement  cppdi
c
c     compute determinant
c
      if (job/10 .eq. 0) go to 70
         det(1) = 1.0e0
         det(2) = 0.0e0
         s = 10.0e0
         ii = 0
         do 50 i = 1, n
            ii = ii + i
            det(1) = real(ap(ii))**2*det(1)
            if (det(1) .eq. 0.0e0) go to 60
   10       if (det(1) .ge. 1.0e0) go to 20
               det(1) = s*det(1)
               det(2) = det(2) - 1.0e0
            go to 10
   20       continue
   30       if (det(1) .lt. s) go to 40
               det(1) = det(1)/s
               det(2) = det(2) + 1.0e0
            go to 30
   40       continue
   50    continue
   60    continue
   70 continue
c
c     compute inverse(r)
c
      if (mod(job,10) .eq. 0) go to 140
         kk = 0
         do 100 k = 1, n
            k1 = kk + 1
            kk = kk + k
            ap(kk) = (1.0e0,0.0e0)/ap(kk)
            t = -ap(kk)
            call cscal(k-1,t,ap(k1),1)
            kp1 = k + 1
            j1 = kk + 1
            kj = kk + k
            if (n .lt. kp1) go to 90
            do 80 j = kp1, n
               t = ap(kj)
               ap(kj) = (0.0e0,0.0e0)
               call caxpy(k,t,ap(k1),1,ap(j1),1)
               j1 = j1 + j
               kj = kj + j
   80       continue
   90       continue
  100    continue
c
c        form  inverse(r) * ctrans(inverse(r))
c
         jj = 0
         do 130 j = 1, n
            j1 = jj + 1
            jj = jj + j
            jm1 = j - 1
            k1 = 1
            kj = j1
            if (jm1 .lt. 1) go to 120
            do 110 k = 1, jm1
               t = conjg(ap(kj))
               call caxpy(k,t,ap(j1),1,ap(k1),1)
               k1 = k1 + k
               kj = kj + 1
  110       continue
  120       continue
            t = conjg(ap(jj))
            call cscal(j,t,ap(j1),1)
  130    continue
  140 continue
      return
      end