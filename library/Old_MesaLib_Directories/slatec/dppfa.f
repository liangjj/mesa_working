*deck dppfa
      subroutine dppfa (ap, n, info)
c***begin prologue  dppfa
c***purpose  factor a real symmetric positive definite matrix stored in
c            packed form.
c***library   slatec (linpack)
c***category  d2b1b
c***type      double precision (sppfa-s, dppfa-d, cppfa-c)
c***keywords  linear algebra, linpack, matrix factorization, packed,
c             positive definite
c***author  moler, c. b., (u. of new mexico)
c***description
c
c     dppfa factors a double precision symmetric positive definite
c     matrix stored in packed form.
c
c     dppfa is usually called by dppco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for dppco) = (1 + 18/n)*(time for dppfa) .
c
c     on entry
c
c        ap      double precision (n*(n+1)/2)
c                the packed form of a symmetric matrix  a .  the
c                columns of the upper triangle are stored sequentially
c                in a one-dimensional array of length  n*(n+1)/2 .
c                see comments below for details.
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        ap      an upper triangular matrix  r , stored in packed
c                form, so that  a = trans(r)*r .
c
c        info    integer
c                = 0  for normal return.
c                = k  if the leading minor of order  k  is not
c                     positive definite.
c
c
c     packed storage
c
c          the following program segment will pack the upper
c          triangle of a symmetric matrix.
c
c                k = 0
c                do 20 j = 1, n
c                   do 10 i = 1, j
c                      k = k + 1
c                      ap(k) = a(i,j)
c             10    continue
c             20 continue
c
c***references  j. j. dongarra, j. r. bunch, c. b. moler, and g. w.
c                 stewart, linpack users' guide, siam, 1979.
c***routines called  ddot
c***revision history  (yymmdd)
c   780814  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900326  removed duplicate information from description section.
c           (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  dppfa
      integer n,info
      double precision ap(*)
c
      double precision ddot,t
      double precision s
      integer j,jj,jm1,k,kj,kk
c***first executable statement  dppfa
         jj = 0
         do 30 j = 1, n
            info = j
            s = 0.0d0
            jm1 = j - 1
            kj = jj
            kk = 0
            if (jm1 .lt. 1) go to 20
            do 10 k = 1, jm1
               kj = kj + 1
               t = ap(kj) - ddot(k-1,ap(kk+1),1,ap(jj+1),1)
               kk = kk + k
               t = t/ap(kk)
               ap(kj) = t
               s = s + t*t
   10       continue
   20       continue
            jj = jj + j
            s = ap(jj) - s
            if (s .le. 0.0d0) go to 40
            ap(jj) = sqrt(s)
   30    continue
         info = 0
   40 continue
      return
      end
