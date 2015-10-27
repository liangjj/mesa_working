*deck dppco
      subroutine dppco (ap, n, rcond, z, info)
c***begin prologue  dppco
c***purpose  factor a symmetric positive definite matrix stored in
c            packed form and estimate the condition number of the
c            matrix.
c***library   slatec (linpack)
c***category  d2b1b
c***type      double precision (sppco-s, dppco-d, cppco-c)
c***keywords  condition number, linear algebra, linpack,
c             matrix factorization, packed, positive definite
c***author  moler, c. b., (u. of new mexico)
c***description
c
c     dppco factors a double precision symmetric positive definite
c     matrix stored in packed form
c     and estimates the condition of the matrix.
c
c     if  rcond  is not needed, dppfa is slightly faster.
c     to solve  a*x = b , follow dppco by dppsl.
c     to compute  inverse(a)*c , follow dppco by dppsl.
c     to compute  determinant(a) , follow dppco by dppdi.
c     to compute  inverse(a) , follow dppco by dppdi.
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
c                if  info .ne. 0 , the factorization is not complete.
c
c        rcond   double precision
c                an estimate of the reciprocal condition of  a .
c                for the system  a*x = b , relative perturbations
c                in  a  and  b  of size  epsilon  may cause
c                relative perturbations in  x  of size  epsilon/rcond .
c                if  rcond  is so small that the logical expression
c                           1.0 + rcond .eq. 1.0
c                is true, then  a  may be singular to working
c                precision.  in particular,  rcond  is zero  if
c                exact singularity is detected or the estimate
c                underflows.  if info .ne. 0 , rcond is unchanged.
c
c        z       double precision(n)
c                a work vector whose contents are usually unimportant.
c                if  a  is singular to working precision, then  z  is
c                an approximate null vector in the sense that
c                norm(a*z) = rcond*norm(a)*norm(z) .
c                if  info .ne. 0 , z  is unchanged.
c
c        info    integer
c                = 0  for normal return.
c                = k  signals an error condition.  the leading minor
c                     of order  k  is not positive definite.
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
c***routines called  dasum, daxpy, ddot, dppfa, dscal
c***revision history  (yymmdd)
c   780814  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900326  removed duplicate information from description section.
c           (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  dppco
      integer n,info
      double precision ap(*),z(*)
      double precision rcond
c
      double precision ddot,ek,t,wk,wkm
      double precision anorm,s,dasum,sm,ynorm
      integer i,ij,j,jm1,j1,k,kb,kj,kk,kp1
c
c     find norm of a
c
c***first executable statement  dppco
      j1 = 1
      do 30 j = 1, n
         z(j) = dasum(j,ap(j1),1)
         ij = j1
         j1 = j1 + j
         jm1 = j - 1
         if (jm1 .lt. 1) go to 20
         do 10 i = 1, jm1
            z(i) = z(i) + abs(ap(ij))
            ij = ij + 1
   10    continue
   20    continue
   30 continue
      anorm = 0.0d0
      do 40 j = 1, n
         anorm = max(anorm,z(j))
   40 continue
c
c     factor
c
      call dppfa(ap,n,info)
      if (info .ne. 0) go to 180
c
c        rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) .
c        estimate = norm(z)/norm(y) where  a*z = y  and  a*y = e .
c        the components of  e  are chosen to cause maximum local
c        growth in the elements of w  where  trans(r)*w = e .
c        the vectors are frequently rescaled to avoid overflow.
c
c        solve trans(r)*w = e
c
         ek = 1.0d0
         do 50 j = 1, n
            z(j) = 0.0d0
   50    continue
         kk = 0
         do 110 k = 1, n
            kk = kk + k
            if (z(k) .ne. 0.0d0) ek = sign(ek,-z(k))
            if (abs(ek-z(k)) .le. ap(kk)) go to 60
               s = ap(kk)/abs(ek-z(k))
               call dscal(n,s,z,1)
               ek = s*ek
   60       continue
            wk = ek - z(k)
            wkm = -ek - z(k)
            s = abs(wk)
            sm = abs(wkm)
            wk = wk/ap(kk)
            wkm = wkm/ap(kk)
            kp1 = k + 1
            kj = kk + k
            if (kp1 .gt. n) go to 100
               do 70 j = kp1, n
                  sm = sm + abs(z(j)+wkm*ap(kj))
                  z(j) = z(j) + wk*ap(kj)
                  s = s + abs(z(j))
                  kj = kj + j
   70          continue
               if (s .ge. sm) go to 90
                  t = wkm - wk
                  wk = wkm
                  kj = kk + k
                  do 80 j = kp1, n
                     z(j) = z(j) + t*ap(kj)
                     kj = kj + j
   80             continue
   90          continue
  100       continue
            z(k) = wk
  110    continue
         s = 1.0d0/dasum(n,z,1)
         call dscal(n,s,z,1)
c
c        solve r*y = w
c
         do 130 kb = 1, n
            k = n + 1 - kb
            if (abs(z(k)) .le. ap(kk)) go to 120
               s = ap(kk)/abs(z(k))
               call dscal(n,s,z,1)
  120       continue
            z(k) = z(k)/ap(kk)
            kk = kk - k
            t = -z(k)
            call daxpy(k-1,t,ap(kk+1),1,z(1),1)
  130    continue
         s = 1.0d0/dasum(n,z,1)
         call dscal(n,s,z,1)
c
         ynorm = 1.0d0
c
c        solve trans(r)*v = y
c
         do 150 k = 1, n
            z(k) = z(k) - ddot(k-1,ap(kk+1),1,z(1),1)
            kk = kk + k
            if (abs(z(k)) .le. ap(kk)) go to 140
               s = ap(kk)/abs(z(k))
               call dscal(n,s,z,1)
               ynorm = s*ynorm
  140       continue
            z(k) = z(k)/ap(kk)
  150    continue
         s = 1.0d0/dasum(n,z,1)
         call dscal(n,s,z,1)
         ynorm = s*ynorm
c
c        solve r*z = v
c
         do 170 kb = 1, n
            k = n + 1 - kb
            if (abs(z(k)) .le. ap(kk)) go to 160
               s = ap(kk)/abs(z(k))
               call dscal(n,s,z,1)
               ynorm = s*ynorm
  160       continue
            z(k) = z(k)/ap(kk)
            kk = kk - k
            t = -z(k)
            call daxpy(k-1,t,ap(kk+1),1,z(1),1)
  170    continue
c        make znorm = 1.0
         s = 1.0d0/dasum(n,z,1)
         call dscal(n,s,z,1)
         ynorm = s*ynorm
c
         if (anorm .ne. 0.0d0) rcond = ynorm/anorm
         if (anorm .eq. 0.0d0) rcond = 0.0d0
  180 continue
      return
      end
