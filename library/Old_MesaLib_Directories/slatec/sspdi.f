*deck sspdi
      subroutine sspdi (ap, n, kpvt, det, inert, work, job)
c***begin prologue  sspdi
c***purpose  compute the determinant, inertia, inverse of a real
c            symmetric matrix stored in packed form using the factors
c            from sspfa.
c***library   slatec (linpack)
c***category  d2b1a, d3b1a
c***type      single precision (sspdi-s, dspdi-d, chpdi-c, cspdi-c)
c***keywords  determinant, inverse, linear algebra, linpack, matrix,
c             packed, symmetric
c***author  bunch, j., (ucsd)
c***description
c
c     sspdi computes the determinant, inertia and inverse
c     of a real symmetric matrix using the factors from sspfa,
c     where the matrix is stored in packed form.
c
c     on entry
c
c        ap      real (n*(n+1)/2)
c                the output from sspfa.
c
c        n       integer
c                the order of the matrix a.
c
c        kpvt    integer(n)
c                the pivot vector from sspfa.
c
c        work    real(n)
c                work vector.  contents ignored.
c
c        job     integer
c                job has the decimal expansion  abc  where
c                   if  c .ne. 0, the inverse is computed,
c                   if  b .ne. 0, the determinant is computed,
c                   if  a .ne. 0, the inertia is computed.
c
c                for example, job = 111  gives all three.
c
c     on return
c
c        variables not requested by job are not used.
c
c        ap     contains the upper triangle of the inverse of
c               the original matrix, stored in packed form.
c               the columns of the upper triangle are stored
c               sequentially in a one-dimensional array.
c
c        det    real(2)
c               determinant of original matrix.
c               determinant = det(1) * 10.0**det(2)
c               with 1.0 .le. abs(det(1)) .lt. 10.0
c               or det(1) = 0.0.
c
c        inert  integer(3)
c               the inertia of the original matrix.
c               inert(1)  =  number of positive eigenvalues.
c               inert(2)  =  number of negative eigenvalues.
c               inert(3)  =  number of zero eigenvalues.
c
c     error condition
c
c        a division by zero will occur if the inverse is requested
c        and  sspco  has set rcond .eq. 0.0
c        or  sspfa  has set  info .ne. 0 .
c
c***references  j. j. dongarra, j. r. bunch, c. b. moler, and g. w.
c                 stewart, linpack users' guide, siam, 1979.
c***routines called  saxpy, scopy, sdot, sswap
c***revision history  (yymmdd)
c   780814  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   891107  modified routine equivalence list.  (wrb)
c   891107  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900326  removed duplicate information from description section.
c           (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  sspdi
      integer n,job
      real ap(*),work(*)
      real det(2)
      integer kpvt(*),inert(3)
c
      real akkp1,sdot,temp
      real ten,d,t,ak,akp1
      integer ij,ik,ikp1,iks,j,jb,jk,jkp1
      integer k,kk,kkp1,km1,ks,ksj,kskp1,kstep
      logical noinv,nodet,noert
c***first executable statement  sspdi
      noinv = mod(job,10) .eq. 0
      nodet = mod(job,100)/10 .eq. 0
      noert = mod(job,1000)/100 .eq. 0
c
      if (nodet .and. noert) go to 140
         if (noert) go to 10
            inert(1) = 0
            inert(2) = 0
            inert(3) = 0
   10    continue
         if (nodet) go to 20
            det(1) = 1.0e0
            det(2) = 0.0e0
            ten = 10.0e0
   20    continue
         t = 0.0e0
         ik = 0
         do 130 k = 1, n
            kk = ik + k
            d = ap(kk)
c
c           check if 1 by 1
c
            if (kpvt(k) .gt. 0) go to 50
c
c              2 by 2 block
c              use det (d  s)  =  (d/t * c - t) * t  ,  t = abs(s)
c                      (s  c)
c              to avoid underflow/overflow troubles.
c              take two passes through scaling.  use  t  for flag.
c
               if (t .ne. 0.0e0) go to 30
                  ikp1 = ik + k
                  kkp1 = ikp1 + k
                  t = abs(ap(kkp1))
                  d = (d/t)*ap(kkp1+1) - t
               go to 40
   30          continue
                  d = t
                  t = 0.0e0
   40          continue
   50       continue
c
            if (noert) go to 60
               if (d .gt. 0.0e0) inert(1) = inert(1) + 1
               if (d .lt. 0.0e0) inert(2) = inert(2) + 1
               if (d .eq. 0.0e0) inert(3) = inert(3) + 1
   60       continue
c
            if (nodet) go to 120
               det(1) = d*det(1)
               if (det(1) .eq. 0.0e0) go to 110
   70             if (abs(det(1)) .ge. 1.0e0) go to 80
                     det(1) = ten*det(1)
                     det(2) = det(2) - 1.0e0
                  go to 70
   80             continue
   90             if (abs(det(1)) .lt. ten) go to 100
                     det(1) = det(1)/ten
                     det(2) = det(2) + 1.0e0
                  go to 90
  100             continue
  110          continue
  120       continue
            ik = ik + k
  130    continue
  140 continue
c
c     compute inverse(a)
c
      if (noinv) go to 270
         k = 1
         ik = 0
  150    if (k .gt. n) go to 260
            km1 = k - 1
            kk = ik + k
            ikp1 = ik + k
            kkp1 = ikp1 + k
            if (kpvt(k) .lt. 0) go to 180
c
c              1 by 1
c
               ap(kk) = 1.0e0/ap(kk)
               if (km1 .lt. 1) go to 170
                  call scopy(km1,ap(ik+1),1,work,1)
                  ij = 0
                  do 160 j = 1, km1
                     jk = ik + j
                     ap(jk) = sdot(j,ap(ij+1),1,work,1)
                     call saxpy(j-1,work(j),ap(ij+1),1,ap(ik+1),1)
                     ij = ij + j
  160             continue
                  ap(kk) = ap(kk) + sdot(km1,work,1,ap(ik+1),1)
  170          continue
               kstep = 1
            go to 220
  180       continue
c
c              2 by 2
c
               t = abs(ap(kkp1))
               ak = ap(kk)/t
               akp1 = ap(kkp1+1)/t
               akkp1 = ap(kkp1)/t
               d = t*(ak*akp1 - 1.0e0)
               ap(kk) = akp1/d
               ap(kkp1+1) = ak/d
               ap(kkp1) = -akkp1/d
               if (km1 .lt. 1) go to 210
                  call scopy(km1,ap(ikp1+1),1,work,1)
                  ij = 0
                  do 190 j = 1, km1
                     jkp1 = ikp1 + j
                     ap(jkp1) = sdot(j,ap(ij+1),1,work,1)
                     call saxpy(j-1,work(j),ap(ij+1),1,ap(ikp1+1),1)
                     ij = ij + j
  190             continue
                  ap(kkp1+1) = ap(kkp1+1)
     1                         + sdot(km1,work,1,ap(ikp1+1),1)
                  ap(kkp1) = ap(kkp1)
     1                       + sdot(km1,ap(ik+1),1,ap(ikp1+1),1)
                  call scopy(km1,ap(ik+1),1,work,1)
                  ij = 0
                  do 200 j = 1, km1
                     jk = ik + j
                     ap(jk) = sdot(j,ap(ij+1),1,work,1)
                     call saxpy(j-1,work(j),ap(ij+1),1,ap(ik+1),1)
                     ij = ij + j
  200             continue
                  ap(kk) = ap(kk) + sdot(km1,work,1,ap(ik+1),1)
  210          continue
               kstep = 2
  220       continue
c
c           swap
c
            ks = abs(kpvt(k))
            if (ks .eq. k) go to 250
               iks = (ks*(ks - 1))/2
               call sswap(ks,ap(iks+1),1,ap(ik+1),1)
               ksj = ik + ks
               do 230 jb = ks, k
                  j = k + ks - jb
                  jk = ik + j
                  temp = ap(jk)
                  ap(jk) = ap(ksj)
                  ap(ksj) = temp
                  ksj = ksj - (j - 1)
  230          continue
               if (kstep .eq. 1) go to 240
                  kskp1 = ikp1 + ks
                  temp = ap(kskp1)
                  ap(kskp1) = ap(kkp1)
                  ap(kkp1) = temp
  240          continue
  250       continue
            ik = ik + k
            if (kstep .eq. 2) ik = ik + k + 1
            k = k + kstep
         go to 150
  260    continue
  270 continue
      return
      end
