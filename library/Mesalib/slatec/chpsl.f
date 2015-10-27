*deck chpsl
      subroutine chpsl (ap, n, kpvt, b)
c***begin prologue  chpsl
c***purpose  solve a complex hermitian system using factors obtained
c            from chpfa.
c***library   slatec (linpack)
c***category  d2d1a
c***type      complex (sspsl-s, dspsl-d, chpsl-c, cspsl-c)
c***keywords  hermitian, linear algebra, linpack, matrix, packed, solve
c***author  bunch, j., (ucsd)
c***description
c
c     chisl solves the complex hermitian system
c     a * x = b
c     using the factors computed by chpfa.
c
c     on entry
c
c        ap      complex(n*(n+1)/2)
c                the output from chpfa.
c
c        n       integer
c                the order of the matrix  a .
c
c        kvpt    integer(n)
c                the pivot vector from chpfa.
c
c        b       complex(n)
c                the right hand side vector.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero may occur if  chpco  has set rcond .eq. 0.0
c        or  chpfa  has set info .ne. 0  .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call chpfa(ap,n,kvpt,info)
c           if (info .ne. 0) go to ...
c           do 10 j = 1, p
c              call chpsl(ap,n,kvpt,c(1,j))
c        10 continue
c
c***references  j. j. dongarra, j. r. bunch, c. b. moler, and g. w.
c                 stewart, linpack users' guide, siam, 1979.
c***routines called  caxpy, cdotc
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
c***end prologue  chpsl
      integer n,kpvt(*)
      complex ap(*),b(*)
c
      complex ak,akm1,bk,bkm1,cdotc,denom,temp
      integer ik,ikm1,ikp1,k,kk,km1k,km1km1,kp
c
c     loop backward applying the transformations and
c     d inverse to b.
c
c***first executable statement  chpsl
      k = n
      ik = (n*(n - 1))/2
   10 if (k .eq. 0) go to 80
         kk = ik + k
         if (kpvt(k) .lt. 0) go to 40
c
c           1 x 1 pivot block.
c
            if (k .eq. 1) go to 30
               kp = kpvt(k)
               if (kp .eq. k) go to 20
c
c                 interchange.
c
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
   20          continue
c
c              apply the transformation.
c
               call caxpy(k-1,b(k),ap(ik+1),1,b(1),1)
   30       continue
c
c           apply d inverse.
c
            b(k) = b(k)/ap(kk)
            k = k - 1
            ik = ik - k
         go to 70
   40    continue
c
c           2 x 2 pivot block.
c
            ikm1 = ik - (k - 1)
            if (k .eq. 2) go to 60
               kp = abs(kpvt(k))
               if (kp .eq. k - 1) go to 50
c
c                 interchange.
c
                  temp = b(k-1)
                  b(k-1) = b(kp)
                  b(kp) = temp
   50          continue
c
c              apply the transformation.
c
               call caxpy(k-2,b(k),ap(ik+1),1,b(1),1)
               call caxpy(k-2,b(k-1),ap(ikm1+1),1,b(1),1)
   60       continue
c
c           apply d inverse.
c
            km1k = ik + k - 1
            kk = ik + k
            ak = ap(kk)/conjg(ap(km1k))
            km1km1 = ikm1 + k - 1
            akm1 = ap(km1km1)/ap(km1k)
            bk = b(k)/conjg(ap(km1k))
            bkm1 = b(k-1)/ap(km1k)
            denom = ak*akm1 - 1.0e0
            b(k) = (akm1*bk - bkm1)/denom
            b(k-1) = (ak*bkm1 - bk)/denom
            k = k - 2
            ik = ik - (k + 1) - k
   70    continue
      go to 10
   80 continue
c
c     loop forward applying the transformations.
c
      k = 1
      ik = 0
   90 if (k .gt. n) go to 160
         if (kpvt(k) .lt. 0) go to 120
c
c           1 x 1 pivot block.
c
            if (k .eq. 1) go to 110
c
c              apply the transformation.
c
               b(k) = b(k) + cdotc(k-1,ap(ik+1),1,b(1),1)
               kp = kpvt(k)
               if (kp .eq. k) go to 100
c
c                 interchange.
c
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
  100          continue
  110       continue
            ik = ik + k
            k = k + 1
         go to 150
  120    continue
c
c           2 x 2 pivot block.
c
            if (k .eq. 1) go to 140
c
c              apply the transformation.
c
               b(k) = b(k) + cdotc(k-1,ap(ik+1),1,b(1),1)
               ikp1 = ik + k
               b(k+1) = b(k+1) + cdotc(k-1,ap(ikp1+1),1,b(1),1)
               kp = abs(kpvt(k))
               if (kp .eq. k) go to 130
c
c                 interchange.
c
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
  130          continue
  140       continue
            ik = ik + k + k + 1
            k = k + 2
  150    continue
      go to 90
  160 continue
      return
      end