*deck dorthr
      subroutine dorthr (a, n, m, nrda, iflag, irank, iscale, diag,
     +   kpivot, scales, rows, rs)
c***begin prologue  dorthr
c***subsidiary
c***purpose  subsidiary to dbvsup and dsuds
c***library   slatec
c***type      double precision (orthor-s, dorthr-d)
c***author  watts, h. a., (snla)
c***description
c
c   reduction of the matrix a to lower triangular form by a sequence of
c   orthogonal householder transformations post-multiplying a.
c
c *********************************************************************
c   input
c *********************************************************************
c
c     a -- contains the matrix to be decomposed, must be dimensioned
c           nrda by n.
c     n -- number of rows in the matrix, n greater or equal to 1.
c     m -- number of columns in the matrix, m greater or equal to n.
c     iflag -- indicates the uncertainty in the matrix data.
c             = 0 when the data is to be treated as exact.
c             =-k when the data is assumed to be accurate to about
c                 k digits.
c     iscale -- scaling indicator.
c               =-1 if the matrix is to be pre-scaled by
c               columns when appropriate.
c               otherwise no scaling will be attempted.
c     nrda -- row dimension of a, nrda greater or equal to n.
c     diag,kpivot,rows, -- arrays of length at least n used internally
c          rs,scales         (except for scales which is m).
c
c *********************************************************************
c   output
c *********************************************************************
c
c     iflag - status indicator
c            =1 for successful decomposition.
c            =2 if improper input is detected.
c            =3 if rank of the matrix is less than n.
c     a -- contains the reduced matrix in the strictly lower triangular
c          part and transformation information.
c     irank -- contains the numerically determined matrix rank.
c     diag -- contains the diagonal elements of the reduced
c             triangular matrix.
c     kpivot -- contains the pivotal information, the column
c               interchanges performed on the original matrix are
c               recorded here.
c     scales -- contains the column scaling parameters.
c
c *********************************************************************
c
c***see also  dbvsup, dsuds
c***references  g. golub, numerical methods for solving linear least
c                 squares problems, numerische mathematik 7, (1965),
c                 pp. 206-216.
c               p. businger and g. golub, linear least squares
c                 solutions by householder transformations, numerische
c                 mathematik  7, (1965), pp. 269-276.
c***routines called  d1mach, dcscal, ddot, xermsg
c***revision history  (yymmdd)
c   750601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900328  added type section.  (wrb)
c   910408  updated the author and references sections.  (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  dorthr
      double precision ddot, d1mach
      integer iflag, irank, iscale, j, jrow, k, kp, kpivot(*), l, m,
     1     mk, n, nrda
      double precision a(nrda,*), acc, akk, anorm, as, asave, diag(*),
     1     diagk, dum, rows(*), rs(*), rss, sad, scales(*), sig, sigma,
     2     sruro, uro
c
c     ******************************************************************
c
c          machine precision (computer unit roundoff value) is defined
c          by the function d1mach.
c
c     ******************************************************************
c
c***first executable statement  dorthr
      uro = d1mach(4)
      if (m .ge. n .and. n .ge. 1 .and. nrda .ge. n) go to 10
         iflag = 2
         call xermsg ('slatec', 'dorthr', 'invalid input parameters.',
     +      2, 1)
      go to 150
   10 continue
c
         acc = 10.0d0*uro
         if (iflag .lt. 0) acc = max(acc,10.0d0**iflag)
         sruro = sqrt(uro)
         iflag = 1
         irank = n
c
c        compute norm**2 of jth row and a matrix norm
c
         anorm = 0.0d0
         do 20 j = 1, n
            kpivot(j) = j
            rows(j) = ddot(m,a(j,1),nrda,a(j,1),nrda)
            rs(j) = rows(j)
            anorm = anorm + rows(j)
   20    continue
c
c        perform column scaling on a when specified
c
         call dcscal(a,nrda,n,m,scales,dum,rows,rs,anorm,scales,iscale,
     1               1)
c
         anorm = sqrt(anorm)
c
c
c        construction of lower triangular matrix and recording of
c        orthogonal transformations
c
c
         do 130 k = 1, n
c           begin block permitting ...exits to 80
               mk = m - k + 1
c           ...exit
               if (k .eq. n) go to 80
               kp = k + 1
c
c              searching for pivotal row
c
               do 60 j = k, n
c                 begin block permitting ...exits to 50
                     if (rows(j) .ge. sruro*rs(j)) go to 30
                        rows(j) = ddot(mk,a(j,k),nrda,a(j,k),nrda)
                        rs(j) = rows(j)
   30                continue
                     if (j .eq. k) go to 40
c                 ......exit
                        if (sigma .ge. 0.99d0*rows(j)) go to 50
   40                continue
                     sigma = rows(j)
                     jrow = j
   50             continue
   60          continue
c           ...exit
               if (jrow .eq. k) go to 80
c
c              perform row interchange
c
               l = kpivot(k)
               kpivot(k) = kpivot(jrow)
               kpivot(jrow) = l
               rows(jrow) = rows(k)
               rows(k) = sigma
               rss = rs(k)
               rs(k) = rs(jrow)
               rs(jrow) = rss
               do 70 l = 1, m
                  asave = a(k,l)
                  a(k,l) = a(jrow,l)
                  a(jrow,l) = asave
   70          continue
   80       continue
c
c           check rank of the matrix
c
            sig = ddot(mk,a(k,k),nrda,a(k,k),nrda)
            diagk = sqrt(sig)
            if (diagk .gt. acc*anorm) go to 90
c
c              rank deficient problem
               iflag = 3
               irank = k - 1
               call xermsg ('slatec', 'dorthr',
     +            'rank of matrix is less than the number of rows.', 1,
     +            1)
c        ......exit
               go to 140
   90       continue
c
c           construct and apply transformation to matrix a
c
            akk = a(k,k)
            if (akk .gt. 0.0d0) diagk = -diagk
            diag(k) = diagk
            a(k,k) = akk - diagk
            if (k .eq. n) go to 120
               sad = diagk*akk - sig
               do 110 j = kp, n
                  as = ddot(mk,a(k,k),nrda,a(j,k),nrda)/sad
                  do 100 l = k, m
                     a(j,l) = a(j,l) + as*a(k,l)
  100             continue
                  rows(j) = rows(j) - a(j,k)**2
  110          continue
  120       continue
  130    continue
  140    continue
  150 continue
c
c
      return
      end
