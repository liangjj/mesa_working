*deck orthor
      subroutine orthor (a, n, m, nrda, iflag, irank, iscale, diag,
     +   kpivot, scales, rows, rs)
c***begin prologue  orthor
c***subsidiary
c***purpose  subsidiary to bvsup
c***library   slatec
c***type      single precision (orthor-s, dorthr-d)
c***author  watts, h. a., (snla)
c***description
c
c   reduction of the matrix a to lower triangular form by a sequence of
c   orthogonal householder transformations post-multiplying a
c
c   modeled after the algol codes in the articles in the references
c   section.
c
c **********************************************************************
c   input
c **********************************************************************
c
c     a -- contains the matrix to be decomposed, must be dimensioned
c           nrda by n
c     n -- number of rows in the matrix, n greater or equal to 1
c     m -- number of columns in the matrix, m greater or equal to n
c     iflag -- indicates the uncertainty in the matrix data
c             = 0 when the data is to be treated as exact
c             =-k when the data is assumed to be accurate to about
c                 k digits
c     iscale -- scaling indicator
c               =-1 if the matrix is to be pre-scaled by
c               columns when appropriate.
c               otherwise no scaling will be attempted
c     nrda -- row dimension of a, nrda greater or equal to n
c     diag,kpivot,rows -- arrays of length at least n used internally
c         ,rs,scales         (except for scales which is m)
c
c **********************************************************************
c   output
c **********************************************************************
c
c     iflag - status indicator
c            =1 for successful decomposition
c            =2 if improper input is detected
c            =3 if rank of the matrix is less than n
c     a -- contains the reduced matrix in the strictly lower triangular
c          part and transformation information
c     irank -- contains the numerically determined matrix rank
c     diag -- contains the diagonal elements of the reduced
c             triangular matrix
c     kpivot -- contains the pivotal information, the column
c               interchanges performed on the original matrix are
c               recorded here.
c     scales -- contains the column scaling parameters
c
c **********************************************************************
c
c***see also  bvsup
c***references  g. golub, numerical methods for solving linear least
c                 squares problems, numerische mathematik 7, (1965),
c                 pp. 206-216.
c               p. businger and g. golub, linear least squares
c                 solutions by householder transformations, numerische
c                 mathematik  7, (1965), pp. 269-276.
c***routines called  cscale, r1mach, sdot, xermsg
c***revision history  (yymmdd)
c   750601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900328  added type section.  (wrb)
c   910408  updated the author and references sections.  (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  orthor
      dimension a(nrda,*),diag(*),kpivot(*),rows(*),rs(*),scales(*)
c
c end of abstract
c
c **********************************************************************
c
c     machine precision (computer unit roundoff value) is defined
c     by the function r1mach.
c
c **********************************************************************
c
c***first executable statement  orthor
      uro = r1mach(4)
      if (m .ge. n  .and.  n .ge. 1  .and.  nrda .ge. n) go to 1
      iflag=2
      call xermsg ('slatec', 'orthor', 'invalid input parameters.', 2,
     +   1)
      return
c
    1 acc=10.*uro
      if (iflag .lt. 0) acc=max(acc,10.**iflag)
      sruro=sqrt(uro)
      iflag=1
      irank=n
c
c     compute norm**2 of jth row and a matrix norm
c
      anorm=0.
      do 2 j=1,n
         kpivot(j)=j
         rows(j)=sdot(m,a(j,1),nrda,a(j,1),nrda)
         rs(j)=rows(j)
         anorm=anorm+rows(j)
    2 continue
c
c     perform column scaling on a when specified
c
      call cscale(a,nrda,n,m,scales,dum,rows,rs,anorm,scales,iscale,1)
c
      anorm=sqrt(anorm)
c
c
c     construction of lower triangular matrix and recording of
c     orthogonal transformations
c
c
      do 50 k=1,n
         mk=m-k+1
         if (k .eq. n) go to 25
         kp=k+1
c
c        searching for pivotal row
c
         do 10 j=k,n
            if (rows(j) .ge. sruro*rs(j)) go to 5
            rows(j)=sdot(mk,a(j,k),nrda,a(j,k),nrda)
            rs(j)=rows(j)
    5       if (j .eq. k) go to 7
            if (sigma .ge. 0.99*rows(j)) go to 10
    7       sigma=rows(j)
            jrow=j
   10    continue
         if (jrow .eq. k) go to 25
c
c        perform row interchange
c
         l=kpivot(k)
         kpivot(k)=kpivot(jrow)
         kpivot(jrow)=l
         rows(jrow)=rows(k)
         rows(k)=sigma
         rss=rs(k)
         rs(k)=rs(jrow)
         rs(jrow)=rss
         do 20 l=1,m
            asave=a(k,l)
            a(k,l)=a(jrow,l)
   20       a(jrow,l)=asave
c
c        check rank of the matrix
c
   25    sig=sdot(mk,a(k,k),nrda,a(k,k),nrda)
         diagk=sqrt(sig)
         if (diagk .gt. acc*anorm) go to 30
c
c        rank deficient problem
         iflag=3
         irank=k-1
         call xermsg ('slatec', 'orthor',
     +      'rank of matrix is less than the number of rows.', 1, 1)
         return
c
c        construct and apply transformation to matrix a
c
   30    akk=a(k,k)
         if (akk .gt. 0.) diagk=-diagk
         diag(k)=diagk
         a(k,k)=akk-diagk
         if (k .eq. n) go to 50
         sad=diagk*akk-sig
         do 40 j=kp,n
            as=sdot(mk,a(k,k),nrda,a(j,k),nrda)/sad
            do 35 l=k,m
   35          a(j,l)=a(j,l)+as*a(k,l)
   40       rows(j)=rows(j)-a(j,k)**2
   50 continue
c
c
      return
      end
