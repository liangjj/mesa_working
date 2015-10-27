*deck orthol
      subroutine orthol (a, m, n, nrda, iflag, irank, iscale, diag,
     +   kpivot, scales, cols, cs)
c***begin prologue  orthol
c***subsidiary
c***purpose  subsidiary to bvsup
c***library   slatec
c***type      single precision (orthol-s)
c***author  watts, h. a., (snla)
c***description
c
c   reduction of the matrix a to upper triangular form by a sequence of
c   orthogonal householder transformations pre-multiplying a
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
c     m -- number of rows in the matrix, m greater or equal to n
c     n -- number of columns in the matrix, n greater or equal to 1
c     iflag -- indicates the uncertainty in the matrix data
c             = 0 when the data is to be treated as exact
c             =-k when the data is assumed to be accurate to about
c                 k digits
c     iscale -- scaling indicator
c               =-1 if the matrix a is to be pre-scaled by
c               columns when appropriate.
c               otherwise no scaling will be attempted
c     nrda -- row dimension of a, nrda greater or equal to m
c     diag,kpivot,cols -- arrays of length at least n used internally
c         ,cs,scales
c
c **********************************************************************
c   output
c **********************************************************************
c
c     iflag - status indicator
c            =1 for successful decomposition
c            =2 if improper input is detected
c            =3 if rank of the matrix is less than n
c     a -- contains the reduced matrix in the strictly upper triangular
c          part and transformation information in the lower part
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
c   900402  added type section.  (wrb)
c   910408  updated the author and references sections.  (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  orthol
      dimension a(nrda,*),diag(*),kpivot(*),cols(*),cs(*),scales(*)
c
c **********************************************************************
c
c     machine precision (computer unit roundoff value) is defined
c     by the function r1mach.
c
c***first executable statement  orthol
      uro = r1mach(3)
c
c **********************************************************************
c
      if (m .ge. n  .and.  n .ge. 1  .and.  nrda .ge. m) go to 1
      iflag=2
      call xermsg ('slatec', 'orthol', 'invalid input parameters.', 2,
     +   1)
      return
c
    1 acc=10.*uro
      if (iflag .lt. 0) acc=max(acc,10.**iflag)
      sruro=sqrt(uro)
      iflag=1
      irank=n
c
c     compute norm**2 of jth column and a matrix norm
c
      anorm=0.
      do 2 j=1,n
         kpivot(j)=j
         cols(j)=sdot(m,a(1,j),1,a(1,j),1)
         cs(j)=cols(j)
         anorm=anorm+cols(j)
    2 continue
c
c     perform column scaling on a when specified
c
      call cscale(a,nrda,m,n,cols,cs,dum,dum,anorm,scales,iscale,0)
c
      anorm=sqrt(anorm)
c
c
c     construction of upper triangular matrix and recording of
c     orthogonal transformations
c
c
      do 50 k=1,n
         mk=m-k+1
         if (k .eq. n) go to 25
         kp=k+1
c
c        searching for pivotal column
c
         do 10 j=k,n
            if (cols(j) .ge. sruro*cs(j)) go to 5
            cols(j)=sdot(mk,a(k,j),1,a(k,j),1)
            cs(j)=cols(j)
    5       if (j .eq. k) go to 7
            if (sigma .ge. 0.99*cols(j)) go to 10
    7       sigma=cols(j)
            jcol=j
   10    continue
         if (jcol .eq. k) go to 25
c
c        perform column interchange
c
         l=kpivot(k)
         kpivot(k)=kpivot(jcol)
         kpivot(jcol)=l
         cols(jcol)=cols(k)
         cols(k)=sigma
         css=cs(k)
         cs(k)=cs(jcol)
         cs(jcol)=css
         sc=scales(k)
         scales(k)=scales(jcol)
         scales(jcol)=sc
         do 20 l=1,m
            asave=a(l,k)
            a(l,k)=a(l,jcol)
   20       a(l,jcol)=asave
c
c        check rank of the matrix
c
   25    sig=sdot(mk,a(k,k),1,a(k,k),1)
         diagk=sqrt(sig)
         if (diagk .gt. acc*anorm) go to 30
c
c        rank deficient problem
         iflag=3
         irank=k-1
         call xermsg ('slatec', 'orthol',
     +      'rank of matrix is less than the number of columns.', 1, 1)
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
            as=sdot(mk,a(k,k),1,a(k,j),1)/sad
            do 35 l=k,m
   35          a(l,j)=a(l,j)+as*a(l,k)
   40       cols(j)=cols(j)-a(k,j)**2
   50 continue
c
c
      return
      end
