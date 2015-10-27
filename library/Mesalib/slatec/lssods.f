*deck lssods
      subroutine lssods (a, x, b, m, n, nrda, iflag, irank, iscale, q,
     +   diag, kpivot, iter, resnrm, xnorm, z, r, div, td, scales)
c***begin prologue  lssods
c***subsidiary
c***purpose  subsidiary to bvsup
c***library   slatec
c***type      single precision (lssods-s)
c***author  (unknown)
c***description
c
c     lssods solves the same problem as sods (in fact, it is called by
c     sods) but is somewhat more flexible in its use. in particular,
c     lssods allows for iterative refinement of the solution, makes the
c     transformation and triangular reduction information more
c     accessible, and enables the user to avoid destruction of the
c     original matrix a.
c
c     modeled after the algol codes in the articles in the references
c     section.
c
c **********************************************************************
c   input
c **********************************************************************
c
c     a -- contains the matrix of m equations in n unknowns and must
c          be dimensioned nrda by n. a remains unchanged
c     x -- solution array of length at least n
c     b -- given constant vector of length m, b remains unchanged
c     m -- number of equations, m greater or equal to 1
c     n -- number of unknowns, n not larger than m
c  nrda -- row dimension of a, nrda greater or equal to m
c iflag -- status indicator
c         = 0 for the first call (and for each new problem defined by
c             a new matrix a) when the matrix data is treated as exact
c         =-k for the first call (and for each new problem defined by
c             a new matrix a) when the matrix data is assumed to be
c             accurate to about k digits
c         = 1 for subsequent calls whenever the matrix a has already
c             been decomposed (problems with new vectors b but
c             same matrix a can be handled efficiently)
c iscale -- scaling indicator
c         =-1 if the matrix a is to be pre-scaled by
c             columns when appropriate
c             if the scaling indicator is not equal to -1
c             no scaling will be attempted
c             for most problems scaling will probably not be necessary
c   iter -- maximum number of iterative improvement steps to be
c           performed,  0 .le. iter .le. 10   (sods uses iter=0)
c      q -- matrix used for the transformation, must be dimensioned
c           nrda by n  (sods puts a in the q location which conserves
c           storage but destroys a)
c           when iterative improvement of the solution is requested,
c           iter .gt. 0, this additional storage for q must be
c           made available
c diag,kpivot,z,r, -- arrays of length n (except for r which is m)
c   div,td,scales     used for internal storage
c
c **********************************************************************
c   output
c **********************************************************************
c
c  iflag -- status indicator
c            =1 if solution was obtained
c            =2 if improper input is detected
c            =3 if rank of matrix is less than n
c               if the minimal length least squares solution is
c               desired, simply reset iflag=1 and call the code again
c
c       the next three iflag values can occur only when
c        the iterative improvement mode is being used.
c            =4 if the problem is ill-conditioned and maximal
c               machine accuracy is not achievable
c            =5 if the problem is very ill-conditioned and the solution
c               is likely to have no correct digits
c            =6 if the allowable number of iterative improvement steps
c               has been completed without getting convergence
c      x -- least squares solution of  a x = b
c  irank -- contains the numerically determined matrix rank
c           the user must not alter this value on succeeding calls
c           with input values of iflag=1
c      q -- contains the strictly upper triangular part of the reduced
c           matrix and the transformation information in the lower
c           triangular part
c   diag -- contains the diagonal elements of the triangular reduced
c           matrix
c kpivot -- contains the pivotal information.  the column interchanges
c           performed on the original matrix are recorded here
c   iter -- the actual number of iterative corrections used
c resnrm -- the euclidean norm of the residual vector  b - a x
c  xnorm -- the euclidean norm of the solution vector
c div,td -- contains transformation information for rank
c           deficient problems
c scales -- contains the column scaling parameters
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
c***routines called  j4save, ohtror, orthol, r1mach, sdot, sdsdot,
c                    xermax, xermsg, xgetf, xsetf
c***revision history  (yymmdd)
c   750601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900402  added type section.  (wrb)
c   910408  updated the references section.  (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  lssods
      dimension a(nrda,*),x(*),b(*),q(nrda,*),diag(*),
     1          z(*),kpivot(*),r(*),div(*),td(*),scales(*)
c
c **********************************************************************
c
c     machine precision (computer unit roundoff value) is defined
c     the function r1mach.
c
c***first executable statement  lssods
      uro = r1mach(3)
c
c **********************************************************************
c
      if (n .lt. 1  .or.  m .lt. n  .or.  nrda .lt. m) go to 1
      if (iter .lt. 0) go to 1
      if (iflag .le. 0) go to 5
      if (iflag .eq. 1) go to 15
c
c     invalid input for lssods
    1 iflag=2
      call xermsg ('slatec', 'lssods', 'invalid input parameters.', 2,
     +   1)
      return
c
    5 call xgetf (nfatal)
      maxmes = j4save (4,0,.false.)
      if (iflag .eq. 0) go to 7
      nfat = -1
      if(nfatal .eq. 0) nfat=0
      call xsetf (nfat)
      call xermax (1)
c
c     copy matrix a into matrix q
c
    7 do 10 j=1,n
         do 10 k=1,m
   10       q(k,j)=a(k,j)
c
c     use orthogonal transformations to reduce q to
c     upper triangular form
c
      call orthol(q,m,n,nrda,iflag,irank,iscale,diag,kpivot,scales,z,td)
c
      call xsetf (nfatal)
      call xermax (maxmes)
      if (irank .eq. n) go to 12
c
c     for rank deficient problems use additional orthogonal
c     transformations to further reduce q
c
      if (irank .ne. 0) call ohtror(q,n,nrda,diag,irank,div,td)
      return
c
c     store divisors for the triangular solution
c
   12 do 13 k=1,n
   13    div(k)=diag(k)
c
   15 irm=irank-1
      irp=irank+1
      iterp=min(iter+1,11)
      acc=10.*uro
c
c     zero out solution array
c
      do 20 k=1,n
   20    x(k)=0.
c
      if (irank .gt. 0) go to 25
c
c     special case for the null matrix
      iter=0
      xnorm=0.
      resnrm=sqrt(sdot(m,b(1),1,b(1),1))
      return
c
c     copy constant vector into r
c
   25 do 30 k=1,m
   30    r(k)=b(k)
c
c **********************************************************************
c     solution section
c     iterative refinement of the residual vector
c **********************************************************************
c
      do 100 it=1,iterp
         iter=it-1
c
c        apply orthogonal transformation to r
c
         do 35 j=1,irank
            mj=m-j+1
            gamma=sdot(mj,q(j,j),1,r(j),1)/(diag(j)*q(j,j))
            do 35 k=j,m
   35          r(k)=r(k)+gamma*q(k,j)
c
c        backward substitution for triangular system solution
c
         z(irank)=r(irank)/div(irank)
         if (irm .eq. 0) go to 45
         do 40 l=1,irm
            k=irank-l
            kp=k+1
   40       z(k)=(r(k)-sdot(l,q(k,kp),nrda,z(kp),1))/div(k)
c
   45    if (irank .eq. n) go to 60
c
c        for rank deficient problems obtain the
c        minimal length solution
c
         nmir=n-irank
         do 50 k=irp,n
   50       z(k)=0.
         do 55 k=1,irank
            gam=((td(k)*z(k))+sdot(nmir,q(k,irp),nrda,z(irp),1))/
     1                (td(k)*div(k))
            z(k)=z(k)+gam*td(k)
            do 55 j=irp,n
   55          z(j)=z(j)+gam*q(k,j)
c
c        reorder solution components according to pivotal points
c        and rescale answers as dictated
c
   60    do 65 k=1,n
            z(k)=z(k)*scales(k)
            l=kpivot(k)
   65       x(l)=x(l)+z(k)
c
c        compute correction vector norm (solution norm)
c
         znorm=sqrt(sdot(n,z(1),1,z(1),1))
         if (it .eq. 1) xnorm=znorm
         if (iterp .gt. 1) go to 80
c
c        no iterative corrections to be performed, so compute
c        the approximate residual norm defined by the equations
c        which are not satisfied by the solution
c        then we are done
c
         mmir=m-irank
         if (mmir .eq. 0) go to 70
         resnrm=sqrt(sdot(mmir,r(irp),1,r(irp),1))
         return
   70    resnrm=0.
         return
c
c        compute residual vector for the iterative improvement process
c
   80    do 85 k=1,m
   85       r(k)=-sdsdot(n,-b(k),a(k,1),nrda,x(1),1)
         resnrm=sqrt(sdot(m,r(1),1,r(1),1))
         if (it .eq. 1) go to 100
c
c        test for convergence
c
         if (znorm .le. acc*xnorm) return
c
c        compare successive refinement vector norms
c        for loop termination criteria
c
         if (znorm .le. 0.25*znrm0) go to 100
         if (it .eq. 2) go to 90
c
         iflag=4
         call xermsg ('slatec', 'lssods',
     +   'problem may be ill-conditioned.  maximal machine accuracy ' //
     +   'is not achievable.', 3, 1)
         return
c
   90    iflag=5
         call xermsg ('slatec', 'lssods',
     +      'problem is very ill-conditioned.  iterative ' //
     +      'improvement is ineffective.', 8, 1)
         return
c
  100    znrm0=znorm
c **********************************************************************
c
c **********************************************************************
      iflag=6
         call xermsg ('slatec', 'lssods',
     +      'convergence has not been obtained with allowable ' //
     +      'number of iterative improvement steps.', 8, 1)
c
      return
      end
