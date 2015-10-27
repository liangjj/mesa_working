*deck lssuds
      subroutine lssuds (a, x, b, n, m, nrda, u, nrdu, iflag, mlso,
     +   irank, iscale, q, diag, kpivot, s, div, td, isflg, scales)
c***begin prologue  lssuds
c***subsidiary
c***purpose  subsidiary to bvsup
c***library   slatec
c***type      single precision (lssuds-s, dlssud-d)
c***author  watts, h. a., (snla)
c***description
c
c    lssuds solves the underdetermined system of equations  a z = b,
c    where a is n by m and n .le. m.  in particular, if rank a equals
c    ira, a vector x and a matrix u are determined such that x is the
c    unique solution of smallest length, satisfying a x = b, and the
c    columns of u form an orthonormal basis for the null space of a,
c    satisfying a u = 0 .  then all solutions z are given by
c              z = x + c(1)*u(1) + ..... + c(m-ira)*u(m-ira)
c    where u(j) represents the j-th column of u and the c(j) are
c    arbitrary constants.
c    if the system of equations are not compatible, only the least
c    squares solution of minimal length is computed.
c
c *********************************************************************
c   input
c *********************************************************************
c
c     a -- contains the matrix of n equations in m unknowns, a remains
c          unchanged, must be dimensioned nrda by m.
c     x -- solution array of length at least m.
c     b -- given constant vector of length n, b remains unchanged.
c     n -- number of equations, n greater or equal to 1.
c     m -- number of unknowns, m greater or equal to n.
c     nrda -- row dimension of a, nrda greater or equal to n.
c     u -- matrix used for solution, must be dimensioned nrdu by
c          (m - rank of a).
c          (storage for u may be ignored when only the minimal length
c           solution x is desired)
c     nrdu -- row dimension of u, nrdu greater or equal to m.
c             (if only the minimal length solution is wanted,
c              nrdu=0 is acceptable)
c     iflag -- status indicator
c           =0  for the first call (and for each new problem defined by
c               a new matrix a) when the matrix data is treated as exact
c           =-k for the first call (and for each new problem defined by
c               a new matrix a) when the matrix data is assumed to be
c               accurate to about k digits.
c           =1  for subsequent calls whenever the matrix a has already
c               been decomposed (problems with new vectors b but
c               same matrix a can be handled efficiently).
c     mlso -- =0 if only the minimal length solution is wanted.
c             =1 if the complete solution is wanted, includes the
c                linear space defined by the matrix u.
c     irank -- variable used for the rank of a, set by the code.
c     iscale -- scaling indicator
c               =-1 if the matrix a is to be pre-scaled by
c               columns when appropriate.
c               if the scaling indicator is not equal to -1
c               no scaling will be attempted.
c            for most problems scaling will probably not be necessary.
c     q -- matrix used for the transformation, must be dimensioned
c            nrda by m.
c     diag,kpivot,s, -- arrays of length at least n used for internal
c      div,td,scales    storage (except for scales which is m).
c     isflg -- storage for an internal variable.
c
c *********************************************************************
c   output
c *********************************************************************
c
c     iflag -- status indicator
c            =1 if solution was obtained.
c            =2 if improper input is detected.
c            =3 if rank of matrix is less than n.
c               to continue, simply reset iflag=1 and call lssuds again.
c            =4 if the system of equations appears to be inconsistent.
c               however, the least squares solution of minimal length
c               was obtained.
c     x -- minimal length least squares solution of a z = b
c     irank -- numerically determined rank of a, must not be altered
c              on succeeding calls with input values of iflag=1.
c     u -- matrix whose m-irank columns are mutually orthogonal unit
c          vectors which span the null space of a. this is to be ignored
c          when mlso was set to zero or iflag=4 on output.
c     q -- contains the strictly upper triangular part of the reduced
c           matrix and transformation information.
c     diag -- contains the diagonal elements of the triangular reduced
c             matrix.
c     kpivot -- contains the pivotal information.  the row interchanges
c               performed on the original matrix are recorded here.
c     s -- contains the solution of the lower triangular system.
c     div,td -- contains transformation information for rank
c               deficient problems.
c     scales -- contains the column scaling parameters.
c
c *********************************************************************
c
c***see also  bvsup
c***references  h. a. watts, solving linear least squares problems
c                 using sods/suds/cods, sandia report sand77-0683,
c                 sandia laboratories, 1977.
c***routines called  j4save, ohtrol, orthor, r1mach, sdot, xermax,
c                    xermsg, xgetf, xsetf
c***revision history  (yymmdd)
c   750601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900328  added type section.  (wrb)
c   900510  fixed an error message.  (rwc)
c   910408  updated the author and references sections.  (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  lssuds
      dimension a(nrda,*),x(*),b(*),u(nrdu,*),q(nrda,*),
     1          diag(*),kpivot(*),s(*),div(*),td(*),scales(*)
c
c **********************************************************************
c
c     machine precision (computer unit roundoff value) is defined
c     by the function r1mach.
c
c***first executable statement  lssuds
      uro = r1mach(4)
c
c **********************************************************************
c
      if (n .lt. 1  .or.  m .lt. n  .or.  nrda .lt. n) go to 1
      if (nrdu .ne. 0  .and.  nrdu .lt. m) go to 1
      if (iflag .le. 0) go to 5
      if (iflag .eq. 1) go to 25
c
c     invalid input for lssuds
    1 iflag=2
      call xermsg ('slatec', 'lssuds', 'invalid input parameters.', 2,
     +   1)
      return
c
    5 call xgetf(nfatal)
      maxmes = j4save (4,0,.false.)
      isflg=-15
      if (iflag .eq. 0) go to 7
      isflg=iflag
      nfat = -1
      if (nfatal .eq. 0) nfat=0
      call xsetf(nfat)
      call xermax(1)
c
c     copy matrix a into matrix q
c
    7 do 10 k=1,m
         do 10 j=1,n
   10       q(j,k)=a(j,k)
c
c     use orthogonal transformations to reduce q to lower
c     triangular form
c
      call orthor(q,n,m,nrda,iflag,irank,iscale,diag,kpivot,scales,
     1            div,td)
c
      call xsetf(nfatal)
      call xermax(maxmes)
      if (irank .eq. n) go to 15
c
c     for rank deficient problems use additional orthogonal
c     transformations to further reduce q
c
      if (irank .ne. 0) call ohtrol(q,n,nrda,diag,irank,div,td)
      return
c
c     store divisors for the triangular solution
c
   15 do 20 k=1,n
   20    div(k)=diag(k)
c
c
   25 if (irank .gt. 0) go to 40
c
c     special case for the null matrix
      do 35 k=1,m
         x(k)=0.
         if (mlso .eq. 0) go to 35
         u(k,k)=1.
         do 30 j=1,m
            if (j .eq. k) go to 30
            u(j,k)=0.
   30    continue
   35 continue
      do 37 k=1,n
         if (b(k) .gt. 0.) iflag=4
   37 continue
      return
c
c     copy constant vector into s after first interchanging
c     the elements according to the pivotal sequence
c
   40 do 45 k=1,n
         kp=kpivot(k)
   45    x(k)=b(kp)
      do 50 k=1,n
   50    s(k)=x(k)
c
      irp=irank+1
      nu=1
      if (mlso .eq. 0) nu=0
      if (irank .eq. n) go to 60
c
c     for rank deficient problems we must apply the
c     orthogonal transformation to s
c     we also check to see if the system appears to be inconsistent
c
      nmir=n-irank
      ss=sdot(n,s(1),1,s(1),1)
      do 55 l=1,irank
         k=irp-l
         gam=((td(k)*s(k))+sdot(nmir,q(irp,k),1,s(irp),1))/
     1             (td(k)*div(k))
         s(k)=s(k)+gam*td(k)
         do 55 j=irp,n
   55       s(j)=s(j)+gam*q(j,k)
      res=sdot(nmir,s(irp),1,s(irp),1)
      if (res .le. ss*(10.*max(10.**isflg,10.*uro))**2) go to 60
c
c     inconsistent system
      iflag=4
      nu=0
c
c     apply forward substitution to solve lower triangular system
c
   60 s(1)=s(1)/div(1)
      if (irank .eq. 1) go to 70
      do 65 k=2,irank
   65    s(k)=(s(k)-sdot(k-1,q(k,1),nrda,s(1),1))/div(k)
c
c     initialize x vector and then apply orthogonal transformation
c
   70 do 75 k=1,m
         x(k)=0.
         if (k .le. irank) x(k)=s(k)
   75 continue
c
      do 80 jr=1,irank
         j=irp-jr
         mj=m-j+1
         gamma=sdot(mj,q(j,j),nrda,x(j),1)/(diag(j)*q(j,j))
         do 80 k=j,m
   80       x(k)=x(k)+gamma*q(j,k)
c
c     rescale answers as dictated
c
      do 85 k=1,m
   85    x(k)=x(k)*scales(k)
c
      if ((nu .eq. 0) .or. (m .eq. irank)) return
c
c     initialize u matrix and then apply orthogonal transformation
c
      l=m-irank
      do 100 k=1,l
         do 90 i=1,m
            u(i,k)=0.
            if (i .eq. irank+k) u(i,k)=1.
   90    continue
c
         do 100 jr=1,irank
            j=irp-jr
            mj=m-j+1
            gamma=sdot(mj,q(j,j),nrda,u(j,k),1)/(diag(j)*q(j,j))
            do 100 i=j,m
  100          u(i,k)=u(i,k)+gamma*q(j,i)
c
      return
      end
