*deck scov
      subroutine scov (fcn, iopt, m, n, x, fvec, r, ldr, info, wa1, wa2,
     +   wa3, wa4)
c***begin prologue  scov
c***purpose  calculate the covariance matrix for a nonlinear data
c            fitting problem.  it is intended to be used after a
c            successful return from either snls1 or snls1e.
c***library   slatec
c***category  k1b1
c***type      single precision (scov-s, dcov-d)
c***keywords  covariance matrix, nonlinear data fitting,
c             nonlinear least squares
c***author  hiebert, k. l., (snla)
c***description
c
c  1. purpose.
c
c     scov calculates the covariance matrix for a nonlinear data
c     fitting problem.  it is intended to be used after a
c     successful return from either snls1 or snls1e. scov
c     and snls1 (and snls1e) have compatible parameters.  the
c     required external subroutine, fcn, is the same
c     for all three codes, scov, snls1, and snls1e.
c
c  2. subroutine and type statements.
c
c     subroutine scov(fcn,iopt,m,n,x,fvec,r,ldr,info,
c                     wa1,wa2,wa3,wa4)
c     integer iopt,m,n,ldr,info
c     real x(n),fvec(m),r(ldr,n),wa1(n),wa2(n),wa3(n),wa4(m)
c     external fcn
c
c  3. parameters.
c
c       fcn is the name of the user-supplied subroutine which calculates
c         the functions.  if the user wants to supply the jacobian
c         (iopt=2 or 3), then fcn must be written to calculate the
c         jacobian, as well as the functions.  see the explanation
c         of the iopt argument below.  fcn must be declared in an
c         external statement in the calling program and should be
c         written as follows.
c
c         subroutine fcn(iflag,m,n,x,fvec,fjac,ldfjac)
c         integer iflag,ldfjac,m,n
c         real x(n),fvec(m)
c         ----------
c         fjac and ldfjac may be ignored     , if iopt=1.
c         real fjac(ldfjac,n)                , if iopt=2.
c         real fjac(n)                       , if iopt=3.
c         ----------
c           iflag will never be zero when fcn is called by scov.
c         return
c         ----------
c           if iflag=1, calculate the functions at x and return
c           this vector in fvec.
c         return
c         ----------
c           if iflag=2, calculate the full jacobian at x and return
c           this matrix in fjac.  note that iflag will never be 2 unless
c           iopt=2.  fvec contains the function values at x and must
c           not be altered.  fjac(i,j) must be set to the derivative
c           of fvec(i) with respect to x(j).
c         return
c         ----------
c           if iflag=3, calculate the ldfjac-th row of the jacobian
c           and return this vector in fjac.  note that iflag will
c           never be 3 unless iopt=3.  fjac(j) must be set to
c           the derivative of fvec(ldfjac) with respect to x(j).
c         return
c         ----------
c         end
c
c
c         the value of iflag should not be changed by fcn unless the
c         user wants to terminate execution of scov.  in this case, set
c         iflag to a negative integer.
c
c
c    iopt is an input variable which specifies how the jacobian will
c         be calculated.  if iopt=2 or 3, then the user must supply the
c         jacobian, as well as the function values, through the
c         subroutine fcn.  if iopt=2, the user supplies the full
c         jacobian with one call to fcn.  if iopt=3, the user supplies
c         one row of the jacobian with each call.  (in this manner,
c         storage can be saved because the full jacobian is not stored.)
c         if iopt=1, the code will approximate the jacobian by forward
c         differencing.
c
c       m is a positive integer input variable set to the number of
c         functions.
c
c       n is a positive integer input variable set to the number of
c         variables.  n must not exceed m.
c
c       x is an array of length n.  on input x must contain the value
c         at which the covariance matrix is to be evaluated.  this is
c         usually the value for x returned from a successful run of
c         snls1 (or snls1e).  the value of x will not be changed.
c
c    fvec is an output array of length m which contains the functions
c         evaluated at x.
c
c       r is an output array.  for iopt=1 and 2, r is an m by n array.
c         for iopt=3, r is an n by n array.  on output, if info=1,
c         the upper n by n submatrix of r contains the covariance
c         matrix evaluated at x.
c
c     ldr is a positive integer input variable which specifies
c         the leading dimension of the array r.  for iopt=1 and 2,
c         ldr must not be less than m.  for iopt=3, ldr must not
c         be less than n.
c
c    info is an integer output variable.  if the user has terminated
c         execution, info is set to the (negative) value of iflag.  see
c         description of fcn. otherwise, info is set as follows.
c
c         info = 0 improper input parameters (m.le.0 or n.le.0).
c
c         info = 1 successful return.  the covariance matrix has been
c                  calculated and stored in the upper n by n
c                  submatrix of r.
c
c         info = 2 the jacobian matrix is singular for the input value
c                  of x.  the covariance matrix cannot be calculated.
c                  the upper n by n submatrix of r contains the qr
c                  factorization of the jacobian (probably not of
c                  interest to the user).
c
c     wa1 is a work array of length n.
c     wa2 is a work array of length n.
c     wa3 is a work array of length n.
c     wa4 is a work array of length m.
c
c***references  (none)
c***routines called  enorm, fdjac3, qrfac, rwupdt, xermsg
c***revision history  (yymmdd)
c   810522  date written
c   890505  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900510  fixed an error message.  (rwc)
c***end prologue  scov
c
c     revised 820707-1100
c     revised yymmdd hhmm
c
      integer i,idum,iflag,info,iopt,j,k,kp1,ldr,m,n,nm1,nrow
      real x(*),r(ldr,*),fvec(*),wa1(*),wa2(*),wa3(*),wa4(*)
      external fcn
      real one,sigma,temp,zero
      logical sing
      save zero, one
      data zero/0.e0/,one/1.e0/
c***first executable statement  scov
      sing=.false.
      iflag=0
      if (m.le.0 .or. n.le.0) go to 300
c
c     calculate sigma = (sum of the squared residuals) / (m-n)
      iflag=1
      call fcn(iflag,m,n,x,fvec,r,ldr)
      if (iflag.lt.0) go to 300
      temp=enorm(m,fvec)
      sigma=one
      if (m.ne.n) sigma=temp*temp/(m-n)
c
c     calculate the jacobian
      if (iopt.eq.3) go to 200
c
c     store the full jacobian using m*n storage
      if (iopt.eq.1) go to 100
c
c     user supplies the jacobian
      iflag=2
      call fcn(iflag,m,n,x,fvec,r,ldr)
      go to 110
c
c     code approximates the jacobian
100   call fdjac3(fcn,m,n,x,fvec,r,ldr,iflag,zero,wa4)
110   if (iflag.lt.0) go to 300
c
c     compute the qr decomposition
      call qrfac(m,n,r,ldr,.false.,idum,1,wa1,wa1,wa1)
      do 120 i=1,n
120   r(i,i)=wa1(i)
      go to 225
c
c     compute the qr factorization of the jacobian matrix calculated one
c     row at a time and stored in the upper triangle of r.
c     ( (q transpose)*fvec is also calculated but not used.)
200   continue
      do 210 j=1,n
      wa2(j)=zero
      do 205 i=1,n
      r(i,j)=zero
205   continue
210   continue
      iflag=3
      do 220 i=1,m
      nrow = i
      call fcn(iflag,m,n,x,fvec,wa1,nrow)
      if (iflag.lt.0) go to 300
      temp=fvec(i)
      call rwupdt(n,r,ldr,wa1,wa2,temp,wa3,wa4)
220   continue
c
c     check if r is singular.
225   continue
      do 230 i=1,n
      if (r(i,i).eq.zero) sing=.true.
230   continue
      if (sing) go to 300
c
c     r is upper triangular.  calculate (r transpose) inverse and store
c     in the upper triangle of r.
      if (n.eq.1) go to 275
      nm1=n-1
      do 270 k=1,nm1
c
c     initialize the right-hand side (wa1(*)) as the k-th column of the
c     identity matrix.
      do 240 i=1,n
      wa1(i)=zero
240   continue
      wa1(k)=one
c
      r(k,k)=wa1(k)/r(k,k)
      kp1=k+1
      do 260 i=kp1,n
c
c     subtract r(k,i-1)*r(i-1,*) from the right-hand side, wa1(*).
      do 250 j=i,n
      wa1(j)=wa1(j)-r(k,i-1)*r(i-1,j)
250   continue
      r(k,i)=wa1(i)/r(i,i)
260   continue
270   continue
275   r(n,n)=one/r(n,n)
c
c     calculate r-inverse * (r transpose) inverse and store in the upper
c     triangle of r.
      do 290 i=1,n
      do 290 j=i,n
      temp=zero
      do 280 k=j,n
      temp=temp+r(i,k)*r(j,k)
280   continue
      r(i,j)=temp*sigma
290   continue
      info=1
c
300   continue
      if (m.le.0 .or. n.le.0) info=0
      if (iflag.lt.0) info=iflag
      if (sing) info=2
      if (info .lt. 0) call xermsg ('slatec', 'scov',
     +   'execution terminated because user set iflag negative.', 1, 1)
      if (info .eq. 0) call xermsg ('slatec', 'scov',
     +   'invalid input parameter.', 2, 1)
      if (info .eq. 2) call xermsg ('slatec', 'scov',
     +   'singular jacobian matrix, covariance matrix cannot be ' //
     +   'calculated.', 1, 1)
      return
      end
