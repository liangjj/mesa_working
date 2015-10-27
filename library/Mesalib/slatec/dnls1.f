*deck dnls1
      subroutine dnls1 (fcn, iopt, m, n, x, fvec, fjac, ldfjac, ftol,
     +   xtol, gtol, maxfev, epsfcn, diag, mode, factor, nprint, info,
     +   nfev, njev, ipvt, qtf, wa1, wa2, wa3, wa4)
c***begin prologue  dnls1
c***purpose  minimize the sum of the squares of m nonlinear functions
c            in n variables by a modification of the levenberg-marquardt
c            algorithm.
c***library   slatec
c***category  k1b1a1, k1b1a2
c***type      double precision (snls1-s, dnls1-d)
c***keywords  levenberg-marquardt, nonlinear data fitting,
c             nonlinear least squares
c***author  hiebert, k. l., (snla)
c***description
c
c 1. purpose.
c
c       the purpose of dnls1 is to minimize the sum of the squares of m
c       nonlinear functions in n variables by a modification of the
c       levenberg-marquardt algorithm.  the user must provide a subrou-
c       tine which calculates the functions.  the user has the option
c       of how the jacobian will be supplied.  the user can supply the
c       full jacobian, or the rows of the jacobian (to avoid storing
c       the full jacobian), or let the code approximate the jacobian by
c       forward-differencing.   this code is the combination of the
c       minpack codes (argonne) lmder, lmdif, and lmstr.
c
c
c 2. subroutine and type statements.
c
c       subroutine dnls1(fcn,iopt,m,n,x,fvec,fjac,ldfjac,ftol,xtol,
c      *                 gtol,maxfev,epsfcn,diag,mode,factor,nprint,info
c      *                 ,nfev,njev,ipvt,qtf,wa1,wa2,wa3,wa4)
c       integer iopt,m,n,ldfjac,maxfev,mode,nprint,info,nfev,njev
c       integer ipvt(n)
c       double precision ftol,xtol,gtol,epsfcn,factor
c       double precision x(n),fvec(m),fjac(ldfjac,n),diag(n),qtf(n),
c      *     wa1(n),wa2(n),wa3(n),wa4(m)
c
c
c 3. parameters.
c
c       parameters designated as input parameters must be specified on
c       entry to dnls1 and are not changed on exit, while parameters
c       designated as output parameters need not be specified on entry
c       and are set to appropriate values on exit from dnls1.
c
c      fcn is the name of the user-supplied subroutine which calculate
c         the functions.  if the user wants to supply the jacobian
c         (iopt=2 or 3), then fcn must be written to calculate the
c         jacobian, as well as the functions.  see the explanation
c         of the iopt argument below.
c         if the user wants the iterates printed (nprint positive), then
c         fcn must do the printing.  see the explanation of nprint
c         below.  fcn must be declared in an external statement in the
c         calling program and should be written as follows.
c
c
c         subroutine fcn(iflag,m,n,x,fvec,fjac,ldfjac)
c         integer iflag,ldfjac,m,n
c         double precision x(n),fvec(m)
c         ----------
c         fjac and ldfjac may be ignored       , if iopt=1.
c         double precision fjac(ldfjac,n)      , if iopt=2.
c         double precision fjac(n)             , if iopt=3.
c         ----------
c           if iflag=0, the values in x and fvec are available
c           for printing.  see the explanation of nprint below.
c           iflag will never be zero unless nprint is positive.
c           the values of x and fvec must not be changed.
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
c           never be 3 unless iopt=3.  fvec contains the function
c           values at x and must not be altered.  fjac(j) must be
c           set to the derivative of fvec(ldfjac) with respect to x(j).
c         return
c         ----------
c         end
c
c
c         the value of iflag should not be changed by fcn unless the
c         user wants to terminate execution of dnls1.  in this case, set
c         iflag to a negative integer.
c
c
c       iopt is an input variable which specifies how the jacobian will
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
c       x is an array of length n.  on input, x must contain an initial
c         estimate of the solution vector.  on output, x contains the
c         final estimate of the solution vector.
c
c       fvec is an output array of length m which contains the functions
c         evaluated at the output x.
c
c       fjac is an output array.  for iopt=1 and 2, fjac is an m by n
c         array.  for iopt=3, fjac is an n by n array.  the upper n by n
c         submatrix of fjac contains an upper triangular matrix r with
c         diagonal elements of nonincreasing magnitude such that
c
c                t     t           t
c               p *(jac *jac)*p = r *r,
c
c         where p is a permutation matrix and jac is the final calcu-
c         lated jacobian.  column j of p is column ipvt(j) (see below)
c         of the identity matrix.  the lower part of fjac contains
c         information generated during the computation of r.
c
c       ldfjac is a positive integer input variable which specifies
c         the leading dimension of the array fjac.  for iopt=1 and 2,
c         ldfjac must not be less than m.  for iopt=3, ldfjac must not
c         be less than n.
c
c       ftol is a non-negative input variable.  termination occurs when
c         both the actual and predicted relative reductions in the sum
c         of squares are at most ftol.  therefore, ftol measures the
c         relative error desired in the sum of squares.  section 4 con-
c         tains more details about ftol.
c
c       xtol is a non-negative input variable.  termination occurs when
c         the relative error between two consecutive iterates is at most
c         xtol.  therefore, xtol measures the relative error desired in
c         the approximate solution.  section 4 contains more details
c         about xtol.
c
c       gtol is a non-negative input variable.  termination occurs when
c         the cosine of the angle between fvec and any column of the
c         jacobian is at most gtol in absolute value.  therefore, gtol
c         measures the orthogonality desired between the function vector
c         and the columns of the jacobian.  section 4 contains more
c         details about gtol.
c
c       maxfev is a positive integer input variable.  termination occurs
c         when the number of calls to fcn to evaluate the functions
c         has reached maxfev.
c
c       epsfcn is an input variable used in determining a suitable step
c         for the forward-difference approximation.  this approximation
c         assumes that the relative errors in the functions are of the
c         order of epsfcn.  if epsfcn is less than the machine preci-
c         sion, it is assumed that the relative errors in the functions
c         are of the order of the machine precision.  if iopt=2 or 3,
c         then epsfcn can be ignored (treat it as a dummy argument).
c
c       diag is an array of length n.  if mode = 1 (see below), diag is
c         internally set.  if mode = 2, diag must contain positive
c         entries that serve as implicit (multiplicative) scale factors
c         for the variables.
c
c       mode is an integer input variable.  if mode = 1, the variables
c         will be scaled internally.  if mode = 2, the scaling is speci-
c         fied by the input diag.  other values of mode are equivalent
c         to mode = 1.
c
c       factor is a positive input variable used in determining the ini-
c         tial step bound.  this bound is set to the product of factor
c         and the euclidean norm of diag*x if nonzero, or else to factor
c         itself.  in most cases factor should lie in the interval
c         (.1,100.).  100. is a generally recommended value.
c
c       nprint is an integer input variable that enables controlled
c         printing of iterates if it is positive.  in this case, fcn is
c         called with iflag = 0 at the beginning of the first iteration
c         and every nprint iterations thereafter and immediately prior
c         to return, with x and fvec available for printing. appropriate
c         print statements must be added to fcn (see example) and
c         fvec should not be altered.  if nprint is not positive, no
c         special calls to fcn with iflag = 0 are made.
c
c       info is an integer output variable.  if the user has terminated
c        execution, info is set to the (negative) value of iflag.  see
c        description of fcn and jac. otherwise, info is set as follows
c
c         info = 0  improper input parameters.
c
c         info = 1  both actual and predicted relative reductions in the
c                   sum of squares are at most ftol.
c
c         info = 2  relative error between two consecutive iterates is
c                   at most xtol.
c
c         info = 3  conditions for info = 1 and info = 2 both hold.
c
c         info = 4  the cosine of the angle between fvec and any column
c                   of the jacobian is at most gtol in absolute value.
c
c         info = 5  number of calls to fcn for function evaluation
c                   has reached maxfev.
c
c         info = 6  ftol is too small.  no further reduction in the sum
c                   of squares is possible.
c
c         info = 7  xtol is too small.  no further improvement in the
c                   approximate solution x is possible.
c
c         info = 8  gtol is too small.  fvec is orthogonal to the
c                   columns of the jacobian to machine precision.
c
c         sections 4 and 5 contain more details about info.
c
c       nfev is an integer output variable set to the number of calls to
c         fcn for function evaluation.
c
c       njev is an integer output variable set to the number of
c         evaluations of the full jacobian.  if iopt=2, only one call to
c         fcn is required for each evaluation of the full jacobian.
c         if iopt=3, the m calls to fcn are required.
c         if iopt=1, then njev is set to zero.
c
c       ipvt is an integer output array of length n.  ipvt defines a
c         permutation matrix p such that jac*p = q*r, where jac is the
c         final calculated jacobian, q is orthogonal (not stored), and r
c         is upper triangular with diagonal elements of nonincreasing
c         magnitude.  column j of p is column ipvt(j) of the identity
c         matrix.
c
c       qtf is an output array of length n which contains the first n
c         elements of the vector (q transpose)*fvec.
c
c       wa1, wa2, and wa3 are work arrays of length n.
c
c       wa4 is a work array of length m.
c
c
c 4. successful completion.
c
c       the accuracy of dnls1 is controlled by the convergence parame-
c       ters ftol, xtol, and gtol.  these parameters are used in tests
c       which make three types of comparisons between the approximation
c       x and a solution xsol.  dnls1 terminates when any of the tests
c       is satisfied.  if any of the convergence parameters is less than
c       the machine precision (as defined by the function r1mach(4)),
c       then dnls1 only attempts to satisfy the test defined by the
c       machine precision.  further progress is not usually possible.
c
c       the tests assume that the functions are reasonably well behaved,
c       and, if the jacobian is supplied by the user, that the functions
c       and the jacobian are coded consistently.  if these conditions
c       are not satisfied, then dnls1 may incorrectly indicate conver-
c       gence.  if the jacobian is coded correctly or iopt=1,
c       then the validity of the answer can be checked, for example, by
c       rerunning dnls1 with tighter tolerances.
c
c       first convergence test.  if enorm(z) denotes the euclidean norm
c         of a vector z, then this test attempts to guarantee that
c
c               enorm(fvec) .le. (1+ftol)*enorm(fvecs),
c
c         where fvecs denotes the functions evaluated at xsol.  if this
c         condition is satisfied with ftol = 10**(-k), then the final
c         residual norm enorm(fvec) has k significant decimal digits and
c         info is set to 1 (or to 3 if the second test is also satis-
c         fied).  unless high precision solutions are required, the
c         recommended value for ftol is the square root of the machine
c         precision.
c
c       second convergence test.  if d is the diagonal matrix whose
c         entries are defined by the array diag, then this test attempts
c         to guarantee that
c
c               enorm(d*(x-xsol)) .le. xtol*enorm(d*xsol).
c
c         if this condition is satisfied with xtol = 10**(-k), then the
c         larger components of d*x have k significant decimal digits and
c         info is set to 2 (or to 3 if the first test is also satis-
c         fied).  there is a danger that the smaller components of d*x
c         may have large relative errors, but if mode = 1, then the
c         accuracy of the components of x is usually related to their
c         sensitivity.  unless high precision solutions are required,
c         the recommended value for xtol is the square root of the
c         machine precision.
c
c       third convergence test.  this test is satisfied when the cosine
c         of the angle between fvec and any column of the jacobian at x
c         is at most gtol in absolute value.  there is no clear rela-
c         tionship between this test and the accuracy of dnls1, and
c         furthermore, the test is equally well satisfied at other crit-
c         ical points, namely maximizers and saddle points.  therefore,
c         termination caused by this test (info = 4) should be examined
c         carefully.  the recommended value for gtol is zero.
c
c
c 5. unsuccessful completion.
c
c       unsuccessful termination of dnls1 can be due to improper input
c       parameters, arithmetic interrupts, or an excessive number of
c       function evaluations.
c
c       improper input parameters.  info is set to 0 if iopt .lt. 1
c         or iopt .gt. 3, or n .le. 0, or m .lt. n, or for iopt=1 or 2
c         ldfjac .lt. m, or for iopt=3 ldfjac .lt. n, or ftol .lt. 0.e0,
c         or xtol .lt. 0.e0, or gtol .lt. 0.e0, or maxfev .le. 0, or
c         factor .le. 0.e0.
c
c       arithmetic interrupts.  if these interrupts occur in the fcn
c         subroutine during an early stage of the computation, they may
c         be caused by an unacceptable choice of x by dnls1.  in this
c         case, it may be possible to remedy the situation by rerunning
c         dnls1 with a smaller value of factor.
c
c       excessive number of function evaluations.  a reasonable value
c         for maxfev is 100*(n+1) for iopt=2 or 3 and 200*(n+1) for
c         iopt=1.  if the number of calls to fcn reaches maxfev, then
c         this indicates that the routine is converging very slowly
c         as measured by the progress of fvec, and info is set to 5.
c         in this case, it may be helpful to restart dnls1 with mode
c         set to 1.
c
c
c 6. characteristics of the algorithm.
c
c       dnls1 is a modification of the levenberg-marquardt algorithm.
c       two of its main characteristics involve the proper use of
c       implicitly scaled variables (if mode = 1) and an optimal choice
c       for the correction.  the use of implicitly scaled variables
c       achieves scale invariance of dnls1 and limits the size of the
c       correction in any direction where the functions are changing
c       rapidly.  the optimal choice of the correction guarantees (under
c       reasonable conditions) global convergence from starting points
c       far from the solution and a fast rate of convergence for
c       problems with small residuals.
c
c       timing.  the time required by dnls1 to solve a given problem
c         depends on m and n, the behavior of the functions, the accu-
c         racy requested, and the starting point.  the number of arith-
c         metic operations needed by dnls1 is about n**3 to process each
c         evaluation of the functions (call to fcn) and to process each
c         evaluation of the jacobian it takes m*n**2 for iopt=2 (one
c         call to fcn), m*n**2 for iopt=1 (n calls to fcn) and
c         1.5*m*n**2 for iopt=3 (m calls to fcn).  unless fcn
c         can be evaluated quickly, the timing of dnls1 will be
c         strongly influenced by the time spent in fcn.
c
c       storage.  dnls1 requires (m*n + 2*m + 6*n) for iopt=1 or 2 and
c         (n**2 + 2*m + 6*n) for iopt=3 single precision storage
c         locations and n integer storage locations, in addition to
c         the storage required by the program.  there are no internally
c         declared storage arrays.
c
c *long description:
c
c 7. example.
c
c       the problem is to determine the values of x(1), x(2), and x(3)
c       which provide the best fit (in the least squares sense) of
c
c             x(1) + u(i)/(v(i)*x(2) + w(i)*x(3)),  i = 1, 15
c
c       to the data
c
c             y = (0.14,0.18,0.22,0.25,0.29,0.32,0.35,0.39,
c                  0.37,0.58,0.73,0.96,1.34,2.10,4.39),
c
c       where u(i) = i, v(i) = 16 - i, and w(i) = min(u(i),v(i)).  the
c       i-th component of fvec is thus defined by
c
c             y(i) - (x(1) + u(i)/(v(i)*x(2) + w(i)*x(3))).
c
c       **********
c
c       program test
c c
c c     driver for dnls1 example.
c c
c       integer j,iopt,m,n,ldfjac,maxfev,mode,nprint,info,nfev,njev,
c      *        nwrite
c       integer ipvt(3)
c       double precision ftol,xtol,gtol,factor,fnorm,epsfcn
c       double precision x(3),fvec(15),fjac(15,3),diag(3),qtf(3),
c      *     wa1(3),wa2(3),wa3(3),wa4(15)
c       double precision denorm,d1mach
c       external fcn
c       data nwrite /6/
c c
c       iopt = 1
c       m = 15
c       n = 3
c c
c c     the following starting values provide a rough fit.
c c
c       x(1) = 1.e0
c       x(2) = 1.e0
c       x(3) = 1.e0
c c
c       ldfjac = 15
c c
c c     set ftol and xtol to the square root of the machine precision
c c     and gtol to zero.  unless high precision solutions are
c c     required, these are the recommended settings.
c c
c       ftol = sqrt(r1mach(4))
c       xtol = sqrt(r1mach(4))
c       gtol = 0.e0
c c
c       maxfev = 400
c       epsfcn = 0.0
c       mode = 1
c       factor = 1.e2
c       nprint = 0
c c
c       call dnls1(fcn,iopt,m,n,x,fvec,fjac,ldfjac,ftol,xtol,
c      *           gtol,maxfev,epsfcn,diag,mode,factor,nprint,
c      *           info,nfev,njev,ipvt,qtf,wa1,wa2,wa3,wa4)
c       fnorm = enorm(m,fvec)
c       write (nwrite,1000) fnorm,nfev,njev,info,(x(j),j=1,n)
c       stop
c  1000 format (5x,' final l2 norm of the residuals',e15.7 //
c      *        5x,' number of function evaluations',i10 //
c      *        5x,' number of jacobian evaluations',i10 //
c      *        5x,' exit parameter',16x,i10 //
c      *        5x,' final approximate solution' // 5x,3e15.7)
c       end
c       subroutine fcn(iflag,m,n,x,fvec,dum,idum)
c c     this is the form of the fcn routine if iopt=1,
c c     that is, if the user does not calculate the jacobian.
c       integer i,m,n,iflag
c       double precision x(n),fvec(m),y(15)
c       double precision tmp1,tmp2,tmp3,tmp4
c       data y(1),y(2),y(3),y(4),y(5),y(6),y(7),y(8),
c      *     y(9),y(10),y(11),y(12),y(13),y(14),y(15)
c      *     /1.4e-1,1.8e-1,2.2e-1,2.5e-1,2.9e-1,3.2e-1,3.5e-1,3.9e-1,
c      *      3.7e-1,5.8e-1,7.3e-1,9.6e-1,1.34e0,2.1e0,4.39e0/
c c
c       if (iflag .ne. 0) go to 5
c c
c c     insert print statements here when nprint is positive.
c c
c       return
c     5 continue
c       do 10 i = 1, m
c          tmp1 = i
c          tmp2 = 16 - i
c          tmp3 = tmp1
c          if (i .gt. 8) tmp3 = tmp2
c          fvec(i) = y(i) - (x(1) + tmp1/(x(2)*tmp2 + x(3)*tmp3))
c    10    continue
c       return
c       end
c
c
c       results obtained with different compilers or machines
c       may be slightly different.
c
c       final l2 norm of the residuals  0.9063596e-01
c
c       number of function evaluations        25
c
c       number of jacobian evaluations         0
c
c       exit parameter                         1
c
c       final approximate solution
c
c        0.8241058e-01  0.1133037e+01  0.2343695e+01
c
c
c       for iopt=2, fcn would be modified as follows to also
c       calculate the full jacobian when iflag=2.
c
c       subroutine fcn(iflag,m,n,x,fvec,fjac,ldfjac)
c c
c c     this is the form of the fcn routine if iopt=2,
c c     that is, if the user calculates the full jacobian.
c c
c       integer i,ldfjac,m,n,iflag
c       double precision x(n),fvec(m),fjac(ldfjac,n),y(15)
c       double precision tmp1,tmp2,tmp3,tmp4
c       data y(1),y(2),y(3),y(4),y(5),y(6),y(7),y(8),
c      *     y(9),y(10),y(11),y(12),y(13),y(14),y(15)
c      *     /1.4e-1,1.8e-1,2.2e-1,2.5e-1,2.9e-1,3.2e-1,3.5e-1,3.9e-1,
c      *      3.7e-1,5.8e-1,7.3e-1,9.6e-1,1.34e0,2.1e0,4.39e0/
c c
c       if (iflag .ne. 0) go to 5
c c
c c     insert print statements here when nprint is positive.
c c
c       return
c     5 continue
c       if(iflag.ne.1) go to 20
c       do 10 i = 1, m
c          tmp1 = i
c          tmp2 = 16 - i
c          tmp3 = tmp1
c          if (i .gt. 8) tmp3 = tmp2
c          fvec(i) = y(i) - (x(1) + tmp1/(x(2)*tmp2 + x(3)*tmp3))
c    10    continue
c       return
c c
c c     below, calculate the full jacobian.
c c
c    20    continue
c c
c       do 30 i = 1, m
c          tmp1 = i
c          tmp2 = 16 - i
c          tmp3 = tmp1
c          if (i .gt. 8) tmp3 = tmp2
c          tmp4 = (x(2)*tmp2 + x(3)*tmp3)**2
c          fjac(i,1) = -1.e0
c          fjac(i,2) = tmp1*tmp2/tmp4
c          fjac(i,3) = tmp1*tmp3/tmp4
c    30    continue
c       return
c       end
c
c
c       for iopt = 3, fjac would be dimensioned as fjac(3,3),
c         ldfjac would be set to 3, and fcn would be written as
c         follows to calculate a row of the jacobian when iflag=3.
c
c       subroutine fcn(iflag,m,n,x,fvec,fjac,ldfjac)
c c     this is the form of the fcn routine if iopt=3,
c c     that is, if the user calculates the jacobian row by row.
c       integer i,m,n,iflag
c       double precision x(n),fvec(m),fjac(n),y(15)
c       double precision tmp1,tmp2,tmp3,tmp4
c       data y(1),y(2),y(3),y(4),y(5),y(6),y(7),y(8),
c      *     y(9),y(10),y(11),y(12),y(13),y(14),y(15)
c      *     /1.4e-1,1.8e-1,2.2e-1,2.5e-1,2.9e-1,3.2e-1,3.5e-1,3.9e-1,
c      *      3.7e-1,5.8e-1,7.3e-1,9.6e-1,1.34e0,2.1e0,4.39e0/
c c
c       if (iflag .ne. 0) go to 5
c c
c c     insert print statements here when nprint is positive.
c c
c       return
c     5 continue
c       if( iflag.ne.1) go to 20
c       do 10 i = 1, m
c          tmp1 = i
c          tmp2 = 16 - i
c          tmp3 = tmp1
c          if (i .gt. 8) tmp3 = tmp2
c          fvec(i) = y(i) - (x(1) + tmp1/(x(2)*tmp2 + x(3)*tmp3))
c    10    continue
c       return
c c
c c     below, calculate the ldfjac-th row of the jacobian.
c c
c    20 continue
c
c       i = ldfjac
c          tmp1 = i
c          tmp2 = 16 - i
c          tmp3 = tmp1
c          if (i .gt. 8) tmp3 = tmp2
c          tmp4 = (x(2)*tmp2 + x(3)*tmp3)**2
c          fjac(1) = -1.e0
c          fjac(2) = tmp1*tmp2/tmp4
c          fjac(3) = tmp1*tmp3/tmp4
c       return
c       end
c
c***references  jorge j. more, the levenberg-marquardt algorithm:
c                 implementation and theory.  in numerical analysis
c                 proceedings (dundee, june 28 - july 1, 1977, g. a.
c                 watson, editor), lecture notes in mathematics 630,
c                 springer-verlag, 1978.
c***routines called  d1mach, dckder, denorm, dfdjc3, dmpar, dqrfac,
c                    dwupdt, xermsg
c***revision history  (yymmdd)
c   800301  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   891006  cosmetic changes to prologue.  (wrb)
c   891006  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900510  convert xerrwv calls to xermsg calls.  (rwc)
c   920205  corrected xern1 declaration.  (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  dnls1
      implicit double precision (a-h,o-z)
      integer iopt,m,n,ldfjac,maxfev,mode,nprint,info,nfev,njev
      integer ijunk,nrow,ipvt(*)
      double precision ftol,xtol,gtol,factor,epsfcn
      double precision x(*),fvec(*),fjac(ldfjac,*),diag(*),qtf(*),
     1     wa1(*),wa2(*),wa3(*),wa4(*)
      logical sing
      external fcn
      integer i,iflag,iter,j,l,modech
      double precision actred,delta,dirder,epsmch,fnorm,fnorm1,gnorm,
     1     one,par,pnorm,prered,p1,p5,p25,p75,p0001,ratio,sum,temp,
     2     temp1,temp2,xnorm,zero
      double precision d1mach,denorm,err,chklim
      character*8 xern1
      character*16 xern3
      save chklim, one, p1, p5, p25, p75, p0001, zero
c
      data chklim/.1d0/
      data one,p1,p5,p25,p75,p0001,zero
     1     /1.0d0,1.0d-1,5.0d-1,2.5d-1,7.5d-1,1.0d-4,0.0d0/
c***first executable statement  dnls1
      epsmch = d1mach(4)
c
      info = 0
      iflag = 0
      nfev = 0
      njev = 0
c
c     check the input parameters for errors.
c
      if (iopt .lt. 1 .or. iopt .gt. 3 .or. n .le. 0 .or.
     1    m .lt. n .or. ldfjac .lt. n .or. ftol .lt. zero
     2    .or. xtol .lt. zero .or. gtol .lt. zero
     3    .or. maxfev .le. 0 .or. factor .le. zero) go to 300
      if (iopt .lt. 3 .and. ldfjac .lt. m) go to 300
      if (mode .ne. 2) go to 20
      do 10 j = 1, n
         if (diag(j) .le. zero) go to 300
   10    continue
   20 continue
c
c     evaluate the function at the starting point
c     and calculate its norm.
c
      iflag = 1
      ijunk = 1
      call fcn(iflag,m,n,x,fvec,fjac,ijunk)
      nfev = 1
      if (iflag .lt. 0) go to 300
      fnorm = denorm(m,fvec)
c
c     initialize levenberg-marquardt parameter and iteration counter.
c
      par = zero
      iter = 1
c
c     beginning of the outer loop.
c
   30 continue
c
c        if requested, call fcn to enable printing of iterates.
c
         if (nprint .le. 0) go to 40
         iflag = 0
         if (mod(iter-1,nprint) .eq. 0)
     1      call fcn(iflag,m,n,x,fvec,fjac,ijunk)
         if (iflag .lt. 0) go to 300
   40    continue
c
c        calculate the jacobian matrix.
c
      if (iopt .eq. 3) go to 475
c
c     store the full jacobian using m*n storage
c
      if (iopt .eq. 1) go to 410
c
c     the user supplies the jacobian
c
         iflag = 2
         call fcn(iflag,m,n,x,fvec,fjac,ldfjac)
         njev = njev + 1
c
c             on the first iteration, check the user supplied jacobian
c
         if (iter .le. 1) then
            if (iflag .lt. 0) go to 300
c
c           get the incremented x-values into wa1(*).
c
            modech = 1
            call dckder(m,n,x,fvec,fjac,ldfjac,wa1,wa4,modech,err)
c
c           evaluate function at incremented value and put in wa4(*).
c
            iflag = 1
            call fcn(iflag,m,n,wa1,wa4,fjac,ldfjac)
            nfev = nfev + 1
            if(iflag .lt. 0) go to 300
            do 350 i = 1, m
               modech = 2
               call dckder(1,n,x,fvec(i),fjac(i,1),ldfjac,wa1,
     1              wa4(i),modech,err)
               if (err .lt. chklim) then
                  write (xern1, '(i8)') i
                  write (xern3, '(1pe15.6)') err
                  call xermsg ('slatec', 'dnls1', 'derivative of ' //
     *               'function ' // xern1 // ' may be wrong, err = ' //
     *               xern3 // ' too close to 0.', 7, 0)
               endif
  350       continue
         endif
c
         go to 420
c
c     the code approximates the jacobian
c
410      iflag = 1
         call dfdjc3(fcn,m,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa4)
         nfev = nfev + n
  420    if (iflag .lt. 0) go to 300
c
c        compute the qr factorization of the jacobian.
c
         call dqrfac(m,n,fjac,ldfjac,.true.,ipvt,n,wa1,wa2,wa3)
c
c        form (q transpose)*fvec and store the first n components in
c        qtf.
c
         do 430 i = 1, m
            wa4(i) = fvec(i)
  430         continue
         do 470 j = 1, n
            if (fjac(j,j) .eq. zero) go to 460
            sum = zero
            do 440 i = j, m
               sum = sum + fjac(i,j)*wa4(i)
  440          continue
            temp = -sum/fjac(j,j)
            do 450 i = j, m
               wa4(i) = wa4(i) + fjac(i,j)*temp
  450          continue
  460       continue
            fjac(j,j) = wa1(j)
            qtf(j) = wa4(j)
  470       continue
         go to 560
c
c        accumulate the jacobian by rows in order to save storage.
c        compute the qr factorization of the jacobian matrix
c        calculated one row at a time, while simultaneously
c        forming (q transpose)*fvec and storing the first
c        n components in qtf.
c
  475    do 490 j = 1, n
            qtf(j) = zero
            do 480 i = 1, n
               fjac(i,j) = zero
  480          continue
  490        continue
         do 500 i = 1, m
            nrow = i
            iflag = 3
            call fcn(iflag,m,n,x,fvec,wa3,nrow)
            if (iflag .lt. 0) go to 300
c
c            on the first iteration, check the user supplied jacobian.
c
            if(iter .gt. 1) go to 498
c
c            get the incremented x-values into wa1(*).
c
            modech = 1
            call dckder(m,n,x,fvec,fjac,ldfjac,wa1,wa4,modech,err)
c
c            evaluate at incremented values, if not already evaluated.
c
            if(i .ne. 1) go to 495
c
c            evaluate function at incremented value and put into wa4(*).
c
            iflag = 1
            call fcn(iflag,m,n,wa1,wa4,fjac,nrow)
            nfev = nfev + 1
            if(iflag .lt. 0) go to 300
495         continue
            modech = 2
            call dckder(1,n,x,fvec(i),wa3,1,wa1,wa4(i),modech,err)
            if (err .lt. chklim) then
               write (xern1, '(i8)') i
               write (xern3, '(1pe15.6)') err
               call xermsg ('slatec', 'dnls1', 'derivative of function '
     *            // xern1 // ' may be wrong, err = ' // xern3 //
     *            ' too close to 0.', 7, 0)
            endif
498         continue
c
            temp = fvec(i)
            call dwupdt(n,fjac,ldfjac,wa3,qtf,temp,wa1,wa2)
  500       continue
         njev = njev + 1
c
c        if the jacobian is rank deficient, call dqrfac to
c        reorder its columns and update the components of qtf.
c
         sing = .false.
         do 510 j = 1, n
            if (fjac(j,j) .eq. zero) sing = .true.
            ipvt(j) = j
            wa2(j) = denorm(j,fjac(1,j))
  510       continue
         if (.not.sing) go to 560
         call dqrfac(n,n,fjac,ldfjac,.true.,ipvt,n,wa1,wa2,wa3)
         do 550 j = 1, n
            if (fjac(j,j) .eq. zero) go to 540
            sum = zero
            do 520 i = j, n
               sum = sum + fjac(i,j)*qtf(i)
  520         continue
            temp = -sum/fjac(j,j)
            do 530 i = j, n
               qtf(i) = qtf(i) + fjac(i,j)*temp
  530          continue
  540       continue
            fjac(j,j) = wa1(j)
  550       continue
  560    continue
c
c        on the first iteration and if mode is 1, scale according
c        to the norms of the columns of the initial jacobian.
c
         if (iter .ne. 1) go to 80
         if (mode .eq. 2) go to 60
         do 50 j = 1, n
            diag(j) = wa2(j)
            if (wa2(j) .eq. zero) diag(j) = one
   50       continue
   60    continue
c
c        on the first iteration, calculate the norm of the scaled x
c        and initialize the step bound delta.
c
         do 70 j = 1, n
            wa3(j) = diag(j)*x(j)
   70       continue
         xnorm = denorm(n,wa3)
         delta = factor*xnorm
         if (delta .eq. zero) delta = factor
   80    continue
c
c        compute the norm of the scaled gradient.
c
         gnorm = zero
         if (fnorm .eq. zero) go to 170
         do 160 j = 1, n
            l = ipvt(j)
            if (wa2(l) .eq. zero) go to 150
            sum = zero
            do 140 i = 1, j
               sum = sum + fjac(i,j)*(qtf(i)/fnorm)
  140          continue
            gnorm = max(gnorm,abs(sum/wa2(l)))
  150       continue
  160       continue
  170    continue
c
c        test for convergence of the gradient norm.
c
         if (gnorm .le. gtol) info = 4
         if (info .ne. 0) go to 300
c
c        rescale if necessary.
c
         if (mode .eq. 2) go to 190
         do 180 j = 1, n
            diag(j) = max(diag(j),wa2(j))
  180       continue
  190    continue
c
c        beginning of the inner loop.
c
  200    continue
c
c           determine the levenberg-marquardt parameter.
c
            call dmpar(n,fjac,ldfjac,ipvt,diag,qtf,delta,par,wa1,wa2,
     1                 wa3,wa4)
c
c           store the direction p and x + p. calculate the norm of p.
c
            do 210 j = 1, n
               wa1(j) = -wa1(j)
               wa2(j) = x(j) + wa1(j)
               wa3(j) = diag(j)*wa1(j)
  210          continue
            pnorm = denorm(n,wa3)
c
c           on the first iteration, adjust the initial step bound.
c
            if (iter .eq. 1) delta = min(delta,pnorm)
c
c           evaluate the function at x + p and calculate its norm.
c
            iflag = 1
            call fcn(iflag,m,n,wa2,wa4,fjac,ijunk)
            nfev = nfev + 1
            if (iflag .lt. 0) go to 300
            fnorm1 = denorm(m,wa4)
c
c           compute the scaled actual reduction.
c
            actred = -one
            if (p1*fnorm1 .lt. fnorm) actred = one - (fnorm1/fnorm)**2
c
c           compute the scaled predicted reduction and
c           the scaled directional derivative.
c
            do 230 j = 1, n
               wa3(j) = zero
               l = ipvt(j)
               temp = wa1(l)
               do 220 i = 1, j
                  wa3(i) = wa3(i) + fjac(i,j)*temp
  220             continue
  230          continue
            temp1 = denorm(n,wa3)/fnorm
            temp2 = (sqrt(par)*pnorm)/fnorm
            prered = temp1**2 + temp2**2/p5
            dirder = -(temp1**2 + temp2**2)
c
c           compute the ratio of the actual to the predicted
c           reduction.
c
            ratio = zero
            if (prered .ne. zero) ratio = actred/prered
c
c           update the step bound.
c
            if (ratio .gt. p25) go to 240
               if (actred .ge. zero) temp = p5
               if (actred .lt. zero)
     1            temp = p5*dirder/(dirder + p5*actred)
               if (p1*fnorm1 .ge. fnorm .or. temp .lt. p1) temp = p1
               delta = temp*min(delta,pnorm/p1)
               par = par/temp
               go to 260
  240       continue
               if (par .ne. zero .and. ratio .lt. p75) go to 250
               delta = pnorm/p5
               par = p5*par
  250          continue
  260       continue
c
c           test for successful iteration.
c
            if (ratio .lt. p0001) go to 290
c
c           successful iteration. update x, fvec, and their norms.
c
            do 270 j = 1, n
               x(j) = wa2(j)
               wa2(j) = diag(j)*x(j)
  270          continue
            do 280 i = 1, m
               fvec(i) = wa4(i)
  280          continue
            xnorm = denorm(n,wa2)
            fnorm = fnorm1
            iter = iter + 1
  290       continue
c
c           tests for convergence.
c
            if (abs(actred) .le. ftol .and. prered .le. ftol
     1          .and. p5*ratio .le. one) info = 1
            if (delta .le. xtol*xnorm) info = 2
            if (abs(actred) .le. ftol .and. prered .le. ftol
     1          .and. p5*ratio .le. one .and. info .eq. 2) info = 3
            if (info .ne. 0) go to 300
c
c           tests for termination and stringent tolerances.
c
            if (nfev .ge. maxfev) info = 5
            if (abs(actred) .le. epsmch .and. prered .le. epsmch
     1          .and. p5*ratio .le. one) info = 6
            if (delta .le. epsmch*xnorm) info = 7
            if (gnorm .le. epsmch) info = 8
            if (info .ne. 0) go to 300
c
c           end of the inner loop. repeat if iteration unsuccessful.
c
            if (ratio .lt. p0001) go to 200
c
c        end of the outer loop.
c
         go to 30
  300 continue
c
c     termination, either normal or user imposed.
c
      if (iflag .lt. 0) info = iflag
      iflag = 0
      if (nprint .gt. 0) call fcn(iflag,m,n,x,fvec,fjac,ijunk)
      if (info .lt. 0) call xermsg ('slatec', 'dnls1',
     +   'execution terminated because user set iflag negative.', 1, 1)
      if (info .eq. 0) call xermsg ('slatec', 'dnls1',
     +   'invalid input parameter.', 2, 1)
      if (info .eq. 4) call xermsg ('slatec', 'dnls1',
     +   'third convergence condition, check results before accepting.',
     +   1, 1)
      if (info .eq. 5) call xermsg ('slatec', 'dnls1',
     +   'too many function evaluations.', 9, 1)
      if (info .ge. 6) call xermsg ('slatec', 'dnls1',
     +   'tolerances too small, no further improvement possible.', 3, 1)
      return
c
c     last card of subroutine dnls1.
c
      end
