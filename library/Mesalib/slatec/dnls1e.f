*deck dnls1e
      subroutine dnls1e (fcn, iopt, m, n, x, fvec, tol, nprint, info,
     +   iw, wa, lwa)
c***begin prologue  dnls1e
c***purpose  an easy-to-use code which minimizes the sum of the squares
c            of m nonlinear functions in n variables by a modification
c            of the levenberg-marquardt algorithm.
c***library   slatec
c***category  k1b1a1, k1b1a2
c***type      double precision (snls1e-s, dnls1e-d)
c***keywords  easy-to-use, levenberg-marquardt, nonlinear data fitting,
c             nonlinear least squares
c***author  hiebert, k. l., (snla)
c***description
c
c 1. purpose.
c
c       the purpose of dnls1e is to minimize the sum of the squares of m
c       nonlinear functions in n variables by a modification of the
c       levenberg-marquardt algorithm.  this is done by using the more
c       general least-squares solver dnls1.  the user must provide a
c       subroutine which calculates the functions.  the user has the
c       option of how the jacobian will be supplied.  the user can
c       supply the full jacobian, or the rows of the jacobian (to avoid
c       storing the full jacobian), or let the code approximate the
c       jacobian by forward-differencing.  this code is the combination
c       of the minpack codes (argonne) lmder1, lmdif1, and lmstr1.
c
c
c 2. subroutine and type statements.
c
c       subroutine dnls1e(fcn,iopt,m,n,x,fvec,tol,nprint,
c      *                  info,iw,wa,lwa)
c       integer iopt,m,n,nprint,info,lwac,iw(n)
c       double precision tol,x(n),fvec(m),wa(lwa)
c       external fcn
c
c
c 3. parameters. all type real parameters are double precision
c
c       parameters designated as input parameters must be specified on
c       entry to dnls1e and are not changed on exit, while parameters
c       designated as output parameters need not be specified on entry
c       and are set to appropriate values on exit from dnls1e.
c
c      fcn is the name of the user-supplied subroutine which calculates
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
c         user wants to terminate execution of dnls1e.  in this case,
c         set iflag to a negative integer.
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
c       tol is a non-negative input variable.  termination occurs when
c         the algorithm estimates either that the relative error in the
c         sum of squares is at most tol or that the relative error
c         between x and the solution is at most tol.  section 4 contains
c         more details about tol.
c
c       nprint is an integer input variable that enables controlled
c         printing of iterates if it is positive.  in this case, fcn is
c         called with iflag = 0 at the beginning of the first iteration
c         and every nprint iterations thereafter and immediately prior
c         to return, with x and fvec available for printing. appropriate
c         print statements must be added to fcn (see example) and
c         fvec should not be altered.  if nprint is not positive, no
c         special calls of fcn with iflag = 0 are made.
c
c       info is an integer output variable.  if the user has terminated
c        execution, info is set to the (negative) value of iflag.  see
c        description of fcn and jac. otherwise, info is set as follows.
c
c         info = 0  improper input parameters.
c
c         info = 1  algorithm estimates that the relative error in the
c                   sum of squares is at most tol.
c
c         info = 2  algorithm estimates that the relative error between
c                   x and the solution is at most tol.
c
c         info = 3  conditions for info = 1 and info = 2 both hold.
c
c         info = 4  fvec is orthogonal to the columns of the jacobian to
c                   machine precision.
c
c         info = 5  number of calls to fcn has reached 100*(n+1)
c                   for iopt=2 or 3 or 200*(n+1) for iopt=1.
c
c         info = 6  tol is too small.  no further reduction in the sum
c                   of squares is possible.
c
c         info = 7  tol is too small.  no further improvement in the
c                   approximate solution x is possible.
c
c         sections 4 and 5 contain more details about info.
c
c       iw is an integer work array of length n.
c
c       wa is a work array of length lwa.
c
c       lwa is a positive integer input variable not less than
c         n*(m+5)+m for iopt=1 and 2 or n*(n+5)+m for iopt=3.
c
c
c 4. successful completion.
c
c       the accuracy of dnls1e is controlled by the convergence parame-
c       ter tol.  this parameter is used in tests which make three types
c       of comparisons between the approximation x and a solution xsol.
c       dnls1e terminates when any of the tests is satisfied.  if tol is
c       less than the machine precision (as defined by the function
c       r1mach(4)), then dnls1e only attempts to satisfy the test
c       defined by the machine precision.  further progress is not usu-
c       ally possible.  unless high precision solutions are required,
c       the recommended value for tol is the square root of the machine
c       precision.
c
c       the tests assume that the functions are reasonably well behaved,
c       and, if the jacobian is supplied by the user, that the functions
c       and the jacobian are coded consistently.  if these conditions
c       are not satisfied, then dnls1e may incorrectly indicate conver-
c       gence.  if the jacobian is coded correctly or iopt=1,
c       then the validity of the answer can be checked, for example, by
c       rerunning dnls1e with tighter tolerances.
c
c       first convergence test.  if enorm(z) denotes the euclidean norm
c         of a vector z, then this test attempts to guarantee that
c
c               enorm(fvec) .le. (1+tol)*enorm(fvecs),
c
c         where fvecs denotes the functions evaluated at xsol.  if this
c         condition is satisfied with tol = 10**(-k), then the final
c         residual norm enorm(fvec) has k significant decimal digits and
c         info is set to 1 (or to 3 if the second test is also satis-
c         fied).
c
c       second convergence test.  if d is a diagonal matrix (implicitly
c         generated by dnls1e) whose entries contain scale factors for
c         the variables, then this test attempts to guarantee that
c
c               enorm(d*(x-xsol)) .le.  tol*enorm(d*xsol).
c
c         if this condition is satisfied with tol = 10**(-k), then the
c         larger components of d*x have k significant decimal digits and
c         info is set to 2 (or to 3 if the first test is also satis-
c         fied).  there is a danger that the smaller components of d*x
c         may have large relative errors, but the choice of d is such
c         that the accuracy of the components of x is usually related to
c         their sensitivity.
c
c       third convergence test.  this test is satisfied when fvec is
c         orthogonal to the columns of the jacobian to machine preci-
c         sion.  there is no clear relationship between this test and
c         the accuracy of dnls1e, and furthermore, the test is equally
c         well satisfied at other critical points, namely maximizers and
c         saddle points.  therefore, termination caused by this test
c         (info = 4) should be examined carefully.
c
c
c 5. unsuccessful completion.
c
c       unsuccessful termination of dnls1e can be due to improper input
c       parameters, arithmetic interrupts, or an excessive number of
c       function evaluations.
c
c       improper input parameters.  info is set to 0 if iopt .lt. 1
c         or iopt .gt. 3, or n .le. 0, or m .lt. n, or tol .lt. 0.e0,
c         or for iopt=1 or 2 lwa .lt. n*(m+5)+m, or for iopt=3
c         lwa .lt. n*(n+5)+m.
c
c       arithmetic interrupts.  if these interrupts occur in the fcn
c         subroutine during an early stage of the computation, they may
c         be caused by an unacceptable choice of x by dnls1e.  in this
c         case, it may be possible to remedy the situation by not evalu-
c         ating the functions here, but instead setting the components
c         of fvec to numbers that exceed those in the initial fvec.
c
c       excessive number of function evaluations.  if the number of
c         calls to fcn reaches 100*(n+1) for iopt=2 or 3 or 200*(n+1)
c         for iopt=1, then this indicates that the routine is converging
c         very slowly as measured by the progress of fvec, and info is
c         set to 5.  in this case, it may be helpful to restart dnls1e,
c         thereby forcing it to disregard old (and possibly harmful)
c         information.
c
c
c 6. characteristics of the algorithm.
c
c       dnls1e is a modification of the levenberg-marquardt algorithm.
c       two of its main characteristics involve the proper use of
c       implicitly scaled variables and an optimal choice for the cor-
c       rection.  the use of implicitly scaled variables achieves scale
c       invariance of dnls1e and limits the size of the correction in
c       any direction where the functions are changing rapidly.  the
c       optimal choice of the correction guarantees (under reasonable
c       conditions) global convergence from starting points far from the
c       solution and a fast rate of convergence for problems with small
c       residuals.
c
c       timing.  the time required by dnls1e to solve a given problem
c         depends on m and n, the behavior of the functions, the accu-
c         racy requested, and the starting point.  the number of arith-
c         metic operations needed by dnls1e is about n**3 to process
c         each evaluation of the functions (call to fcn) and to process
c         each evaluation of the jacobian dnls1e takes m*n**2 for iopt=2
c         (one call to jac), m*n**2 for iopt=1 (n calls to fcn) and
c         1.5*m*n**2 for iopt=3 (m calls to fcn).  unless fcn
c         can be evaluated quickly, the timing of dnls1e will be
c         strongly influenced by the time spent in fcn.
c
c       storage.  dnls1e requires (m*n + 2*m + 6*n) for iopt=1 or 2 and
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
c c     driver for dnls1e example.
c c
c       integer i,iopt,m,n,nprint,jnfo,lwa,nwrite
c       integer iw(3)
c       double precision tol,fnorm,x(3),fvec(15),wa(75)
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
c       lwa = 75
c       nprint = 0
c c
c c     set tol to the square root of the machine precision.
c c     unless high precision solutions are required,
c c     this is the recommended setting.
c c
c       tol = sqrt(r1mach(4))
c c
c       call dnls1e(fcn,iopt,m,n,x,fvec,tol,nprint,
c      *            info,iw,wa,lwa)
c       fnorm = enorm(m,fvec)
c       write (nwrite,1000) fnorm,info,(x(j),j=1,n)
c       stop
c  1000 format (5x,' final l2 norm of the residuals',e15.7 //
c      *        5x,' exit
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
c***routines called  dnls1, xermsg
c***revision history  (yymmdd)
c   800301  date written
c   890831  modified array declarations.  (wrb)
c   891006  cosmetic changes to prologue.  (wrb)
c   891006  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   920501  reformatted the references section.  (wrb)
c***end prologue  dnls1e
      implicit double precision (a-h,o-z)
      integer m,n,nprint,info,lwa,iopt
      integer index,iw(*)
      double precision tol
      double precision x(*),fvec(*),wa(*)
      external fcn
      integer maxfev,mode,nfev,njev
      double precision factor,ftol,gtol,xtol,zero,epsfcn
      save factor, zero
      data factor,zero /1.0d2,0.0d0/
c***first executable statement  dnls1e
      info = 0
c
c     check the input parameters for errors.
c
      if (iopt .lt. 1 .or. iopt .gt. 3 .or.
     1    n .le. 0 .or. m .lt. n .or. tol .lt. zero
     2    .or. lwa .lt. n*(n+5) + m) go to 10
      if (iopt .lt. 3 .and. lwa .lt. n*(m+5) + m) go to 10
c
c     call dnls1.
c
      maxfev = 100*(n + 1)
      if (iopt .eq. 1) maxfev = 2*maxfev
      ftol = tol
      xtol = tol
      gtol = zero
      epsfcn = zero
      mode = 1
      index = 5*n+m
      call dnls1(fcn,iopt,m,n,x,fvec,wa(index+1),m,ftol,xtol,gtol,
     1           maxfev,epsfcn,wa(1),mode,factor,nprint,info,nfev,njev,
     2           iw,wa(n+1),wa(2*n+1),wa(3*n+1),wa(4*n+1),wa(5*n+1))
      if (info .eq. 8) info = 4
   10 continue
      if (info .eq. 0) call xermsg ('slatec', 'dnls1e',
     +   'invalid input parameter.', 2, 1)
      return
c
c     last card of subroutine dnls1e.
c
      end
