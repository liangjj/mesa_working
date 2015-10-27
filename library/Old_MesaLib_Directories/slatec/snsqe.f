*deck snsqe
      subroutine snsqe (fcn, jac, iopt, n, x, fvec, tol, nprint, info,
     +   wa, lwa)
c***begin prologue  snsqe
c***purpose  an easy-to-use code to find a zero of a system of n
c            nonlinear functions in n variables by a modification of
c            the powell hybrid method.
c***library   slatec
c***category  f2a
c***type      single precision (snsqe-s, dnsqe-d)
c***keywords  easy-to-use, nonlinear square system,
c             powell hybrid method, zeros
c***author  hiebert, k. l., (snla)
c***description
c
c 1. purpose.
c
c
c       the purpose of snsqe is to find a zero of a system of n non-
c       linear functions in n variables by a modification of the powell
c       hybrid method.  this is done by using the more general nonlinear
c       equation solver snsq.  the user must provide a subroutine which
c       calculates the functions.  the user has the option of either to
c       provide a subroutine which calculates the jacobian or to let the
c       code calculate it by a forward-difference approximation.  this
c       code is the combination of the minpack codes (argonne) hybrd1
c       and hybrj1.
c
c
c 2. subroutine and type statements.
c
c       subroutine snsqe(fcn,jac,iopt,n,x,fvec,tol,nprint,info,
c      *                  wa,lwa)
c       integer iopt,n,nprint,info,lwa
c       real tol
c       real x(n),fvec(n),wa(lwa)
c       external fcn,jac
c
c
c 3. parameters.
c
c       parameters designated as input parameters must be specified on
c       entry to snsqe and are not changed on exit, while parameters
c       designated as output parameters need not be specified on entry
c       and are set to appropriate values on exit from snsqe.
c
c       fcn is the name of the user-supplied subroutine which calculates
c         the functions.  fcn must be declared in an external statement
c         in the user calling program, and should be written as follows.
c
c         subroutine fcn(n,x,fvec,iflag)
c         integer n,iflag
c         real x(n),fvec(n)
c         ----------
c         calculate the functions at x and
c         return this vector in fvec.
c         ----------
c         return
c         end
c
c         the value of iflag should not be changed by fcn unless the
c         user wants to terminate execution of snsqe.  in this case, set
c         iflag to a negative integer.
c
c       jac is the name of the user-supplied subroutine which calculates
c         the jacobian.  if iopt=1, then jac must be declared in an
c         external statement in the user calling program, and should be
c         written as follows.
c
c         subroutine jac(n,x,fvec,fjac,ldfjac,iflag)
c         integer n,ldfjac,iflag
c         real x(n),fvec(n),fjac(ldfjac,n)
c         ----------
c         calculate the jacobian at x and return this
c         matrix in fjac.  fvec contains the function
c         values at x and should not be altered.
c         ----------
c         return
c         end
c
c         the value of iflag should not be changed by jac unless the
c         user wants to terminate execution of snsqe.  in this case, set
c         iflag to a negative integer.
c
c         if iopt=2, jac can be ignored (treat it as a dummy argument).
c
c       iopt is an input variable which specifies how the jacobian will
c         be calculated.  if iopt=1, then the user must supply the
c         jacobian through the subroutine jac.  if iopt=2, then the
c         code will approximate the jacobian by forward-differencing.
c
c       n is a positive integer input variable set to the number of
c         functions and variables.
c
c       x is an array of length n.  on input, x must contain an initial
c         estimate of the solution vector.  on output, x contains the
c         final estimate of the solution vector.
c
c       fvec is an output array of length n which contains the functions
c         evaluated at the output x.
c
c       tol is a non-negative input variable.  termination occurs when
c         the algorithm estimates that the relative error between x and
c         the solution is at most tol.  section 4 contains more details
c         about tol.
c
c       nprint is an integer input variable that enables controlled
c         printing of iterates if it is positive.  in this case, fcn is
c         called with iflag = 0 at the beginning of the first iteration
c         and every nprint iteration thereafter and immediately prior
c         to return, with x and fvec available for printing. appropriate
c         print statements must be added to fcn (see example). if nprint
c         is not positive, no special calls of fcn with iflag = 0 are
c         made.
c
c       info is an integer output variable.  if the user has terminated
c         execution, info is set to the (negative) value of iflag.  see
c         description of fcn and jac. otherwise, info is set as follows.
c
c         info = 0  improper input parameters.
c
c         info = 1  algorithm estimates that the relative error between
c                   x and the solution is at most tol.
c
c         info = 2  number of calls to fcn has reached or exceeded
c                   100*(n+1) for iopt=1 or 200*(n+1) for iopt=2.
c
c         info = 3  tol is too small.  no further improvement in the
c                   approximate solution x is possible.
c
c         info = 4  iteration is not making good progress.
c
c         sections 4 and 5 contain more details about info.
c
c       wa is a work array of length lwa.
c
c       lwa is a positive integer input variable not less than
c         (3*n**2+13*n))/2.
c
c
c 4. successful completion.
c
c       the accuracy of snsqe is controlled by the convergence parame-
c       ter tol.  this parameter is used in a test which makes a compar-
c       ison between the approximation x and a solution xsol.  snsqe
c       terminates when the test is satisfied.  if tol is less than the
c       machine precision (as defined by the function r1mach(4)), then
c       snsqe attempts only to satisfy the test defined by the machine
c       precision.  further progress is not usually possible.  unless
c       high precision solutions are required, the recommended value
c       for tol is the square root of the machine precision.
c
c       the test assumes that the functions are reasonably well behaved,
c       and, if the jacobian is supplied by the user, that the functions
c       and the jacobian  coded consistently.  if these conditions
c       are not satisfied, snsqe may incorrectly indicate convergence.
c       the coding of the jacobian can be checked by the subroutine
c       chkder.  if the jacobian is coded correctly or iopt=2, then
c       the validity of the answer can be checked, for example, by
c       rerunning snsqe with a tighter tolerance.
c
c       convergence test.  if enorm(z) denotes the euclidean norm of a
c         vector z, then this test attempts to guarantee that
c
c               enorm(x-xsol) .le.  tol*enorm(xsol).
c
c         if this condition is satisfied with tol = 10**(-k), then the
c         larger components of x have k significant decimal digits and
c         info is set to 1.  there is a danger that the smaller compo-
c         nents of x may have large relative errors, but the fast rate
c         of convergence of snsqe usually avoids this possibility.
c
c
c 5. unsuccessful completion.
c
c       unsuccessful termination of snsqe can be due to improper input
c       parameters, arithmetic interrupts, an excessive number of func-
c       tion evaluations, errors in the functions, or lack of good prog-
c       ress.
c
c       improper input parameters.  info is set to 0 if iopt .lt. 1, or
c         iopt .gt. 2, or n .le. 0, or tol .lt. 0.e0, or
c         lwa .lt. (3*n**2+13*n)/2.
c
c       arithmetic interrupts.  if these interrupts occur in the fcn
c         subroutine during an early stage of the computation, they may
c         be caused by an unacceptable choice of x by snsqe.  in this
c         case, it may be possible to remedy the situation by not evalu-
c         ating the functions here, but instead setting the components
c         of fvec to numbers that exceed those in the initial fvec.
c
c       excessive number of function evaluations.  if the number of
c         calls to fcn reaches 100*(n+1) for iopt=1 or 200*(n+1) for
c         iopt=2, then this indicates that the routine is converging
c         very slowly as measured by the progress of fvec, and info is
c         set to 2.  this situation should be unusual because, as
c         indicated below, lack of good progress is usually diagnosed
c         earlier by snsqe, causing termination with info = 4.
c
c       errors in the functions.  when iopt=2, the choice of step length
c         in the forward-difference approximation to the jacobian
c         assumes that the relative errors in the functions are of the
c         order of the machine precision.  if this is not the case,
c         snsqe may fail (usually with info = 4).  the user should
c         then either use snsq and set the step length or use iopt=1
c         and supply the jacobian.
c
c       lack of good progress.  snsqe searches for a zero of the system
c         by minimizing the sum of the squares of the functions.  in so
c         doing, it can become trapped in a region where the minimum
c         does not correspond to a zero of the system and, in this situ-
c         ation, the iteration eventually fails to make good progress.
c         in particular, this will happen if the system does not have a
c         zero.  if the system has a zero, rerunning snsqe from a dif-
c         ferent starting point may be helpful.
c
c
c 6. characteristics of the algorithm.
c
c       snsqe is a modification of the powell hybrid method.  two of
c       its main characteristics involve the choice of the correction as
c       a convex combination of the newton and scaled gradient direc-
c       tions, and the updating of the jacobian by the rank-1 method of
c       broyden.  the choice of the correction guarantees (under reason-
c       able conditions) global convergence for starting points far from
c       the solution and a fast rate of convergence.  the jacobian is
c       calculated at the starting point by either the user-supplied
c       subroutine or a forward-difference approximation, but it is not
c       recalculated until the rank-1 method fails to produce satis-
c       factory progress.
c
c       timing.  the time required by snsqe to solve a given problem
c         depends on n, the behavior of the functions, the accuracy
c         requested, and the starting point.  the number of arithmetic
c         operations needed by snsqe is about 11.5*(n**2) to process
c         each evaluation of the functions (call to fcn) and 1.3*(n**3)
c         to process each evaluation of the jacobian (call to jac,
c         if iopt = 1).  unless fcn and jac can be evaluated quickly,
c         the timing of snsqe will be strongly influenced by the time
c         spent in fcn and jac.
c
c       storage.  snsqe requires (3*n**2 + 17*n)/2 single precision
c         storage locations, in addition to the storage required by the
c         program.  there are no internally declared storage arrays.
c
c
c 7. example.
c
c       the problem is to determine the values of x(1), x(2), ..., x(9),
c       which solve the system of tridiagonal equations
c
c       (3-2*x(1))*x(1)           -2*x(2)                   = -1
c               -x(i-1) + (3-2*x(i))*x(i)         -2*x(i+1) = -1, i=2-8
c                                   -x(8) + (3-2*x(9))*x(9) = -1
c
c       **********
c
c       program test
c c
c c     driver for snsqe example.
c c
c       integer j,n,iopt,nprint,info,lwa,nwrite
c       real tol,fnorm
c       real x(9),fvec(9),wa(180)
c       real enorm,r1mach
c       external fcn
c       data nwrite /6/
c c
c       iopt = 2
c       n = 9
c c
c c     the following starting values provide a rough solution.
c c
c       do 10 j = 1, 9
c          x(j) = -1.e0
c    10    continue
c
c       lwa = 180
c       nprint = 0
c c
c c     set tol to the square root of the machine precision.
c c     unless high precision solutions are required,
c c     this is the recommended setting.
c c
c       tol = sqrt(r1mach(4))
c c
c       call snsqe(fcn,jac,iopt,n,x,fvec,tol,nprint,info,wa,lwa)
c       fnorm = enorm(n,fvec)
c       write (nwrite,1000) fnorm,info,(x(j),j=1,n)
c       stop
c  1000 format (5x,' final l2 norm of the residuals',e15.7 //
c      *        5x,' exit parameter',16x,i10 //
c      *        5x,' final approximate solution' // (5x,3e15.7))
c       end
c       subroutine fcn(n,x,fvec,iflag)
c       integer n,iflag
c       real x(n),fvec(n)
c       integer k
c       real one,temp,temp1,temp2,three,two,zero
c       data zero,one,two,three /0.e0,1.e0,2.e0,3.e0/
c c
c       do 10 k = 1, n
c          temp = (three - two*x(k))*x(k)
c          temp1 = zero
c          if (k .ne. 1) temp1 = x(k-1)
c          temp2 = zero
c          if (k .ne. n) temp2 = x(k+1)
c          fvec(k) = temp - temp1 - two*temp2 + one
c    10    continue
c       return
c       end
c
c       results obtained with different compilers or machines
c       may be slightly different.
c
c       final l2 norm of the residuals  0.1192636e-07
c
c       exit parameter                         1
c
c       final approximate solution
c
c       -0.5706545e+00 -0.6816283e+00 -0.7017325e+00
c       -0.7042129e+00 -0.7013690e+00 -0.6918656e+00
c       -0.6657920e+00 -0.5960342e+00 -0.4164121e+00
c
c***references  m. j. d. powell, a hybrid method for nonlinear equa-
c                 tions. in numerical methods for nonlinear algebraic
c                 equations, p. rabinowitz, editor.  gordon and breach,
c                 1988.
c***routines called  snsq, xermsg
c***revision history  (yymmdd)
c   800301  date written
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   920501  reformatted the references section.  (wrb)
c***end prologue  snsqe
      integer iopt,n,nprint,info,lwa
      real tol
      real x(*),fvec(*),wa(lwa)
      external fcn, jac
      integer index,j,lr,maxfev,ml,mode,mu,nfev,njev
      real epsfcn,factor,one,xtol,zero
      save factor, one, zero
      data factor,one,zero /1.0e2,1.0e0,0.0e0/
c***first executable statement  snsqe
      info = 0
c
c     check the input parameters for errors.
c
      if (iopt .lt. 1 .or. iopt .gt. 2 .or. n .le. 0
     1    .or. tol .lt. zero .or. lwa .lt. (3*n**2 +13*n)/2)
     2   go to 20
c
c     call snsq.
c
      maxfev = 100*(n + 1)
      if (iopt .eq. 2) maxfev = 2*maxfev
      xtol = tol
      ml = n - 1
      mu = n - 1
      epsfcn = zero
      mode = 2
      do 10 j = 1, n
         wa(j) = one
   10    continue
      lr = (n*(n + 1))/2
      index=6*n+lr
      call snsq(fcn,jac,iopt,n,x,fvec,wa(index+1),n,xtol,maxfev,ml,mu,
     1           epsfcn,wa(1),mode,factor,nprint,info,nfev,njev,
     2           wa(6*n+1),lr,wa(n+1),wa(2*n+1),wa(3*n+1),wa(4*n+1),
     3           wa(5*n+1))
      if (info .eq. 5) info = 4
   20 continue
      if (info .eq. 0) call xermsg ('slatec', 'snsqe',
     +   'invalid input parameter.', 2, 1)
      return
c
c     last card of subroutine snsqe.
c
      end
