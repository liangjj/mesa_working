*deck dnsq
      subroutine dnsq (fcn, jac, iopt, n, x, fvec, fjac, ldfjac, xtol,
     +   maxfev, ml, mu, epsfcn, diag, mode, factor, nprint, info, nfev,
     +   njev, r, lr, qtf, wa1, wa2, wa3, wa4)
c***begin prologue  dnsq
c***purpose  find a zero of a system of a n nonlinear functions in n
c            variables by a modification of the powell hybrid method.
c***library   slatec
c***category  f2a
c***type      double precision (snsq-s, dnsq-d)
c***keywords  nonlinear square system, powell hybrid method, zeros
c***author  hiebert, k. l. (snla)
c***description
c
c 1. purpose.
c
c       the purpose of dnsq is to find a zero of a system of n nonlinear
c       functions in n variables by a modification of the powell
c       hybrid method.  the user must provide a subroutine which
c       calculates the functions.  the user has the option of either to
c       provide a subroutine which calculates the jacobian or to let the
c       code calculate it by a forward-difference approximation.
c       this code is the combination of the minpack codes (argonne)
c       hybrd and hybrdj.
c
c 2. subroutine and type statements.
c
c       subroutine dnsq(fcn,jac,iopt,n,x,fvec,fjac,ldfjac,xtol,maxfev,
c      *                 ml,mu,epsfcn,diag,mode,factor,nprint,info,nfev,
c      *                 njev,r,lr,qtf,wa1,wa2,wa3,wa4)
c       integer iopt,n,maxfev,ml,mu,mode,nprint,info,nfev,ldfjac,njev,lr
c       double precision xtol,epsfcn,factor
c       double precision
c       x(n),fvec(n),diag(n),fjac(ldfjac,n),r(lr),qtf(n),
c      *     wa1(n),wa2(n),wa3(n),wa4(n)
c       external fcn,jac
c
c 3. parameters.
c
c       parameters designated as input parameters must be specified on
c       entry to dnsq and are not changed on exit, while parameters
c       designated as output parameters need not be specified on entry
c       and are set to appropriate values on exit from dnsq.
c
c       fcn is the name of the user-supplied subroutine which calculates
c         the functions.  fcn must be declared in an external statement
c         in the user calling program, and should be written as follows.
c
c         subroutine fcn(n,x,fvec,iflag)
c         integer n,iflag
c         double precision x(n),fvec(n)
c         ----------
c         calculate the functions at x and
c         return this vector in fvec.
c         ----------
c         return
c         end
c
c         the value of iflag should not be changed by fcn unless the
c         user wants to terminate execution of dnsq.  in this case set
c         iflag to a negative integer.
c
c       jac is the name of the user-supplied subroutine which calculates
c         the jacobian.  if iopt=1, then jac must be declared in an
c         external statement in the user calling program, and should be
c         written as follows.
c
c         subroutine jac(n,x,fvec,fjac,ldfjac,iflag)
c         integer n,ldfjac,iflag
c         double precision x(n),fvec(n),fjac(ldfjac,n)
c         ----------
c         calculate the jacobian at x and return this
c         matrix in fjac.  fvec contains the function
c         values at x and should not be altered.
c         ----------
c         return
c         end
c
c         the value of iflag should not be changed by jac unless the
c         user wants to terminate execution of dnsq.  in this case set
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
c       x is an array of length n.  on input x must contain an initial
c         estimate of the solution vector.  on output x contains the
c         final estimate of the solution vector.
c
c       fvec is an output array of length n which contains the functions
c         evaluated at the output x.
c
c       fjac is an output n by n array which contains the orthogonal
c         matrix q produced by the qr factorization of the final
c         approximate jacobian.
c
c       ldfjac is a positive integer input variable not less than n
c         which specifies the leading dimension of the array fjac.
c
c       xtol is a nonnegative input variable.  termination occurs when
c         the relative error between two consecutive iterates is at most
c         xtol.  therefore, xtol measures the relative error desired in
c         the approximate solution.  section 4 contains more details
c         about xtol.
c
c       maxfev is a positive integer input variable.  termination occurs
c         when the number of calls to fcn is at least maxfev by the end
c         of an iteration.
c
c       ml is a nonnegative integer input variable which specifies the
c         number of subdiagonals within the band of the jacobian matrix.
c         if the jacobian is not banded or iopt=1, set ml to at
c         least n - 1.
c
c       mu is a nonnegative integer input variable which specifies the
c         number of superdiagonals within the band of the jacobian
c         matrix.  if the jacobian is not banded or iopt=1, set mu to at
c         least n - 1.
c
c       epsfcn is an input variable used in determining a suitable step
c         for the forward-difference approximation.  this approximation
c         assumes that the relative errors in the functions are of the
c         order of epsfcn.  if epsfcn is less than the machine
c         precision, it is assumed that the relative errors in the
c         functions are of the order of the machine precision.  if
c         iopt=1, then epsfcn can be ignored (treat it as a dummy
c         argument).
c
c       diag is an array of length n.  if mode = 1 (see below), diag is
c         internally set.  if mode = 2, diag must contain positive
c         entries that serve as implicit (multiplicative) scale factors
c         for the variables.
c
c       mode is an integer input variable.  if mode = 1, the variables
c         will be scaled internally.  if mode = 2, the scaling is
c         specified by the input diag.  other values of mode are
c         equivalent to mode = 1.
c
c       factor is a positive input variable used in determining the
c         initial step bound.  this bound is set to the product of
c         factor and the euclidean norm of diag*x if nonzero, or else to
c         factor itself.  in most cases factor should lie in the
c         interval (.1,100.).  100. is a generally recommended value.
c
c       nprint is an integer input variable that enables controlled
c         printing of iterates if it is positive.  in this case, fcn is
c         called with iflag = 0 at the beginning of the first iteration
c         and every nprint iterations thereafter and immediately prior
c         to return, with x and fvec available for printing. appropriate
c         print statements must be added to fcn(see example).  if nprint
c         is not positive, no special calls of fcn with iflag = 0 are
c         made.
c
c       info is an integer output variable.  if the user has terminated
c         execution, info is set to the (negative) value of iflag.  see
c         description of fcn and jac. otherwise, info is set as follows.
c
c         info = 0  improper input parameters.
c
c         info = 1  relative error between two consecutive iterates is
c                   at most xtol.
c
c         info = 2  number of calls to fcn has reached or exceeded
c                   maxfev.
c
c         info = 3  xtol is too small.  no further improvement in the
c                   approximate solution x is possible.
c
c         info = 4  iteration is not making good progress, as measured
c                   by the improvement from the last five jacobian
c                   evaluations.
c
c         info = 5  iteration is not making good progress, as measured
c                   by the improvement from the last ten iterations.
c
c         sections 4 and 5 contain more details about info.
c
c       nfev is an integer output variable set to the number of calls to
c         fcn.
c
c       njev is an integer output variable set to the number of calls to
c         jac. (if iopt=2, then njev is set to zero.)
c
c       r is an output array of length lr which contains the upper
c         triangular matrix produced by the qr factorization of the
c         final approximate jacobian, stored rowwise.
c
c       lr is a positive integer input variable not less than
c         (n*(n+1))/2.
c
c       qtf is an output array of length n which contains the vector
c         (q transpose)*fvec.
c
c       wa1, wa2, wa3, and wa4 are work arrays of length n.
c
c
c 4. successful completion.
c
c       the accuracy of dnsq is controlled by the convergence parameter
c       xtol.  this parameter is used in a test which makes a comparison
c       between the approximation x and a solution xsol.  dnsq
c       terminates when the test is satisfied.  if the convergence
c       parameter is less than the machine precision (as defined by the
c       function d1mach(4)), then dnsq only attempts to satisfy the test
c       defined by the machine precision.  further progress is not
c       usually possible.
c
c       the test assumes that the functions are reasonably well behaved,
c       and, if the jacobian is supplied by the user, that the functions
c       and the jacobian are coded consistently.  if these conditions
c       are not satisfied, then dnsq may incorrectly indicate
c       convergence.  the coding of the jacobian can be checked by the
c       subroutine dckder. if the jacobian is coded correctly or iopt=2,
c       then the validity of the answer can be checked, for example, by
c       rerunning dnsq with a tighter tolerance.
c
c       convergence test.  if denorm(z) denotes the euclidean norm of a
c         vector z and d is the diagonal matrix whose entries are
c         defined by the array diag, then this test attempts to
c         guarantee that
c
c               denorm(d*(x-xsol)) .le. xtol*denorm(d*xsol).
c
c         if this condition is satisfied with xtol = 10**(-k), then the
c         larger components of d*x have k significant decimal digits and
c         info is set to 1.  there is a danger that the smaller
c         components of d*x may have large relative errors, but the fast
c         rate of convergence of dnsq usually avoids this possibility.
c         unless high precision solutions are required, the recommended
c         value for xtol is the square root of the machine precision.
c
c
c 5. unsuccessful completion.
c
c       unsuccessful termination of dnsq can be due to improper input
c       parameters, arithmetic interrupts, an excessive number of
c       function evaluations, or lack of good progress.
c
c       improper input parameters.  info is set to 0 if iopt .lt .1,
c         or iopt .gt. 2, or n .le. 0, or ldfjac .lt. n, or
c         xtol .lt. 0.e0, or maxfev .le. 0, or ml .lt. 0, or mu .lt. 0,
c         or factor .le. 0.e0, or lr .lt. (n*(n+1))/2.
c
c       arithmetic interrupts.  if these interrupts occur in the fcn
c         subroutine during an early stage of the computation, they may
c         be caused by an unacceptable choice of x by dnsq.  in this
c         case, it may be possible to remedy the situation by rerunning
c         dnsq with a smaller value of factor.
c
c       excessive number of function evaluations.  a reasonable value
c         for maxfev is 100*(n+1) for iopt=1 and 200*(n+1) for iopt=2.
c         if the number of calls to fcn reaches maxfev, then this
c         indicates that the routine is converging very slowly as
c         measured by the progress of fvec, and info is set to 2. this
c         situation should be unusual because, as indicated below, lack
c         of good progress is usually diagnosed earlier by dnsq,
c         causing termination with info = 4 or info = 5.
c
c       lack of good progress.  dnsq searches for a zero of the system
c         by minimizing the sum of the squares of the functions.  in so
c         doing, it can become trapped in a region where the minimum
c         does not correspond to a zero of the system and, in this
c         situation, the iteration eventually fails to make good
c         progress.  in particular, this will happen if the system does
c         not have a zero.  if the system has a zero, rerunning dnsq
c         from a different starting point may be helpful.
c
c
c 6. characteristics of the algorithm.
c
c       dnsq is a modification of the powell hybrid method.  two of its
c       main characteristics involve the choice of the correction as a
c       convex combination of the newton and scaled gradient directions,
c       and the updating of the jacobian by the rank-1 method of
c       broyden.  the choice of the correction guarantees (under
c       reasonable conditions) global convergence for starting points
c       far from the solution and a fast rate of convergence.  the
c       jacobian is calculated at the starting point by either the
c       user-supplied subroutine or a forward-difference approximation,
c       but it is not recalculated until the rank-1 method fails to
c       produce satisfactory progress.
c
c       timing.  the time required by dnsq to solve a given problem
c         depends on n, the behavior of the functions, the accuracy
c         requested, and the starting point.  the number of arithmetic
c         operations needed by dnsq is about 11.5*(n**2) to process
c         each evaluation of the functions (call to fcn) and 1.3*(n**3)
c         to process each evaluation of the jacobian (call to jac,
c         if iopt = 1).  unless fcn and jac can be evaluated quickly,
c         the timing of dnsq will be strongly influenced by the time
c         spent in fcn and jac.
c
c       storage.  dnsq requires (3*n**2 + 17*n)/2 single precision
c         storage locations, in addition to the storage required by the
c         program.  there are no internally declared storage arrays.
c
c *long description:
c
c 7. example.
c
c       the problem is to determine the values of x(1), x(2), ..., x(9),
c       which solve the system of tridiagonal equations
c
c       (3-2*x(1))*x(1)           -2*x(2)                   = -1
c               -x(i-1) + (3-2*x(i))*x(i)         -2*x(i+1) = -1, i=2-8
c                                   -x(8) + (3-2*x(9))*x(9) = -1
c c     **********
c
c       program test
c c
c c     driver for dnsq example.
c c
c       integer j,iopt,n,maxfev,ml,mu,mode,nprint,info,nfev,ldfjac,lr,
c      *        nwrite
c       double precision xtol,epsfcn,factor,fnorm
c       double precision x(9),fvec(9),diag(9),fjac(9,9),r(45),qtf(9),
c      *     wa1(9),wa2(9),wa3(9),wa4(9)
c       double precision denorm,d1mach
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
c c
c       ldfjac = 9
c       lr = 45
c c
c c     set xtol to the square root of the machine precision.
c c     unless high precision solutions are required,
c c     this is the recommended setting.
c c
c       xtol = sqrt(d1mach(4))
c c
c       maxfev = 2000
c       ml = 1
c       mu = 1
c       epsfcn = 0.e0
c       mode = 2
c       do 20 j = 1, 9
c          diag(j) = 1.e0
c    20    continue
c       factor = 1.e2
c       nprint = 0
c c
c       call dnsq(fcn,jac,iopt,n,x,fvec,fjac,ldfjac,xtol,maxfev,ml,mu,
c      *           epsfcn,diag,mode,factor,nprint,info,nfev,njev,
c      *           r,lr,qtf,wa1,wa2,wa3,wa4)
c       fnorm = denorm(n,fvec)
c       write (nwrite,1000) fnorm,nfev,info,(x(j),j=1,n)
c       stop
c  1000 format (5x,' final l2 norm of the residuals',e15.7 //
c      *        5x,' number of function evaluations',i10 //
c      *        5x,' exit parameter',16x,i10 //
c      *        5x,' final approximate solution' // (5x,3e15.7))
c       end
c       subroutine fcn(n,x,fvec,iflag)
c       integer n,iflag
c       double precision x(n),fvec(n)
c       integer k
c       double precision one,temp,temp1,temp2,three,two,zero
c       data zero,one,two,three /0.e0,1.e0,2.e0,3.e0/
c c
c       if (iflag .ne. 0) go to 5
c c
c c     insert print statements here when nprint is positive.
c c
c       return
c     5 continue
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
c       number of function evaluations        14
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
c***routines called  d1mach, d1mpyq, d1updt, ddoglg, denorm, dfdjc1,
c                    dqform, dqrfac, xermsg
c***revision history  (yymmdd)
c   800301  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   920501  reformatted the references section.  (wrb)
c***end prologue  dnsq
      double precision d1mach,denorm
      integer i, iflag, info, iopt, iter, iwa(1), j, jm1, l, ldfjac,
     1     lr, maxfev, ml, mode, mu, n, ncfail, ncsuc, nfev, njev,
     2     nprint, nslow1, nslow2
      double precision actred, delta, diag(*), epsfcn, epsmch, factor,
     1     fjac(ldfjac,*), fnorm, fnorm1, fvec(*), one, p0001, p001,
     2     p1, p5, pnorm, prered, qtf(*), r(*), ratio, sum, temp,
     3     wa1(*), wa2(*), wa3(*), wa4(*), x(*), xnorm, xtol, zero
      external fcn
      logical jeval,sing
      save one, p1, p5, p001, p0001, zero
      data one,p1,p5,p001,p0001,zero
     1     /1.0d0,1.0d-1,5.0d-1,1.0d-3,1.0d-4,0.0d0/
c
c     begin block permitting ...exits to 320
c***first executable statement  dnsq
         epsmch = d1mach(4)
c
         info = 0
         iflag = 0
         nfev = 0
         njev = 0
c
c        check the input parameters for errors.
c
c     ...exit
         if (iopt .lt. 1 .or. iopt .gt. 2 .or. n .le. 0
     1       .or. xtol .lt. zero .or. maxfev .le. 0 .or. ml .lt. 0
     2       .or. mu .lt. 0 .or. factor .le. zero .or. ldfjac .lt. n
     3       .or. lr .lt. (n*(n + 1))/2) go to 320
         if (mode .ne. 2) go to 20
            do 10 j = 1, n
c     .........exit
               if (diag(j) .le. zero) go to 320
   10       continue
   20    continue
c
c        evaluate the function at the starting point
c        and calculate its norm.
c
         iflag = 1
         call fcn(n,x,fvec,iflag)
         nfev = 1
c     ...exit
         if (iflag .lt. 0) go to 320
         fnorm = denorm(n,fvec)
c
c        initialize iteration counter and monitors.
c
         iter = 1
         ncsuc = 0
         ncfail = 0
         nslow1 = 0
         nslow2 = 0
c
c        beginning of the outer loop.
c
   30    continue
c           begin block permitting ...exits to 90
               jeval = .true.
c
c              calculate the jacobian matrix.
c
               if (iopt .eq. 2) go to 40
c
c                 user supplies jacobian
c
                  call jac(n,x,fvec,fjac,ldfjac,iflag)
                  njev = njev + 1
               go to 50
   40          continue
c
c                 code approximates the jacobian
c
                  iflag = 2
                  call dfdjc1(fcn,n,x,fvec,fjac,ldfjac,iflag,ml,mu,
     1                        epsfcn,wa1,wa2)
                  nfev = nfev + min(ml+mu+1,n)
   50          continue
c
c     .........exit
               if (iflag .lt. 0) go to 320
c
c              compute the qr factorization of the jacobian.
c
               call dqrfac(n,n,fjac,ldfjac,.false.,iwa,1,wa1,wa2,wa3)
c
c              on the first iteration and if mode is 1, scale according
c              to the norms of the columns of the initial jacobian.
c
c           ...exit
               if (iter .ne. 1) go to 90
               if (mode .eq. 2) go to 70
                  do 60 j = 1, n
                     diag(j) = wa2(j)
                     if (wa2(j) .eq. zero) diag(j) = one
   60             continue
   70          continue
c
c              on the first iteration, calculate the norm of the scaled
c              x and initialize the step bound delta.
c
               do 80 j = 1, n
                  wa3(j) = diag(j)*x(j)
   80          continue
               xnorm = denorm(n,wa3)
               delta = factor*xnorm
               if (delta .eq. zero) delta = factor
   90       continue
c
c           form (q transpose)*fvec and store in qtf.
c
            do 100 i = 1, n
               qtf(i) = fvec(i)
  100       continue
            do 140 j = 1, n
               if (fjac(j,j) .eq. zero) go to 130
                  sum = zero
                  do 110 i = j, n
                     sum = sum + fjac(i,j)*qtf(i)
  110             continue
                  temp = -sum/fjac(j,j)
                  do 120 i = j, n
                     qtf(i) = qtf(i) + fjac(i,j)*temp
  120             continue
  130          continue
  140       continue
c
c           copy the triangular factor of the qr factorization into r.
c
            sing = .false.
            do 170 j = 1, n
               l = j
               jm1 = j - 1
               if (jm1 .lt. 1) go to 160
               do 150 i = 1, jm1
                  r(l) = fjac(i,j)
                  l = l + n - i
  150          continue
  160          continue
               r(l) = wa1(j)
               if (wa1(j) .eq. zero) sing = .true.
  170       continue
c
c           accumulate the orthogonal factor in fjac.
c
            call dqform(n,n,fjac,ldfjac,wa1)
c
c           rescale if necessary.
c
            if (mode .eq. 2) go to 190
               do 180 j = 1, n
                  diag(j) = max(diag(j),wa2(j))
  180          continue
  190       continue
c
c           beginning of the inner loop.
c
  200       continue
c
c              if requested, call fcn to enable printing of iterates.
c
               if (nprint .le. 0) go to 210
                  iflag = 0
                  if (mod(iter-1,nprint) .eq. 0)
     1               call fcn(n,x,fvec,iflag)
c     ............exit
                  if (iflag .lt. 0) go to 320
  210          continue
c
c              determine the direction p.
c
               call ddoglg(n,r,lr,diag,qtf,delta,wa1,wa2,wa3)
c
c              store the direction p and x + p. calculate the norm of p.
c
               do 220 j = 1, n
                  wa1(j) = -wa1(j)
                  wa2(j) = x(j) + wa1(j)
                  wa3(j) = diag(j)*wa1(j)
  220          continue
               pnorm = denorm(n,wa3)
c
c              on the first iteration, adjust the initial step bound.
c
               if (iter .eq. 1) delta = min(delta,pnorm)
c
c              evaluate the function at x + p and calculate its norm.
c
               iflag = 1
               call fcn(n,wa2,wa4,iflag)
               nfev = nfev + 1
c     .........exit
               if (iflag .lt. 0) go to 320
               fnorm1 = denorm(n,wa4)
c
c              compute the scaled actual reduction.
c
               actred = -one
               if (fnorm1 .lt. fnorm) actred = one - (fnorm1/fnorm)**2
c
c              compute the scaled predicted reduction.
c
               l = 1
               do 240 i = 1, n
                  sum = zero
                  do 230 j = i, n
                     sum = sum + r(l)*wa1(j)
                     l = l + 1
  230             continue
                  wa3(i) = qtf(i) + sum
  240          continue
               temp = denorm(n,wa3)
               prered = zero
               if (temp .lt. fnorm) prered = one - (temp/fnorm)**2
c
c              compute the ratio of the actual to the predicted
c              reduction.
c
               ratio = zero
               if (prered .gt. zero) ratio = actred/prered
c
c              update the step bound.
c
               if (ratio .ge. p1) go to 250
                  ncsuc = 0
                  ncfail = ncfail + 1
                  delta = p5*delta
               go to 260
  250          continue
                  ncfail = 0
                  ncsuc = ncsuc + 1
                  if (ratio .ge. p5 .or. ncsuc .gt. 1)
     1               delta = max(delta,pnorm/p5)
                  if (abs(ratio-one) .le. p1) delta = pnorm/p5
  260          continue
c
c              test for successful iteration.
c
               if (ratio .lt. p0001) go to 280
c
c                 successful iteration. update x, fvec, and their norms.
c
                  do 270 j = 1, n
                     x(j) = wa2(j)
                     wa2(j) = diag(j)*x(j)
                     fvec(j) = wa4(j)
  270             continue
                  xnorm = denorm(n,wa2)
                  fnorm = fnorm1
                  iter = iter + 1
  280          continue
c
c              determine the progress of the iteration.
c
               nslow1 = nslow1 + 1
               if (actred .ge. p001) nslow1 = 0
               if (jeval) nslow2 = nslow2 + 1
               if (actred .ge. p1) nslow2 = 0
c
c              test for convergence.
c
               if (delta .le. xtol*xnorm .or. fnorm .eq. zero) info = 1
c     .........exit
               if (info .ne. 0) go to 320
c
c              tests for termination and stringent tolerances.
c
               if (nfev .ge. maxfev) info = 2
               if (p1*max(p1*delta,pnorm) .le. epsmch*xnorm) info = 3
               if (nslow2 .eq. 5) info = 4
               if (nslow1 .eq. 10) info = 5
c     .........exit
               if (info .ne. 0) go to 320
c
c              criterion for recalculating jacobian
c
c           ...exit
               if (ncfail .eq. 2) go to 310
c
c              calculate the rank one modification to the jacobian
c              and update qtf if necessary.
c
               do 300 j = 1, n
                  sum = zero
                  do 290 i = 1, n
                     sum = sum + fjac(i,j)*wa4(i)
  290             continue
                  wa2(j) = (sum - wa3(j))/pnorm
                  wa1(j) = diag(j)*((diag(j)*wa1(j))/pnorm)
                  if (ratio .ge. p0001) qtf(j) = sum
  300          continue
c
c              compute the qr factorization of the updated jacobian.
c
               call d1updt(n,n,r,lr,wa1,wa2,wa3,sing)
               call d1mpyq(n,n,fjac,ldfjac,wa2,wa3)
               call d1mpyq(1,n,qtf,1,wa2,wa3)
c
c              end of the inner loop.
c
               jeval = .false.
            go to 200
  310       continue
c
c           end of the outer loop.
c
         go to 30
  320 continue
c
c     termination, either normal or user imposed.
c
      if (iflag .lt. 0) info = iflag
      iflag = 0
      if (nprint .gt. 0) call fcn(n,x,fvec,iflag)
      if (info .lt. 0) call xermsg ('slatec', 'dnsq',
     +   'execution terminated because user set iflag negative.', 1, 1)
      if (info .eq. 0) call xermsg ('slatec', 'dnsq',
     +   'invalid input parameter.', 2, 1)
      if (info .eq. 2) call xermsg ('slatec', 'dnsq',
     +   'too many function evaluations.', 9, 1)
      if (info .eq. 3) call xermsg ('slatec', 'dnsq',
     +   'xtol too small. no further improvement possible.', 3, 1)
      if (info .gt. 4) call xermsg ('slatec', 'dnsq',
     +   'iteration not making good progress.', 1, 1)
      return
c
c     last card of subroutine dnsq.
c
      end
