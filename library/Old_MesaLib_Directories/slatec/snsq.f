*deck snsq
      subroutine snsq (fcn, jac, iopt, n, x, fvec, fjac, ldfjac, xtol,
     +   maxfev, ml, mu, epsfcn, diag, mode, factor, nprint, info, nfev,
     +   njev, r, lr, qtf, wa1, wa2, wa3, wa4)
c***begin prologue  snsq
c***purpose  find a zero of a system of a n nonlinear functions in n
c            variables by a modification of the powell hybrid method.
c***library   slatec
c***category  f2a
c***type      single precision (snsq-s, dnsq-d)
c***keywords  nonlinear square system, powell hybrid method, zeros
c***author  hiebert, k. l., (snla)
c***description
c
c 1. purpose.
c
c       the purpose of snsq is to find a zero of a system of n non-
c       linear functions in n variables by a modification of the powell
c       hybrid method.  the user must provide a subroutine which calcu-
c       lates the functions.  the user has the option of either to
c       provide a subroutine which calculates the jacobian or to let the
c       code calculate it by a forward-difference approximation.
c       this code is the combination of the minpack codes (argonne)
c       hybrd and hybrdj.
c
c
c 2. subroutine and type statements.
c
c       subroutine snsq(fcn,jac,iopt,n,x,fvec,fjac,ldfjac,xtol,maxfev,
c      *                 ml,mu,epsfcn,diag,mode,factor,nprint,info,nfev,
c      *                 njev,r,lr,qtf,wa1,wa2,wa3,wa4)
c       integer iopt,n,maxfev,ml,mu,mode,nprint,info,nfev,ldfjac,njev,lr
c       real xtol,epsfcn,factor
c       real x(n),fvec(n),diag(n),fjac(ldfjac,n),r(lr),qtf(n),
c      *     wa1(n),wa2(n),wa3(n),wa4(n)
c       external fcn,jac
c
c
c 3. parameters.
c
c       parameters designated as input parameters must be specified on
c       entry to snsq and are not changed on exit, while parameters
c       designated as output parameters need not be specified on entry
c       and are set to appropriate values on exit from snsq.
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
c         user wants to terminate execution of snsq.  in this case, set
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
c         user wants to terminate execution of snsq.  in this case, set
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
c       fjac is an output n by n array which contains the orthogonal
c         matrix q produced by the qr factorization of the final approx-
c         imate jacobian.
c
c       ldfjac is a positive integer input variable not less than n
c         which specifies the leading dimension of the array fjac.
c
c       xtol is a non-negative input variable.  termination occurs when
c         the relative error between two consecutive iterates is at most
c         xtol.  therefore, xtol measures the relative error desired in
c         the approximate solution.  section 4 contains more details
c         about xtol.
c
c       maxfev is a positive integer input variable.  termination occurs
c         when the number of calls to fcn is at least maxfev by the end
c         of an iteration.
c
c       ml is a non-negative integer input variable which specifies the
c         number of subdiagonals within the band of the jacobian matrix.
c         if the jacobian is not banded or iopt=1, set ml to at
c         least n - 1.
c
c       mu is a non-negative integer input variable which specifies the
c         number of superdiagonals within the band of the jacobian
c         matrix.  if the jacobian is not banded or iopt=1, set mu to at
c         least n - 1.
c
c       epsfcn is an input variable used in determining a suitable step
c         for the forward-difference approximation.  this approximation
c         assumes that the relative errors in the functions are of the
c         order of epsfcn.  if epsfcn is less than the machine preci-
c         sion, it is assumed that the relative errors in the functions
c         are of the order of the machine precision.  if iopt=1, then
c         epsfcn can be ignored (treat it as a dummy argument).
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
c         and every nprint iteration thereafter and immediately prior
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
c                   by the improvement from the last five jacobian eval-
c                   uations.
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
c       the accuracy of snsq is controlled by the convergence parameter
c       xtol.  this parameter is used in a test which makes a comparison
c       between the approximation x and a solution xsol.  snsq termi-
c       nates when the test is satisfied.  if the convergence parameter
c       is less than the machine precision (as defined by the function
c       r1mach(4)), then snsq only attempts to satisfy the test
c       defined by the machine precision.  further progress is not
c       usually possible.
c
c       the test assumes that the functions are reasonably well behaved,
c       and, if the jacobian is supplied by the user, that the functions
c       and the jacobian are coded consistently.  if these conditions
c       are not satisfied, then snsq may incorrectly indicate conver-
c       gence.  the coding of the jacobian can be checked by the
c       subroutine chkder. if the jacobian is coded correctly or iopt=2,
c       then the validity of the answer can be checked, for example, by
c       rerunning snsq with a tighter tolerance.
c
c       convergence test.  if enorm(z) denotes the euclidean norm of a
c         vector z and d is the diagonal matrix whose entries are
c         defined by the array diag, then this test attempts to guaran-
c         tee that
c
c               enorm(d*(x-xsol)) .le. xtol*enorm(d*xsol).
c
c         if this condition is satisfied with xtol = 10**(-k), then the
c         larger components of d*x have k significant decimal digits and
c         info is set to 1.  there is a danger that the smaller compo-
c         nents of d*x may have large relative errors, but the fast rate
c         of convergence of snsq usually avoids this possibility.
c         unless high precision solutions are required, the recommended
c         value for xtol is the square root of the machine precision.
c
c
c 5. unsuccessful completion.
c
c       unsuccessful termination of snsq can be due to improper input
c       parameters, arithmetic interrupts, an excessive number of func-
c       tion evaluations, or lack of good progress.
c
c       improper input parameters.  info is set to 0 if iopt .lt. 1,
c         or iopt .gt. 2, or n .le. 0, or ldfjac .lt. n, or
c         xtol .lt. 0.e0, or maxfev .le. 0, or ml .lt. 0, or mu .lt. 0,
c         or factor .le. 0.e0, or lr .lt. (n*(n+1))/2.
c
c       arithmetic interrupts.  if these interrupts occur in the fcn
c         subroutine during an early stage of the computation, they may
c         be caused by an unacceptable choice of x by snsq.  in this
c         case, it may be possible to remedy the situation by rerunning
c         snsq with a smaller value of factor.
c
c       excessive number of function evaluations.  a reasonable value
c         for maxfev is 100*(n+1) for iopt=1 and 200*(n+1) for iopt=2.
c         if the number of calls to fcn reaches maxfev, then this
c         indicates that the routine is converging very slowly as
c         measured by the progress of fvec, and info is set to 2.  this
c         situation should be unusual because, as indicated below, lack
c         of good progress is usually diagnosed earlier by snsq,
c         causing termination with info = 4 or info = 5.
c
c       lack of good progress.  snsq searches for a zero of the system
c         by minimizing the sum of the squares of the functions.  in so
c         doing, it can become trapped in a region where the minimum
c         does not correspond to a zero of the system and, in this situ-
c         ation, the iteration eventually fails to make good progress.
c         in particular, this will happen if the system does not have a
c         zero.  if the system has a zero, rerunning snsq from a dif-
c         ferent starting point may be helpful.
c
c
c 6. characteristics of the algorithm.
c
c       snsq is a modification of the powell hybrid method.  two of its
c       main characteristics involve the choice of the correction as a
c       convex combination of the newton and scaled gradient directions,
c       and the updating of the jacobian by the rank-1 method of broy-
c       den.  the choice of the correction guarantees (under reasonable
c       conditions) global convergence for starting points far from the
c       solution and a fast rate of convergence.  the jacobian is
c       calculated at the starting point by either the user-supplied
c       subroutine or a forward-difference approximation, but it is not
c       recalculated until the rank-1 method fails to produce satis-
c       factory progress.
c
c       timing.  the time required by snsq to solve a given problem
c         depends on n, the behavior of the functions, the accuracy
c         requested, and the starting point.  the number of arithmetic
c         operations needed by snsq is about 11.5*(n**2) to process
c         each evaluation of the functions (call to fcn) and 1.3*(n**3)
c         to process each evaluation of the jacobian (call to jac,
c         if iopt = 1).  unless fcn and jac can be evaluated quickly,
c         the timing of snsq will be strongly influenced by the time
c         spent in fcn and jac.
c
c       storage.  snsq requires (3*n**2 + 17*n)/2 single precision
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
c c     **********
c
c       program test
c c
c c     driver for snsq example.
c c
c       integer j,iopt,n,maxfev,ml,mu,mode,nprint,info,nfev,ldfjac,lr,
c      *        nwrite
c       real xtol,epsfcn,factor,fnorm
c       real x(9),fvec(9),diag(9),fjac(9,9),r(45),qtf(9),
c      *     wa1(9),wa2(9),wa3(9),wa4(9)
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
c c
c       ldfjac = 9
c       lr = 45
c c
c c     set xtol to the square root of the machine precision.
c c     unless high precision solutions are required,
c c     this is the recommended setting.
c c
c       xtol = sqrt(r1mach(4))
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
c       call snsq(fcn,jac,iopt,n,x,fvec,fjac,ldfjac,xtol,maxfev,ml,mu,
c      *           epsfcn,diag,mode,factor,nprint,info,nfev,njev,
c      *           r,lr,qtf,wa1,wa2,wa3,wa4)
c       fnorm = enorm(n,fvec)
c       write (nwrite,1000) fnorm,nfev,info,(x(j),j=1,n)
c       stop
c  1000 format (5x,' final l2 norm of the residuals',e15.7 //
c      *        5x,' number of function evaluations',i10 //
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
c***routines called  dogleg, enorm, fdjac1, qform, qrfac, r1mach,
c                    r1mpyq, r1updt, xermsg
c***revision history  (yymmdd)
c   800301  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   920501  reformatted the references section.  (wrb)
c***end prologue  snsq
      integer iopt,n,maxfev,ml,mu,mode,nprint,info,nfev,ldfjac,lr,njev
      real xtol,epsfcn,factor
      real x(*),fvec(*),diag(*),fjac(ldfjac,*),r(lr),qtf(*),wa1(*),
     1     wa2(*),wa3(*),wa4(*)
      external fcn
      integer i,iflag,iter,j,jm1,l,ncfail,ncsuc,nslow1,nslow2
      integer iwa(1)
      logical jeval,sing
      real actred,delta,epsmch,fnorm,fnorm1,one,pnorm,prered,p1,p5,
     1     p001,p0001,ratio,sum,temp,xnorm,zero
      real r1mach,enorm
      save one, p1, p5, p001, p0001, zero
      data one,p1,p5,p001,p0001,zero
     1     /1.0e0,1.0e-1,5.0e-1,1.0e-3,1.0e-4,0.0e0/
c
c***first executable statement  snsq
      epsmch = r1mach(4)
c
      info = 0
      iflag = 0
      nfev = 0
      njev = 0
c
c     check the input parameters for errors.
c
      if (iopt .lt. 1 .or. iopt .gt. 2 .or.
     1    n .le. 0 .or. xtol .lt. zero .or. maxfev .le. 0
     2    .or. ml .lt. 0 .or. mu .lt. 0 .or. factor .le. zero
     3    .or. ldfjac .lt. n .or. lr .lt. (n*(n + 1))/2) go to 300
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
      call fcn(n,x,fvec,iflag)
      nfev = 1
      if (iflag .lt. 0) go to 300
      fnorm = enorm(n,fvec)
c
c     initialize iteration counter and monitors.
c
      iter = 1
      ncsuc = 0
      ncfail = 0
      nslow1 = 0
      nslow2 = 0
c
c     beginning of the outer loop.
c
   30 continue
         jeval = .true.
c
c        calculate the jacobian matrix.
c
         if (iopt .eq. 2) go to 31
c
c        user supplies jacobian
c
            call jac(n,x,fvec,fjac,ldfjac,iflag)
            njev = njev+1
            go to 32
c
c        code approximates the jacobian
c
   31       iflag = 2
            call fdjac1(fcn,n,x,fvec,fjac,ldfjac,iflag,ml,mu,epsfcn,wa1,
     1               wa2)
            nfev = nfev + min(ml+mu+1,n)
c
   32    if (iflag .lt. 0) go to 300
c
c        compute the qr factorization of the jacobian.
c
         call qrfac(n,n,fjac,ldfjac,.false.,iwa,1,wa1,wa2,wa3)
c
c        on the first iteration and if mode is 1, scale according
c        to the norms of the columns of the initial jacobian.
c
         if (iter .ne. 1) go to 70
         if (mode .eq. 2) go to 50
         do 40 j = 1, n
            diag(j) = wa2(j)
            if (wa2(j) .eq. zero) diag(j) = one
   40       continue
   50    continue
c
c        on the first iteration, calculate the norm of the scaled x
c        and initialize the step bound delta.
c
         do 60 j = 1, n
            wa3(j) = diag(j)*x(j)
   60       continue
         xnorm = enorm(n,wa3)
         delta = factor*xnorm
         if (delta .eq. zero) delta = factor
   70    continue
c
c        form (q transpose)*fvec and store in qtf.
c
         do 80 i = 1, n
            qtf(i) = fvec(i)
   80       continue
         do 120 j = 1, n
            if (fjac(j,j) .eq. zero) go to 110
            sum = zero
            do 90 i = j, n
               sum = sum + fjac(i,j)*qtf(i)
   90          continue
            temp = -sum/fjac(j,j)
            do 100 i = j, n
               qtf(i) = qtf(i) + fjac(i,j)*temp
  100          continue
  110       continue
  120       continue
c
c        copy the triangular factor of the qr factorization into r.
c
         sing = .false.
         do 150 j = 1, n
            l = j
            jm1 = j - 1
            if (jm1 .lt. 1) go to 140
            do 130 i = 1, jm1
               r(l) = fjac(i,j)
               l = l + n - i
  130          continue
  140       continue
            r(l) = wa1(j)
            if (wa1(j) .eq. zero) sing = .true.
  150       continue
c
c        accumulate the orthogonal factor in fjac.
c
         call qform(n,n,fjac,ldfjac,wa1)
c
c        rescale if necessary.
c
         if (mode .eq. 2) go to 170
         do 160 j = 1, n
            diag(j) = max(diag(j),wa2(j))
  160       continue
  170    continue
c
c        beginning of the inner loop.
c
  180    continue
c
c           if requested, call fcn to enable printing of iterates.
c
            if (nprint .le. 0) go to 190
            iflag = 0
            if (mod(iter-1,nprint) .eq. 0) call fcn(n,x,fvec,iflag)
            if (iflag .lt. 0) go to 300
  190       continue
c
c           determine the direction p.
c
            call dogleg(n,r,lr,diag,qtf,delta,wa1,wa2,wa3)
c
c           store the direction p and x + p. calculate the norm of p.
c
            do 200 j = 1, n
               wa1(j) = -wa1(j)
               wa2(j) = x(j) + wa1(j)
               wa3(j) = diag(j)*wa1(j)
  200          continue
            pnorm = enorm(n,wa3)
c
c           on the first iteration, adjust the initial step bound.
c
            if (iter .eq. 1) delta = min(delta,pnorm)
c
c           evaluate the function at x + p and calculate its norm.
c
            iflag = 1
            call fcn(n,wa2,wa4,iflag)
            nfev = nfev + 1
            if (iflag .lt. 0) go to 300
            fnorm1 = enorm(n,wa4)
c
c           compute the scaled actual reduction.
c
            actred = -one
            if (fnorm1 .lt. fnorm) actred = one - (fnorm1/fnorm)**2
c
c           compute the scaled predicted reduction.
c
            l = 1
            do 220 i = 1, n
               sum = zero
               do 210 j = i, n
                  sum = sum + r(l)*wa1(j)
                  l = l + 1
  210             continue
               wa3(i) = qtf(i) + sum
  220          continue
            temp = enorm(n,wa3)
            prered = zero
            if (temp .lt. fnorm) prered = one - (temp/fnorm)**2
c
c           compute the ratio of the actual to the predicted
c           reduction.
c
            ratio = zero
            if (prered .gt. zero) ratio = actred/prered
c
c           update the step bound.
c
            if (ratio .ge. p1) go to 230
               ncsuc = 0
               ncfail = ncfail + 1
               delta = p5*delta
               go to 240
  230       continue
               ncfail = 0
               ncsuc = ncsuc + 1
               if (ratio .ge. p5 .or. ncsuc .gt. 1)
     1            delta = max(delta,pnorm/p5)
               if (abs(ratio-one) .le. p1) delta = pnorm/p5
  240       continue
c
c           test for successful iteration.
c
            if (ratio .lt. p0001) go to 260
c
c           successful iteration. update x, fvec, and their norms.
c
            do 250 j = 1, n
               x(j) = wa2(j)
               wa2(j) = diag(j)*x(j)
               fvec(j) = wa4(j)
  250          continue
            xnorm = enorm(n,wa2)
            fnorm = fnorm1
            iter = iter + 1
  260       continue
c
c           determine the progress of the iteration.
c
            nslow1 = nslow1 + 1
            if (actred .ge. p001) nslow1 = 0
            if (jeval) nslow2 = nslow2 + 1
            if (actred .ge. p1) nslow2 = 0
c
c           test for convergence.
c
            if (delta .le. xtol*xnorm .or. fnorm .eq. zero) info = 1
            if (info .ne. 0) go to 300
c
c           tests for termination and stringent tolerances.
c
            if (nfev .ge. maxfev) info = 2
            if (p1*max(p1*delta,pnorm) .le. epsmch*xnorm) info = 3
            if (nslow2 .eq. 5) info = 4
            if (nslow1 .eq. 10) info = 5
            if (info .ne. 0) go to 300
c
c           criterion for recalculating jacobian
c
            if (ncfail .eq. 2) go to 290
c
c           calculate the rank one modification to the jacobian
c           and update qtf if necessary.
c
            do 280 j = 1, n
               sum = zero
               do 270 i = 1, n
                  sum = sum + fjac(i,j)*wa4(i)
  270             continue
               wa2(j) = (sum - wa3(j))/pnorm
               wa1(j) = diag(j)*((diag(j)*wa1(j))/pnorm)
               if (ratio .ge. p0001) qtf(j) = sum
  280          continue
c
c           compute the qr factorization of the updated jacobian.
c
            call r1updt(n,n,r,lr,wa1,wa2,wa3,sing)
            call r1mpyq(n,n,fjac,ldfjac,wa2,wa3)
            call r1mpyq(1,n,qtf,1,wa2,wa3)
c
c           end of the inner loop.
c
            jeval = .false.
            go to 180
  290    continue
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
      if (nprint .gt. 0) call fcn(n,x,fvec,iflag)
      if (info .lt. 0) call xermsg ('slatec', 'snsq',
     +   'execution terminated because user set iflag negative.', 1, 1)
      if (info .eq. 0) call xermsg ('slatec', 'snsq',
     +   'invalid input parameter.', 2, 1)
      if (info .eq. 2) call xermsg ('slatec', 'snsq',
     +   'too many function evaluations.', 9, 1)
      if (info .eq. 3) call xermsg ('slatec', 'snsq',
     +   'xtol too small. no further improvement possible.', 3, 1)
      if (info .gt. 4) call xermsg ('slatec', 'snsq',
     +   'iteration not making good progress.', 1, 1)
      return
c
c     last card of subroutine snsq.
c
      end
