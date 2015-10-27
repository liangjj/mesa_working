*deck sos
      subroutine sos (fnc, neq, x, rtolx, atolx, tolf, iflag, rw, lrw,
     +   iw, liw)
c***begin prologue  sos
c***purpose  solve a square system of nonlinear equations.
c***library   slatec
c***category  f2a
c***type      single precision (sos-s, dsos-d)
c***keywords  brown's method, newton's method, nonlinear equations,
c             roots, solutions
c***author  watts, h. a., (snla)
c***description
c
c     sos solves a system of neq simultaneous nonlinear equations in
c     neq unknowns.  that is, it solves the problem   f(x)=0
c     where x is a vector with components  x(1),...,x(neq)  and  f
c     is a vector of nonlinear functions.  each equation is of the form
c
c               f (x(1),...,x(neq))=0     for k=1,...,neq.
c                k
c
c     the algorithm is based on an iterative method which is a
c     variation of newton's method using gaussian elimination
c     in a manner similar to the gauss-seidel process.  convergence
c     is roughly quadratic.  all partial derivatives required by
c     the algorithm are approximated by first difference quotients.
c     the convergence behavior of this code is affected by the
c     ordering of the equations, and it is advantageous to place linear
c     and mildly nonlinear equations first in the ordering.
c
c     actually, sos is merely an interfacing routine for
c     calling subroutine soseqs which embodies the solution
c     algorithm.  the purpose of this is to add greater
c     flexibility and ease of use for the prospective user.
c
c     soseqs calls the accompanying routine sossol, which solves special
c     triangular linear systems by back-substitution.
c
c     the user must supply a function subprogram which evaluates the
c     k-th equation only (k specified by soseqs) for each call
c     to the subprogram.
c
c     sos represents an implementation of the mathematical algorithm
c     described in the references below.  it is a modification of the
c     code sosnle written by h. a. watts in 1973.
c
c **********************************************************************
c   -input-
c
c     fnc -name of the function program which evaluates the equations.
c          this name must be in an external statement in the calling
c          program.  the user must supply fnc in the form fnc(x,k),
c          where x is the solution vector (which must be dimensioned
c          in fnc) and fnc returns the value of the k-th function.
c
c     neq -number of equations to be solved.
c
c     x   -solution vector.  initial guesses must be supplied.
c
c     rtolx -relative error tolerance used in the convergence criteria.
c          each solution component x(i) is checked by an accuracy test
c          of the form   abs(x(i)-xold(i)) .le. rtolx*abs(x(i))+atolx,
c          where xold(i) represents the previous iteration value.
c          rtolx must be non-negative.
c
c     atolx -absolute error tolerance used in the convergence criteria.
c          atolx must be non-negative.  if the user suspects some
c          solution component may be zero, he should set atolx to an
c          appropriate (depends on the scale of the remaining variables)
c          positive value for better efficiency.
c
c     tolf -residual error tolerance used in the convergence criteria.
c          convergence will be indicated if all residuals (values of the
c          functions or equations) are not bigger than tolf in
c          magnitude.  note that extreme care must be given in assigning
c          an appropriate value for tolf because this convergence test
c          is dependent on the scaling of the equations.  an
c          inappropriate value can cause premature termination of the
c          iteration process.
c
c     iflag -optional input indicator.  you must set  iflag=-1  if you
c          want to use any of the optional input items listed below.
c          otherwise set it to zero.
c
c     rw  -a real work array which is split apart by sos and used
c          internally by soseqs.
c
c     lrw -dimension of the rw array.  lrw must be at least
c                    1 + 6*neq + neq*(neq+1)/2
c
c     iw  -an integer work array which is split apart by sos and used
c          internally by soseqs.
c
c     liw -dimension of the iw array. liw must be at least  3 + neq.
c
c   -optional input-
c
c     iw(1) -internal printing parameter.  you must set  iw(1)=-1  if
c          you want the intermediate solution iterates to be printed.
c
c     iw(2) -iteration limit.  the maximum number of allowable
c          iterations can be specified, if desired.  to override the
c          default value of 50, set iw(2) to the number wanted.
c
c     remember, if you tell the code that you are using one of the
c               options (by setting iflag=-1), you must supply values
c               for both iw(1) and iw(2).
c
c **********************************************************************
c   -output-
c
c     x   -solution vector.
c
c     iflag -status indicator
c
c                         *** convergence to a solution ***
c
c          1 means satisfactory convergence to a solution was achieved.
c            each solution component x(i) satisfies the error tolerance
c            test   abs(x(i)-xold(i)) .le. rtolx*abs(x(i))+atolx.
c
c          2 means procedure converged to a solution such that all
c            residuals are at most tolf in magnitude,
c            abs(fnc(x,i)) .le. tolf.
c
c          3 means that conditions for both iflag=1 and iflag=2 hold.
c
c          4 means possible numerical convergence.  behavior indicates
c            limiting precision calculations as a result of user asking
c            for too much accuracy or else convergence is very slow.
c            residual norms and solution increment norms have
c            remained roughly constant over several consecutive
c            iterations.
c
c                         *** task interrupted ***
c
c          5 means the allowable number of iterations has been met
c            without obtaining a solution to the specified accuracy.
c            very slow convergence may be indicated.  examine the
c            approximate solution returned and see if the error
c            tolerances seem appropriate.
c
c          6 means the allowable number of iterations has been met and
c            the iterative process does not appear to be converging.
c            a local minimum may have been encountered or there may be
c            limiting precision difficulties.
c
c          7 means that the iterative scheme appears to be diverging.
c            residual norms and solution increment norms have
c            increased over several consecutive iterations.
c
c                         *** task cannot be continued ***
c
c          8 means that a jacobian-related matrix was singular.
c
c          9 means improper input parameters.
c
c          *** iflag should be examined after each call to   ***
c          *** sos with the appropriate action being taken.  ***
c
c
c     rw(1) -contains a norm of the residual.
c
c     iw(3) -contains the number of iterations used by the process.
c
c **********************************************************************
c***references  k. m. brown, solution of simultaneous nonlinear
c                 equations, algorithm 316, communications of the
c                 a.c.m. 10, (1967), pp. 728-729.
c               k. m. brown, a quadratically convergent newton-like
c                 method based upon gaussian elimination, siam journal
c                 on numerical analysis 6, (1969), pp. 560-569.
c***routines called  soseqs, xermsg
c***revision history  (yymmdd)
c   801001  date written
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900510  convert xerrwv calls to xermsg calls, changed prologue
c           comments to agree with dsos.  (rwc)
c   920501  reformatted the references section.  (wrb)
c***end prologue  sos
      dimension x(*), rw(*), iw(*)
      character*8 xern1
      character*16 xern3, xern4
      external fnc
c***first executable statement  sos
      inpflg = iflag
c
c     check for valid input
c
      if (neq .le. 0) then
         write (xern1, '(i8)') neq
         call xermsg ('slatec', 'sos', 'the number of equations ' //
     *      'must be a positive integer.  you have called the ' //
     *      'code with neq = ' // xern1, 1, 1)
         iflag = 9
      endif
c
      if (rtolx .lt. 0.0d0 .or. atolx .lt. 0.0d0) then
         write (xern3, '(1pe15.6)') atolx
         write (xern4, '(1pe15.6)') rtolx
         call xermsg ('slatec', 'sos', 'the error tolerances for ' //
     *      'the solution iterates cannot be negative. you have ' //
     *      'called the code with  rtolx = ' // xern3 //
     *      ' and atolx = ' // xern4,2, 1)
            iflag = 9
      endif
c
      if (tolf .lt. 0.0d0) then
         write (xern3, '(1pe15.6)') tolf
         call xermsg ('slatec', 'sos', 'the residual error ' //
     *      'tolerance must be non-negative.  you have called the ' //
     *      'code with tolf = ' // xern3, 3, 1)
            iflag = 9
      endif
c
      iprint = 0
      mxit = 50
      if (inpflg .eq. (-1)) then
         if (iw(1) .eq. (-1)) iprint = -1
         mxit = iw(2)
         if (mxit .le. 0) then
            write (xern1, '(i8)') mxit
            call xermsg ('slatec', 'sos', 'you have told the code ' //
     *         'to use optional in put items by setting  iflag=-1. ' //
     *         'however you have called the code with the maximum ' //
     *         'allowable number of iterations set to  iw(2) = ' //
     *         xern1, 4, 1)
            iflag = 9
         endif
      endif
c
      nc = (neq*(neq+1))/2
      if (lrw .lt. 1 + 6*neq + nc) then
         write (xern1, '(i8)') lrw
         call xermsg ('slatec', 'sos', 'dimension of the rw array ' //
     *      'must be at least 1 + 6*neq + neq*(neq+1)/2 .  you have ' //
     *      'called the code with lrw = ' // xern1, 5, 1)
         iflag = 9
      endif
c
      if (liw .lt. 3 + neq) then
         write (xern1, '(i8)') liw
         call xermsg ('slatec', 'sos', 'dimension of the iw array ' //
     *      'must be at least   3 + neq.  you have called the code ' //
     *      'with  liw = ' // xern1, 6, 1)
         iflag = 9
      endif
c
      if (iflag .ne. 9) then
         ncjs = 6
         nsrrc = 4
         nsri = 5
c
         k1 = nc + 2
         k2 = k1 + neq
         k3 = k2 + neq
         k4 = k3 + neq
         k5 = k4 + neq
         k6 = k5 + neq
c
         call soseqs(fnc, neq, x, rtolx, atolx, tolf, iflag, mxit, ncjs,
     1               nsrrc, nsri, iprint, rw(1), rw(2), nc, rw(k1),
     2               rw(k2), rw(k3), rw(k4), rw(k5), rw(k6), iw(4))
c
         iw(3) = mxit
      endif
      return
      end
