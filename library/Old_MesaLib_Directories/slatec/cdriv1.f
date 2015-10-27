*deck cdriv1
      subroutine cdriv1 (n, t, y, f, tout, mstate, eps, work, lenw,
     8   ierflg)
c***begin prologue  cdriv1
c***purpose  the function of cdriv1 is to solve n (200 or fewer)
c            ordinary differential equations of the form
c            dy(i)/dt = f(y(i),t), given the initial conditions
c            y(i) = yi.  cdriv1 allows complex-valued differential
c            equations.
c***library   slatec (sdrive)
c***category  i1a2, i1a1b
c***type      complex (sdriv1-s, ddriv1-d, cdriv1-c)
c***keywords  complex valued, gear's method, initial value problems,
c             ode, ordinary differential equations, sdrive, stiff
c***author  kahaner, d. k., (nist)
c             national institute of standards and technology
c             gaithersburg, md  20899
c           sutherland, c. d., (lanl)
c             mail stop d466
c             los alamos national laboratory
c             los alamos, nm  87545
c***description
c
c   version 92.1
c
c  i.  choosing the correct routine  ...................................
c
c     sdriv
c     ddriv
c     cdriv
c           these are the generic names for three packages for solving
c           initial value problems for ordinary differential equations.
c           sdriv uses single precision arithmetic.  ddriv uses double
c           precision arithmetic.  cdriv allows complex-valued
c           differential equations, integrated with respect to a single,
c           real, independent variable.
c
c    as an aid in selecting the proper program, the following is a
c    discussion of the important options or restrictions associated with
c    each program:
c
c      a. cdriv1 should be tried first for those routine problems with
c         no more than 200 differential equations (cdriv2 and cdriv3
c         have no such restriction.)  internally this routine has two
c         important technical defaults:
c           1. numerical approximation of the jacobian matrix of the
c              right hand side is used.
c           2. the stiff solver option is used.
c         most users of cdriv1 should not have to concern themselves
c         with these details.
c
c      b. cdriv2 should be considered for those problems for which
c         cdriv1 is inadequate.  for example, cdriv1 may have difficulty
c         with problems having zero initial conditions and zero
c         derivatives.  in this case cdriv2, with an appropriate value
c         of the parameter ewt, should perform more efficiently.  cdriv2
c         provides three important additional options:
c           1. the nonstiff equation solver (as well as the stiff
c              solver) is available.
c           2. the root-finding option is available.
c           3. the program can dynamically select either the non-stiff
c              or the stiff methods.
c         internally this routine also defaults to the numerical
c         approximation of the jacobian matrix of the right hand side.
c
c      c. cdriv3 is the most flexible, and hence the most complex, of
c         the programs.  its important additional features include:
c           1. the ability to exploit band structure in the jacobian
c              matrix.
c           2. the ability to solve some implicit differential
c              equations, i.e., those having the form:
c                   a(y,t)*dy/dt = f(y,t).
c           3. the option of integrating in the one step mode.
c           4. the option of allowing the user to provide a routine
c              which computes the analytic jacobian matrix of the right
c              hand side.
c           5. the option of allowing the user to provide a routine
c              which does all the matrix algebra associated with
c              corrections to the solution components.
c
c  ii.  parameters  ....................................................
c
c    the user should use parameter names in the call sequence of cdriv1
c    for those quantities whose value may be altered by cdriv1.  the
c    parameters in the call sequence are:
c
c    n      = (input) the number of differential equations, n .le. 200
c
c    t      = (real) the independent variable.  on input for the first
c             call, t is the initial point.  on output, t is the point
c             at which the solution is given.
c
c    y      = (complex) the vector of dependent variables.  y is used as
c             input on the first call, to set the initial values.  on
c             output, y is the computed solution vector.  this array y
c             is passed in the call sequence of the user-provided
c             routine f.  thus parameters required by f can be stored in
c             this array in components n+1 and above.  (note: changes by
c             the user to the first n components of this array will take
c             effect only after a restart, i.e., after setting mstate to
c             +1(-1).)
c
c    f      = a subroutine supplied by the user.  the name must be
c             declared external in the user's calling program.  this
c             subroutine is of the form:
c                   subroutine f (n, t, y, ydot)
c                   complex y(*), ydot(*)
c                     .
c                     .
c                   ydot(1) = ...
c                     .
c                     .
c                   ydot(n) = ...
c                   end (sample)
c             this computes ydot = f(y,t), the right hand side of the
c             differential equations.  here y is a vector of length at
c             least n.  the actual length of y is determined by the
c             user's declaration in the program which calls cdriv1.
c             thus the dimensioning of y in f, while required by fortran
c             convention, does not actually allocate any storage.  when
c             this subroutine is called, the first n components of y are
c             intermediate approximations to the solution components.
c             the user should not alter these values.  here ydot is a
c             vector of length n.  the user should only compute ydot(i)
c             for i from 1 to n.  normally a return from f passes
c             control back to  cdriv1.  however, if the user would like
c             to abort the calculation, i.e., return control to the
c             program which calls cdriv1, he should set n to zero.
c             cdriv1 will signal this by returning a value of mstate
c             equal to +5(-5).  altering the value of n in f has no
c             effect on the value of n in the call sequence of cdriv1.
c
c    tout   = (input, real) the point at which the solution is desired.
c
c    mstate = an integer describing the status of integration.  the user
c             must initialize mstate to +1 or -1.  if mstate is
c             positive, the routine will integrate past tout and
c             interpolate the solution.  this is the most efficient
c             mode.  if mstate is negative, the routine will adjust its
c             internal step to reach tout exactly (useful if a
c             singularity exists beyond tout.)  the meaning of the
c             magnitude of mstate:
c               1  (input) means the first call to the routine.  this
c                  value must be set by the user.  on all subsequent
c                  calls the value of mstate should be tested by the
c                  user.  unless cdriv1 is to be reinitialized, only the
c                  sign of mstate may be changed by the user.  (as a
c                  convenience to the user who may wish to put out the
c                  initial conditions, cdriv1 can be called with
c                  mstate=+1(-1), and tout=t.  in this case the program
c                  will return with mstate unchanged, i.e.,
c                  mstate=+1(-1).)
c               2  (output) means a successful integration.  if a normal
c                  continuation is desired (i.e., a further integration
c                  in the same direction), simply advance tout and call
c                  again.  all other parameters are automatically set.
c               3  (output)(unsuccessful) means the integrator has taken
c                  1000 steps without reaching tout.  the user can
c                  continue the integration by simply calling cdriv1
c                  again.
c               4  (output)(unsuccessful) means too much accuracy has
c                  been requested.  eps has been increased to a value
c                  the program estimates is appropriate.  the user can
c                  continue the integration by simply calling cdriv1
c                  again.
c               5  (output)(unsuccessful) n has been set to zero in
c                  subroutine f.
c               6  (output)(successful) for mstate negative, t is beyond
c                  tout.  the solution was obtained by interpolation.
c                  the user can continue the integration by simply
c                  advancing tout and calling cdriv1 again.
c               7  (output)(unsuccessful) the solution could not be
c                  obtained.  the value of ierflg (see description
c                  below) for a "recoverable" situation indicates the
c                  type of difficulty encountered: either an illegal
c                  value for a parameter or an inability to continue the
c                  solution.  for this condition the user should take
c                  corrective action and reset mstate to +1(-1) before
c                  calling cdriv1 again.  otherwise the program will
c                  terminate the run.
c
c    eps    = (real) on input, the requested relative accuracy in all
c             solution components.  on output, the adjusted relative
c             accuracy if the input value was too small.  the value of
c             eps should be set as large as is reasonable, because the
c             amount of work done by cdriv1 increases as eps decreases.
c
c    work
c    lenw   = (input)
c             work is an array of lenw complex words used
c             internally for temporary storage.  the user must allocate
c             space for this array in the calling program by a statement
c             such as
c                       complex work(...)
c             the length of work should be at least n*n + 11*n + 300
c             and lenw should be set to the value used.  the contents of
c             work should not be disturbed between calls to cdriv1.
c
c    ierflg = an error flag.  the error number associated with a
c             diagnostic message (see section iv-a below) is the same as
c             the corresponding value of ierflg.  the meaning of ierflg:
c               0  the routine completed successfully. (no message is
c                  issued.)
c               3  (warning) the number of steps required to reach tout
c                  exceeds 1000 .
c               4  (warning) the value of eps is too small.
c              11  (warning) for mstate negative, t is beyond tout.
c                  the solution was obtained by interpolation.
c              15  (warning) the integration step size is below the
c                  roundoff level of t.  (the program issues this
c                  message as a warning but does not return control to
c                  the user.)
c              21  (recoverable) n is greater than 200 .
c              22  (recoverable) n is not positive.
c              26  (recoverable) the magnitude of mstate is either 0 or
c                  greater than 7 .
c              27  (recoverable) eps is less than zero.
c              32  (recoverable) insufficient storage has been allocated
c                  for the work array.
c              41  (recoverable) the integration step size has gone
c                  to zero.
c              42  (recoverable) the integration step size has been
c                  reduced about 50 times without advancing the
c                  solution.  the problem setup may not be correct.
c             999  (fatal) the magnitude of mstate is 7 .
c
c  iii.  usage  ........................................................
c
c                program sample
c                external f
c                complex alfa
c                real eps, t, tout
c          c                                n is the number of equations
c                parameter(alfa = (1.e0, 1.e0), n = 3,
c               8          lenw = n*n + 11*n + 300)
c                complex work(lenw), y(n+1)
c          c                                               initial point
c                t = 0.00001e0
c          c                                      set initial conditions
c                y(1) = 10.e0
c                y(2) = 0.e0
c                y(3) = 10.e0
c          c                                              pass parameter
c                y(4) = alfa
c                tout = t
c                mstate = 1
c                eps = .001e0
c           10   call cdriv1 (n, t, y, f, tout, mstate, eps, work, lenw,
c               8             ierflg)
c                if (mstate .gt. 2) stop
c                write(*, '(5e12.3)') tout, (y(i), i=1,3)
c                tout = 10.e0*tout
c                if (tout .lt. 50.e0) go to 10
c                end
c
c                subroutine f (n, t, y, ydot)
c                complex alfa, y(*), ydot(*)
c                real t
c                alfa = y(n+1)
c                ydot(1) = 1.e0 + alfa*(y(2) - y(1)) - y(1)*y(3)
c                ydot(2) = alfa*(y(1) - y(2)) - y(2)*y(3)
c                ydot(3) = 1.e0 - y(3)*(y(1) + y(2))
c                end
c
c  iv.  other communication to the user  ...............................
c
c    a. the solver communicates to the user through the parameters
c       above.  in addition it writes diagnostic messages through the
c       standard error handling program xermsg.  a complete description
c       of xermsg is given in "guide to the slatec common mathematical
c       library" by kirby w. fong et al..  at installations which do not
c       have this error handling package the short but serviceable
c       routine, xermsg, available with this package, can be used.  that
c       program uses the file named output to transmit messages.
c
c    b. the number of evaluations of the right hand side can be found
c       in the work array in the location determined by:
c            lenw - (n + 50) + 4
c
c  v.  remarks  ........................................................
c
c    for other information, see section iv of the writeup for cdriv3.
c
c***references  c. w. gear, numerical initial value problems in
c                 ordinary differential equations, prentice-hall, 1971.
c***routines called  cdriv3, xermsg
c***revision history  (yymmdd)
c   790601  date written
c   900329  initial submission to slatec.
c***end prologue  cdriv1
      external f
      complex work(*), y(*)
      real eps, ewtcom(1), hmax, t, tout
      integer i, idliw, ierflg, ierror, impl, leniw, lenw, lenwcm,
     8        lnwchk, mint, miter, ml, mstate, mu, mxn, mxord, mxstep,
     8        n, nde, nroot, nstate, ntask
      parameter(mxn = 200, idliw = 50)
      integer iwork(idliw+mxn)
      character intgr1*8
      parameter(nroot = 0, ierror = 2, mint = 2, miter = 2, impl = 0,
     8          mxord = 5, mxstep = 1000)
      data ewtcom(1) /1.e0/
c***first executable statement  cdriv1
      if (abs(mstate) .eq. 0 .or. abs(mstate) .gt. 7) then
        write(intgr1, '(i8)') mstate
        ierflg = 26
        call xermsg('slatec', 'cdriv1',
     8  'illegal input.  the magnitude of mstate, '//intgr1//
     8  ', is not in the range 1 to 6 .', ierflg, 1)
        mstate = sign(7, mstate)
        return
      else if (abs(mstate) .eq. 7) then
        ierflg = 999
        call xermsg('slatec', 'cdriv1',
     8  'illegal input.  the magnitude of mstate is 7 .', ierflg, 2)
        return
      end if
      if (n .gt. mxn) then
        write(intgr1, '(i8)') n
        ierflg = 21
        call xermsg('slatec', 'cdriv1',
     8  'illegal input.  the number of equations, '//intgr1//
     8  ', is greater than the maximum allowed: 200 .', ierflg, 1)
        mstate = sign(7, mstate)
        return
      end if
      if (mstate .gt. 0) then
        nstate = mstate
        ntask = 1
      else
        nstate = - mstate
        ntask = 3
      end if
      hmax = 2.e0*abs(tout - t)
      leniw = n + idliw
      lenwcm = lenw - leniw
      if (lenwcm .lt. (n*n + 10*n + 250)) then
        lnwchk = n*n + 10*n + 250 + leniw
        write(intgr1, '(i8)') lnwchk
        ierflg = 32
        call xermsg('slatec', 'cdriv1',
     8  'insufficient storage allocated for the work array.  '//
     8  'the required storage is at least '//intgr1//' .', ierflg, 1)
        mstate = sign(7, mstate)
        return
      end if
      if (nstate .ne. 1) then
        do 20 i = 1,leniw
 20       iwork(i) = work(i+lenwcm)
      end if
      call cdriv3 (n, t, y, f, nstate, tout, ntask, nroot, eps, ewtcom,
     8             ierror, mint, miter, impl, ml, mu, mxord, hmax, work,
     8             lenwcm, iwork, leniw, f, f, nde, mxstep, f, f,
     8             ierflg)
      do 40 i = 1,leniw
 40     work(i+lenwcm) = iwork(i)
      if (nstate .le. 4) then
        mstate = sign(nstate, mstate)
      else if (nstate .eq. 6) then
        mstate = sign(5, mstate)
      else if (ierflg .eq. 11) then
        mstate = sign(6, mstate)
      else if (ierflg .gt. 11) then
        mstate = sign(7, mstate)
      end if
      return
      end
