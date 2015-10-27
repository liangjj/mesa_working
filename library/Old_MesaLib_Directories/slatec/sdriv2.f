*deck sdriv2
      subroutine sdriv2 (n, t, y, f, tout, mstate, nroot, eps, ewt,
     8   mint, work, lenw, iwork, leniw, g, ierflg)
c***begin prologue  sdriv2
c***purpose  the function of sdriv2 is to solve n ordinary differential
c            equations of the form dy(i)/dt = f(y(i),t), given the
c            initial conditions y(i) = yi.  the program has options to
c            allow the solution of both stiff and non-stiff differential
c            equations.  sdriv2 uses single precision arithmetic.
c***library   slatec (sdrive)
c***category  i1a2, i1a1b
c***type      single precision (sdriv2-s, ddriv2-d, cdriv2-c)
c***keywords  gear's method, initial value problems, ode,
c             ordinary differential equations, sdrive, single precision,
c             stiff
c***author  kahaner, d. k., (nist)
c             national institute of standards and technology
c             gaithersburg, md  20899
c           sutherland, c. d., (lanl)
c             mail stop d466
c             los alamos national laboratory
c             los alamos, nm  87545
c***description
c
c  i.  parameters  .....................................................
c
c    the user should use parameter names in the call sequence of sdriv2
c    for those quantities whose value may be altered by sdriv2.  the
c    parameters in the call sequence are:
c
c    n      = (input) the number of differential equations.
c
c    t      = the independent variable.  on input for the first call, t
c             is the initial point.  on output, t is the point at which
c             the solution is given.
c
c    y      = the vector of dependent variables.  y is used as input on
c             the first call, to set the initial values.  on output, y
c             is the computed solution vector.  this array y is passed
c             in the call sequence of the user-provided routines f and
c             g.  thus parameters required by f and g can be stored in
c             this array in components n+1 and above.  (note: changes
c             by the user to the first n components of this array will
c             take effect only after a restart, i.e., after setting
c             mstate to +1(-1).)
c
c    f      = a subroutine supplied by the user.  the name must be
c             declared external in the user's calling program.  this
c             subroutine is of the form:
c                   subroutine f (n, t, y, ydot)
c                   real y(*), ydot(*)
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
c             user's declaration in the program which calls sdriv2.
c             thus the dimensioning of y in f, while required by fortran
c             convention, does not actually allocate any storage.  when
c             this subroutine is called, the first n components of y are
c             intermediate approximations to the solution components.
c             the user should not alter these values.  here ydot is a
c             vector of length n.  the user should only compute ydot(i)
c             for i from 1 to n.  normally a return from f passes
c             control back to  sdriv2.  however, if the user would like
c             to abort the calculation, i.e., return control to the
c             program which calls sdriv2, he should set n to zero.
c             sdriv2 will signal this by returning a value of mstate
c             equal to +6(-6).  altering the value of n in f has no
c             effect on the value of n in the call sequence of sdriv2.
c
c    tout   = (input) the point at which the solution is desired.
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
c                  user.  unless sdriv2 is to be reinitialized, only the
c                  sign of mstate may be changed by the user.  (as a
c                  convenience to the user who may wish to put out the
c                  initial conditions, sdriv2 can be called with
c                  mstate=+1(-1), and tout=t.  in this case the program
c                  will return with mstate unchanged, i.e.,
c                  mstate=+1(-1).)
c               2  (output) means a successful integration.  if a normal
c                  continuation is desired (i.e., a further integration
c                  in the same direction), simply advance tout and call
c                  again.  all other parameters are automatically set.
c               3  (output)(unsuccessful) means the integrator has taken
c                  1000 steps without reaching tout.  the user can
c                  continue the integration by simply calling sdriv2
c                  again.  other than an error in problem setup, the
c                  most likely cause for this condition is trying to
c                  integrate a stiff set of equations with the non-stiff
c                  integrator option. (see description of mint below.)
c               4  (output)(unsuccessful) means too much accuracy has
c                  been requested.  eps has been increased to a value
c                  the program estimates is appropriate.  the user can
c                  continue the integration by simply calling sdriv2
c                  again.
c               5  (output) a root was found at a point less than tout.
c                  the user can continue the integration toward tout by
c                  simply calling sdriv2 again.
c               6  (output)(unsuccessful) n has been set to zero in
c                  subroutine f.
c               7  (output)(unsuccessful) n has been set to zero in
c                  function g.  see description of g below.
c               8  (output)(successful) for mstate negative, t is beyond
c                  tout.  the solution was obtained by interpolation.
c                  the user can continue the integration by simply
c                  advancing tout and calling sdriv2 again.
c               9  (output)(unsuccessful) the solution could not be
c                  obtained.  the value of ierflg (see description
c                  below) for a "recoverable" situation indicates the
c                  type of difficulty encountered: either an illegal
c                  value for a parameter or an inability to continue the
c                  solution.  for this condition the user should take
c                  corrective action and reset mstate to +1(-1) before
c                  calling sdriv2 again.  otherwise the program will
c                  terminate the run.
c
c    nroot  = (input) the number of equations whose roots are desired.
c             if nroot is zero, the root search is not active.  this
c             option is useful for obtaining output at points which are
c             not known in advance, but depend upon the solution, e.g.,
c             when some solution component takes on a specified value.
c             the root search is carried out using the user-written
c             function g (see description of g below.)  sdriv2 attempts
c             to find the value of t at which one of the equations
c             changes sign.  sdriv2 can find at most one root per
c             equation per internal integration step, and will then
c             return the solution either at tout or at a root, whichever
c             occurs first in the direction of integration.  the initial
c             point is never reported as a root.  the index of the
c             equation whose root is being reported is stored in the
c             sixth element of iwork.
c             note: nroot is never altered by this program.
c
c    eps    = on input, the requested relative accuracy in all solution
c             components.  eps = 0 is allowed.  on output, the adjusted
c             relative accuracy if the input value was too small.  the
c             value of eps should be set as large as is reasonable,
c             because the amount of work done by sdriv2 increases as
c             eps decreases.
c
c    ewt    = (input) problem zero, i.e., the smallest physically
c             meaningful value for the solution.  this is used inter-
c             nally to compute an array ywt(i) = max(abs(y(i)), ewt).
c             one step error estimates divided by ywt(i) are kept less
c             than eps.  setting ewt to zero provides pure relative
c             error control.  however, setting ewt smaller than
c             necessary can adversely affect the running time.
c
c    mint   = (input) the integration method flag.
c               mint = 1  means the adams methods, and is used for
c                         non-stiff problems.
c               mint = 2  means the stiff methods of gear (i.e., the
c                         backward differentiation formulas), and is
c                         used for stiff problems.
c               mint = 3  means the program dynamically selects the
c                         adams methods when the problem is non-stiff
c                         and the gear methods when the problem is
c                         stiff.
c             mint may not be changed without restarting, i.e., setting
c             the magnitude of mstate to 1.
c
c    work
c    lenw   = (input)
c             work is an array of lenw real words used
c             internally for temporary storage.  the user must allocate
c             space for this array in the calling program by a statement
c             such as
c                       real work(...)
c             the length of work should be at least
c               16*n + 2*nroot + 250         if mint is 1, or
c               n*n + 10*n + 2*nroot + 250   if mint is 2, or
c               n*n + 17*n + 2*nroot + 250   if mint is 3,
c             and lenw should be set to the value used.  the contents of
c             work should not be disturbed between calls to sdriv2.
c
c    iwork
c    leniw  = (input)
c             iwork is an integer array of length leniw used internally
c             for temporary storage.  the user must allocate space for
c             this array in the calling program by a statement such as
c                       integer iwork(...)
c             the length of iwork should be at least
c               50      if mint is 1, or
c               n+50    if mint is 2 or 3,
c             and leniw should be set to the value used.  the contents
c             of iwork should not be disturbed between calls to sdriv2.
c
c    g      = a real fortran function supplied by the user
c             if nroot is not 0.  in this case, the name must be
c             declared external in the user's calling program.  g is
c             repeatedly called with different values of iroot to
c             obtain the value of each of the nroot equations for which
c             a root is desired.  g is of the form:
c                   real function g (n, t, y, iroot)
c                   real y(*)
c                   go to (10, ...), iroot
c              10   g = ...
c                     .
c                     .
c                   end (sample)
c             here, y is a vector of length at least n, whose first n
c             components are the solution components at the point t.
c             the user should not alter these values.  the actual length
c             of y is determined by the user's declaration in the
c             program which calls sdriv2.  thus the dimensioning of y in
c             g, while required by fortran convention, does not actually
c             allocate any storage.  normally a return from g passes
c             control back to  sdriv2.  however, if the user would like
c             to abort the calculation, i.e., return control to the
c             program which calls sdriv2, he should set n to zero.
c             sdriv2 will signal this by returning a value of mstate
c             equal to +7(-7).  in this case, the index of the equation
c             being evaluated is stored in the sixth element of iwork.
c             altering the value of n in g has no effect on the value of
c             n in the call sequence of sdriv2.
c
c    ierflg = an error flag.  the error number associated with a
c             diagnostic message (see section ii-a below) is the same as
c             the corresponding value of ierflg.  the meaning of ierflg:
c               0  the routine completed successfully. (no message is
c                  issued.)
c               3  (warning) the number of steps required to reach tout
c                  exceeds mxstep.
c               4  (warning) the value of eps is too small.
c              11  (warning) for mstate negative, t is beyond tout.
c                  the solution was obtained by interpolation.
c              15  (warning) the integration step size is below the
c                  roundoff level of t.  (the program issues this
c                  message as a warning but does not return control to
c                  the user.)
c              22  (recoverable) n is not positive.
c              23  (recoverable) mint is less than 1 or greater than 3 .
c              26  (recoverable) the magnitude of mstate is either 0 or
c                  greater than 9 .
c              27  (recoverable) eps is less than zero.
c              32  (recoverable) insufficient storage has been allocated
c                  for the work array.
c              33  (recoverable) insufficient storage has been allocated
c                  for the iwork array.
c              41  (recoverable) the integration step size has gone
c                  to zero.
c              42  (recoverable) the integration step size has been
c                  reduced about 50 times without advancing the
c                  solution.  the problem setup may not be correct.
c             999  (fatal) the magnitude of mstate is 9 .
c
c  ii.  other communication to the user  ...............................
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
c    b. the first three elements of work and the first five elements of
c       iwork will contain the following statistical data:
c         avgh     the average step size used.
c         hused    the step size last used (successfully).
c         avgord   the average order used.
c         imxerr   the index of the element of the solution vector that
c                  contributed most to the last error test.
c         nqused   the order last used (successfully).
c         nstep    the number of steps taken since last initialization.
c         nfe      the number of evaluations of the right hand side.
c         nje      the number of evaluations of the jacobian matrix.
c
c  iii.  remarks  ......................................................
c
c    a. on any return from sdriv2 all information necessary to continue
c       the calculation is contained in the call sequence parameters,
c       including the work arrays.  thus it is possible to suspend one
c       problem, integrate another, and then return to the first.
c
c    b. if this package is to be used in an overlay situation, the user
c       must declare in the primary overlay the variables in the call
c       sequence to sdriv2.
c
c    c. when the routine g is not required, difficulties associated with
c       an unsatisfied external can be avoided by using the name of the
c       routine which calculates the right hand side of the differential
c       equations in place of g in the call sequence of sdriv2.
c
c  iv.  usage  .........................................................
c
c               program sample
c               external f
c               parameter(mint = 1, nroot = 0, n = ...,
c              8          lenw = 16*n + 2*nroot + 250, leniw = 50)
c         c                                 n is the number of equations
c               real eps, ewt, t, tout, work(lenw), y(n)
c               integer iwork(leniw)
c               open(file='tape6', unit=6, status='new')
c         c                                                initial point
c               t = 0.
c         c                                       set initial conditions
c               do 10 i = 1,n
c          10     y(i) = ...
c               tout = t
c               ewt = ...
c               mstate = 1
c               eps = ...
c          20   call sdriv2 (n, t, y, f, tout, mstate, nroot, eps, ewt,
c              8             mint, work, lenw, iwork, leniw, f, ierflg)
c         c                                 next to last argument is not
c         c                                    f if rootfinding is used.
c               if (mstate .gt. 2) stop
c               write(6, 100) tout, (y(i), i=1,n)
c               tout = tout + 1.
c               if (tout .le. 10.) go to 20
c          100  format(...)
c               end (sample)
c
c***references  c. w. gear, numerical initial value problems in
c                 ordinary differential equations, prentice-hall, 1971.
c***routines called  sdriv3, xermsg
c***revision history  (yymmdd)
c   790601  date written
c   900329  initial submission to slatec.
c***end prologue  sdriv2
      external f, g
      real eps, ewt, ewtcom(1), g, hmax, t, tout,
     8     work(*), y(*)
      integer iwork(*)
      integer ierflg, ierror, impl, leniw, lenw, mint, miter, ml,
     8        mstate, mu, mxord, mxstep, n, nde, nroot, nstate, ntask
      character intgr1*8
      parameter(impl = 0, mxstep = 1000)
c***first executable statement  sdriv2
      if (abs(mstate) .eq. 9) then
        ierflg = 999
        call xermsg('slatec', 'sdriv2',
     8  'illegal input.  the magnitude of mstate is 9 .',
     8  ierflg, 2)
        return
      else if (abs(mstate) .eq. 0 .or. abs(mstate) .gt. 9) then
        write(intgr1, '(i8)') mstate
        ierflg = 26
        call xermsg('slatec', 'sdriv2',
     8  'illegal input.  the magnitude of mstate, '//intgr1//
     8  ' is not in the range 1 to 8 .', ierflg, 1)
        mstate = sign(9, mstate)
        return
      end if
      if (mint .lt. 1 .or. mint .gt. 3) then
        write(intgr1, '(i8)') mint
        ierflg = 23
        call xermsg('slatec', 'sdriv2',
     8  'illegal input.  improper value for the integration method '//
     8  'flag, '//intgr1//' .', ierflg, 1)
        mstate = sign(9, mstate)
        return
      end if
      if (mstate .ge. 0) then
        nstate = mstate
        ntask = 1
      else
        nstate = - mstate
        ntask = 3
      end if
      ewtcom(1) = ewt
      if (ewt .ne. 0.e0) then
        ierror = 3
      else
        ierror = 2
      end if
      if (mint .eq. 1) then
        miter = 0
        mxord = 12
      else if (mint .eq. 2) then
        miter = 2
        mxord = 5
      else if (mint .eq. 3) then
        miter = 2
        mxord = 12
      end if
      hmax = 2.e0*abs(tout - t)
      call sdriv3 (n, t, y, f, nstate, tout, ntask, nroot, eps, ewtcom,
     8             ierror, mint, miter, impl, ml, mu, mxord, hmax, work,
     8             lenw, iwork, leniw, f, f, nde, mxstep, g, f, ierflg)
      if (nstate .le. 7) then
        mstate = sign(nstate, mstate)
      else if (nstate .eq. 11) then
        mstate = sign(8, mstate)
      else if (nstate .gt. 11) then
        mstate = sign(9, mstate)
      end if
      return
      end
