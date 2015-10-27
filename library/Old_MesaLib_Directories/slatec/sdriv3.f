*deck sdriv3
      subroutine sdriv3 (n, t, y, f, nstate, tout, ntask, nroot, eps,
     8   ewt, ierror, mint, miter, impl, ml, mu, mxord, hmax, work,
     8   lenw, iwork, leniw, jacobn, fa, nde, mxstep, g, users, ierflg)
c***begin prologue  sdriv3
c***purpose  the function of sdriv3 is to solve n ordinary differential
c            equations of the form dy(i)/dt = f(y(i),t), given the
c            initial conditions y(i) = yi.  the program has options to
c            allow the solution of both stiff and non-stiff differential
c            equations.  other important options are available.  sdriv3
c            uses single precision arithmetic.
c***library   slatec (sdrive)
c***category  i1a2, i1a1b
c***type      single precision (sdriv3-s, ddriv3-d, cdriv3-c)
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
c  i.  abstract  .......................................................
c
c    the primary function of sdriv3 is to solve n ordinary differential
c    equations of the form dy(i)/dt = f(y(i),t), given the initial
c    conditions y(i) = yi.  the program has options to allow the
c    solution of both stiff and non-stiff differential equations.  in
c    addition, sdriv3 may be used to solve:
c      1. the initial value problem, a*dy(i)/dt = f(y(i),t), where a is
c         a non-singular matrix depending on y and t.
c      2. the hybrid differential/algebraic initial value problem,
c         a*dy(i)/dt = f(y(i),t), where a is a vector (whose values may
c         depend upon y and t) some of whose components will be zero
c         corresponding to those equations which are algebraic rather
c         than differential.
c    sdriv3 is to be called once for each output point of t.
c
c  ii.  parameters  ....................................................
c
c    the user should use parameter names in the call sequence of sdriv3
c    for those quantities whose value may be altered by sdriv3.  the
c    parameters in the call sequence are:
c
c    n      = (input) the number of dependent functions whose solution
c             is desired.  n must not be altered during a problem.
c
c    t      = the independent variable.  on input for the first call, t
c             is the initial point.  on output, t is the point at which
c             the solution is given.
c
c    y      = the vector of dependent variables.  y is used as input on
c             the first call, to set the initial values.  on output, y
c             is the computed solution vector.  this array y is passed
c             in the call sequence of the user-provided routines f,
c             jacobn, fa, users, and g.  thus parameters required by
c             those routines can be stored in this array in components
c             n+1 and above.  (note: changes by the user to the first
c             n components of this array will take effect only after a
c             restart, i.e., after setting nstate to 1 .)
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
c             user's declaration in the program which calls sdriv3.
c             thus the dimensioning of y in f, while required by fortran
c             convention, does not actually allocate any storage.  when
c             this subroutine is called, the first n components of y are
c             intermediate approximations to the solution components.
c             the user should not alter these values.  here ydot is a
c             vector of length n.  the user should only compute ydot(i)
c             for i from 1 to n.  normally a return from f passes
c             control back to  sdriv3.  however, if the user would like
c             to abort the calculation, i.e., return control to the
c             program which calls sdriv3, he should set n to zero.
c             sdriv3 will signal this by returning a value of nstate
c             equal to 6 .  altering the value of n in f has no effect
c             on the value of n in the call sequence of sdriv3.
c
c    nstate = an integer describing the status of integration.  the
c             meaning of nstate is as follows:
c               1  (input) means the first call to the routine.  this
c                  value must be set by the user.  on all subsequent
c                  calls the value of nstate should be tested by the
c                  user, but must not be altered.  (as a convenience to
c                  the user who may wish to put out the initial
c                  conditions, sdriv3 can be called with nstate=1, and
c                  tout=t.  in this case the program will return with
c                  nstate unchanged, i.e., nstate=1.)
c               2  (output) means a successful integration.  if a normal
c                  continuation is desired (i.e., a further integration
c                  in the same direction), simply advance tout and call
c                  again.  all other parameters are automatically set.
c               3  (output)(unsuccessful) means the integrator has taken
c                  mxstep steps without reaching tout.  the user can
c                  continue the integration by simply calling sdriv3
c                  again.
c               4  (output)(unsuccessful) means too much accuracy has
c                  been requested.  eps has been increased to a value
c                  the program estimates is appropriate.  the user can
c                  continue the integration by simply calling sdriv3
c                  again.
c               5  (output) a root was found at a point less than tout.
c                  the user can continue the integration toward tout by
c                  simply calling sdriv3 again.
c               6  (output)(unsuccessful) n has been set to zero in
c                  subroutine f.
c               7  (output)(unsuccessful) n has been set to zero in
c                  function g.  see description of g below.
c               8  (output)(unsuccessful) n has been set to zero in
c                  subroutine jacobn.  see description of jacobn below.
c               9  (output)(unsuccessful) n has been set to zero in
c                  subroutine fa.  see description of fa below.
c              10  (output)(unsuccessful) n has been set to zero in
c                  subroutine users.  see description of users below.
c              11  (output)(successful) for ntask = 2 or 3, t is beyond
c                  tout.  the solution was obtained by interpolation.
c                  the user can continue the integration by simply
c                  advancing tout and calling sdriv3 again.
c              12  (output)(unsuccessful) the solution could not be
c                  obtained.  the value of ierflg (see description
c                  below) for a "recoverable" situation indicates the
c                  type of difficulty encountered: either an illegal
c                  value for a parameter or an inability to continue the
c                  solution.  for this condition the user should take
c                  corrective action and reset nstate to 1 before
c                  calling sdriv3 again.  otherwise the program will
c                  terminate the run.
c
c    tout   = (input) the point at which the solution is desired.  the
c             position of tout relative to t on the first call
c             determines the direction of integration.
c
c    ntask  = (input) an index specifying the manner of returning the
c             solution, according to the following:
c               ntask = 1  means sdriv3 will integrate past tout and
c                          interpolate the solution.  this is the most
c                          efficient mode.
c               ntask = 2  means sdriv3 will return the solution after
c                          each internal integration step, or at tout,
c                          whichever comes first.  in the latter case,
c                          the program integrates exactly to tout.
c               ntask = 3  means sdriv3 will adjust its internal step to
c                          reach tout exactly (useful if a singularity
c                          exists beyond tout.)
c
c    nroot  = (input) the number of equations whose roots are desired.
c             if nroot is zero, the root search is not active.  this
c             option is useful for obtaining output at points which are
c             not known in advance, but depend upon the solution, e.g.,
c             when some solution component takes on a specified value.
c             the root search is carried out using the user-written
c             function g (see description of g below.)  sdriv3 attempts
c             to find the value of t at which one of the equations
c             changes sign.  sdriv3 can find at most one root per
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
c             because the amount of work done by sdriv3 increases as eps
c             decreases.
c
c    ewt    = (input) problem zero, i.e., the smallest, nonzero,
c             physically meaningful value for the solution.  (array,
c             possibly of length one.  see following description of
c             ierror.)  setting ewt smaller than necessary can adversely
c             affect the running time.
c
c    ierror = (input) error control indicator.  a value of 3 is
c             suggested for most problems.  other choices and detailed
c             explanations of ewt and ierror are given below for those
c             who may need extra flexibility.
c
c             these last three input quantities eps, ewt and ierror
c             control the accuracy of the computed solution.  ewt and
c             ierror are used internally to compute an array ywt.  one
c             step error estimates divided by ywt(i) are kept less than
c             eps in root mean square norm.
c                 ierror (set by the user) =
c                 1  means ywt(i) = 1. (absolute error control)
c                                   ewt is ignored.
c                 2  means ywt(i) = abs(y(i)),  (relative error control)
c                                   ewt is ignored.
c                 3  means ywt(i) = max(abs(y(i)), ewt(1)).
c                 4  means ywt(i) = max(abs(y(i)), ewt(i)).
c                    this choice is useful when the solution components
c                    have differing scales.
c                 5  means ywt(i) = ewt(i).
c             if ierror is 3, ewt need only be dimensioned one.
c             if ierror is 4 or 5, the user must dimension ewt at least
c             n, and set its values.
c
c    mint   = (input) the integration method indicator.
c               mint = 1  means the adams methods, and is used for
c                         non-stiff problems.
c               mint = 2  means the stiff methods of gear (i.e., the
c                         backward differentiation formulas), and is
c                         used for stiff problems.
c               mint = 3  means the program dynamically selects the
c                         adams methods when the problem is non-stiff
c                         and the gear methods when the problem is
c                         stiff.  when using the adams methods, the
c                         program uses a value of miter=0; when using
c                         the gear methods, the program uses the value
c                         of miter provided by the user.  only a value
c                         of impl = 0 and a value of miter = 1, 2, 4, or
c                         5 is allowed for this option.  the user may
c                         not alter the value of mint or miter without
c                         restarting, i.e., setting nstate to 1.
c
c    miter  = (input) the iteration method indicator.
c               miter = 0  means functional iteration.  this value is
c                          suggested for non-stiff problems.
c               miter = 1  means chord method with analytic jacobian.
c                          in this case, the user supplies subroutine
c                          jacobn (see description below).
c               miter = 2  means chord method with jacobian calculated
c                          internally by finite differences.
c               miter = 3  means chord method with corrections computed
c                          by the user-written routine users (see
c                          description of users below.)  this option
c                          allows all matrix algebra and storage
c                          decisions to be made by the user.  when using
c                          a value of miter = 3, the subroutine fa is
c                          not required, even if impl is not 0.  for
c                          further information on using this option, see
c                          section iv-e below.
c               miter = 4  means the same as miter = 1 but the a and
c                          jacobian matrices are assumed to be banded.
c               miter = 5  means the same as miter = 2 but the a and
c                          jacobian matrices are assumed to be banded.
c
c    impl   = (input) the implicit method indicator.
c               impl = 0    means solving dy(i)/dt = f(y(i),t).
c               impl = 1    means solving a*dy(i)/dt = f(y(i),t), non-
c                           singular a (see description of fa below.)
c                           only mint = 1 or 2, and miter = 1, 2, 3, 4,
c                           or 5 are allowed for this option.
c               impl = 2,3  means solving certain systems of hybrid
c                           differential/algebraic equations (see
c                           description of fa below.)  only mint = 2 and
c                           miter = 1, 2, 3, 4, or 5, are allowed for
c                           this option.
c               the value of impl must not be changed during a problem.
c
c    ml     = (input) the lower half-bandwidth in the case of a banded
c             a or jacobian matrix.  (i.e., maximum(r-c) for nonzero
c             a(r,c).)
c
c    mu     = (input) the upper half-bandwidth in the case of a banded
c             a or jacobian matrix.  (i.e., maximum(c-r).)
c
c    mxord  = (input) the maximum order desired. this is .le. 12 for
c             the adams methods and .le. 5 for the gear methods.  normal
c             value is 12 and 5, respectively.  if mint is 3, the
c             maximum order used will be min(mxord, 12) when using the
c             adams methods, and min(mxord, 5) when using the gear
c             methods.  mxord must not be altered during a problem.
c
c    hmax   = (input) the maximum magnitude of the step size that will
c             be used for the problem.  this is useful for ensuring that
c             important details are not missed.  if this is not the
c             case, a large value, such as the interval length, is
c             suggested.
c
c    work
c    lenw   = (input)
c             work is an array of lenw real words used
c             internally for temporary storage.  the user must allocate
c             space for this array in the calling program by a statement
c             such as
c                       real work(...)
c             the following table gives the required minimum value for
c             the length of work, depending on the value of impl and
c             miter.  lenw should be set to the value used.  the
c             contents of work should not be disturbed between calls to
c             sdriv3.
c
c      impl =   0            1               2             3
c              ---------------------------------------------------------
c miter =  0   (mxord+4)*n   not allowed     not allowed   not allowed
c              + 2*nroot
c              + 250
c
c         1,2  n*n +         2*n*n +         n*n +         n*(n + nde)
c              (mxord+5)*n   (mxord+5)*n     (mxord+6)*n   + (mxord+5)*n
c              + 2*nroot     + 2*nroot       + 2*nroot     + 2*nroot
c              + 250         + 250           + 250         + 250
c
c          3   (mxord+4)*n   (mxord+4)*n     (mxord+4)*n   (mxord+4)*n
c              + 2*nroot     + 2*nroot       + 2*nroot     + 2*nroot
c              + 250         + 250           + 250         + 250
c
c         4,5  (2*ml+mu+1)   2*(2*ml+mu+1)   (2*ml+mu+1)   (2*ml+mu+1)*
c              *n +          *n +            *n +          (n+nde) +
c              (mxord+5)*n   (mxord+5)*n     (mxord+6)*n   + (mxord+5)*n
c              + 2*nroot     + 2*nroot       + 2*nroot     + 2*nroot
c              + 250         + 250           + 250         + 250
c              ---------------------------------------------------------
c
c    iwork
c    leniw  = (input)
c             iwork is an integer array of length leniw used internally
c             for temporary storage.  the user must allocate space for
c             this array in the calling program by a statement such as
c                       integer iwork(...)
c             the length of iwork should be at least
c               50      if miter is 0 or 3, or
c               n+50    if miter is 1, 2, 4, or 5, or mint is 3,
c             and leniw should be set to the value used.  the contents
c             of iwork should not be disturbed between calls to sdriv3.
c
c    jacobn = a subroutine supplied by the user, if miter is 1 or 4.
c             if this is the case, the name must be declared external in
c             the user's calling program.  given a system of n
c             differential equations, it is meaningful to speak about
c             the partial derivative of the i-th right hand side with
c             respect to the j-th dependent variable.  in general there
c             are n*n such quantities.  often however the equations can
c             be ordered so that the i-th differential equation only
c             involves dependent variables with index near i, e.g., i+1,
c             i-2.  such a system is called banded.  if, for all i, the
c             i-th equation depends on at most the variables
c               y(i-ml), y(i-ml+1), ... , y(i), y(i+1), ... , y(i+mu)
c             then we call ml+mu+1 the bandwidth of the system.  in a
c             banded system many of the partial derivatives above are
c             automatically zero.  for the cases miter = 1, 2, 4, and 5,
c             some of these partials are needed.  for the cases
c             miter = 2 and 5 the necessary derivatives are
c             approximated numerically by sdriv3, and we only ask the
c             user to tell sdriv3 the value of ml and mu if the system
c             is banded.  for the cases miter = 1 and 4 the user must
c             derive these partials algebraically and encode them in
c             subroutine jacobn.  by computing these derivatives the
c             user can often save 20-30 per cent of the computing time.
c             usually, however, the accuracy is not much affected and
c             most users will probably forego this option.  the optional
c             user-written subroutine jacobn has the form:
c                   subroutine jacobn (n, t, y, dfdy, matdim, ml, mu)
c                   real y(*), dfdy(matdim,*)
c                     .
c                     .
c                     calculate values of dfdy
c                     .
c                     .
c                   end (sample)
c             here y is a vector of length at least n.  the actual
c             length of y is determined by the user's declaration in the
c             program which calls sdriv3.  thus the dimensioning of y in
c             jacobn, while required by fortran convention, does not
c             actually allocate any storage.  when this subroutine is
c             called, the first n components of y are intermediate
c             approximations to the solution components.  the user
c             should not alter these values.  if the system is not
c             banded (miter=1), the partials of the i-th equation with
c             respect to the j-th dependent function are to be stored in
c             dfdy(i,j).  thus partials of the i-th equation are stored
c             in the i-th row of dfdy.  if the system is banded
c             (miter=4), then the partials of the i-th equation with
c             respect to y(j) are to be stored in dfdy(k,j), where
c             k=i-j+mu+1 .  normally a return from jacobn passes control
c             back to sdriv3.  however, if the user would like to abort
c             the calculation, i.e., return control to the program which
c             calls sdriv3, he should set n to zero.  sdriv3 will signal
c             this by returning a value of nstate equal to +8(-8).
c             altering the value of n in jacobn has no effect on the
c             value of n in the call sequence of sdriv3.
c
c    fa     = a subroutine supplied by the user if impl is not zero, and
c             miter is not 3.  if so, the name must be declared external
c             in the user's calling program.  this subroutine computes
c             the array a, where a*dy(i)/dt = f(y(i),t).
c             there are three cases:
c
c               impl=1.
c               subroutine fa is of the form:
c                   subroutine fa (n, t, y, a, matdim, ml, mu, nde)
c                   real y(*), a(matdim,*)
c                     .
c                     .
c                     calculate all values of a
c                     .
c                     .
c                   end (sample)
c               in this case a is assumed to be a nonsingular matrix,
c               with the same structure as dfdy (see jacobn description
c               above).  programming considerations prevent complete
c               generality.  if miter is 1 or 2, a is assumed to be full
c               and the user must compute and store all values of
c               a(i,j), i,j=1, ... ,n.  if miter is 4 or 5, a is assumed
c               to be banded with lower and upper half bandwidth ml and
c               mu.  the left hand side of the i-th equation is a linear
c               combination of dy(i-ml)/dt, dy(i-ml+1)/dt, ... ,
c               dy(i)/dt, ... , dy(i+mu-1)/dt, dy(i+mu)/dt.  thus in the
c               i-th equation, the coefficient of dy(j)/dt is to be
c               stored in a(k,j), where k=i-j+mu+1.
c               note: the array a will be altered between calls to fa.
c
c               impl=2.
c               subroutine fa is of the form:
c                   subroutine fa (n, t, y, a, matdim, ml, mu, nde)
c                   real y(*), a(*)
c                     .
c                     .
c                     calculate non-zero values of a(1),...,a(nde)
c                     .
c                     .
c                   end (sample)
c               in this case it is assumed that the system is ordered by
c               the user so that the differential equations appear
c               first, and the algebraic equations appear last.  the
c               algebraic equations must be written in the form:
c               0 = f(y(i),t).  when using this option it is up to the
c               user to provide initial values for the y(i) that satisfy
c               the algebraic equations as well as possible.  it is
c               further assumed that a is a vector of length nde.  all
c               of the components of a, which may depend on t, y(i),
c               etc., must be set by the user to non-zero values.
c
c               impl=3.
c               subroutine fa is of the form:
c                   subroutine fa (n, t, y, a, matdim, ml, mu, nde)
c                   real y(*), a(matdim,*)
c                     .
c                     .
c                     calculate all values of a
c                     .
c                     .
c                   end (sample)
c               in this case a is assumed to be a nonsingular nde by nde
c               matrix with the same structure as dfdy (see jacobn
c               description above).  programming considerations prevent
c               complete generality.  if miter is 1 or 2, a is assumed
c               to be full and the user must compute and store all
c               values of a(i,j), i,j=1, ... ,nde.  if miter is 4 or 5,
c               a is assumed to be banded with lower and upper half
c               bandwidths ml and mu.  the left hand side of the i-th
c               equation is a linear combination of dy(i-ml)/dt,
c               dy(i-ml+1)/dt, ... , dy(i)/dt, ... , dy(i+mu-1)/dt,
c               dy(i+mu)/dt.  thus in the i-th equation, the coefficient
c               of dy(j)/dt is to be stored in a(k,j), where k=i-j+mu+1.
c               it is assumed that the system is ordered by the user so
c               that the differential equations appear first, and the
c               algebraic equations appear last.  the algebraic
c               equations must be written in the form 0 = f(y(i),t).
c               when using this option it is up to the user to provide
c               initial values for the y(i) that satisfy the algebraic
c               equations as well as possible.
c               note: for impl = 3, the array a will be altered between
c               calls to fa.
c             here y is a vector of length at least n.  the actual
c             length of y is determined by the user's declaration in the
c             program which calls sdriv3.  thus the dimensioning of y in
c             fa, while required by fortran convention, does not
c             actually allocate any storage.  when this subroutine is
c             called, the first n components of y are intermediate
c             approximations to the solution components.  the user
c             should not alter these values.  fa is always called
c             immediately after calling f, with the same values of t
c             and y.  normally a return from fa passes control back to
c             sdriv3.  however, if the user would like to abort the
c             calculation, i.e., return control to the program which
c             calls sdriv3, he should set n to zero.  sdriv3 will signal
c             this by returning a value of nstate equal to +9(-9).
c             altering the value of n in fa has no effect on the value
c             of n in the call sequence of sdriv3.
c
c    nde    = (input) the number of differential equations.  this is
c             required only for impl = 2 or 3, with nde .lt. n.
c
c    mxstep = (input) the maximum number of internal steps allowed on
c             one call to sdriv3.
c
c    g      = a real fortran function supplied by the user
c             if nroot is not 0.  in this case, the name must be
c             declared external in the user's calling program.  g is
c             repeatedly called with different values of iroot to obtain
c             the value of each of the nroot equations for which a root
c             is desired.  g is of the form:
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
c             program which calls sdriv3.  thus the dimensioning of y in
c             g, while required by fortran convention, does not actually
c             allocate any storage.  normally a return from g passes
c             control back to  sdriv3.  however, if the user would like
c             to abort the calculation, i.e., return control to the
c             program which calls sdriv3, he should set n to zero.
c             sdriv3 will signal this by returning a value of nstate
c             equal to +7(-7).  in this case, the index of the equation
c             being evaluated is stored in the sixth element of iwork.
c             altering the value of n in g has no effect on the value of
c             n in the call sequence of sdriv3.
c
c    users  = a subroutine supplied by the user, if miter is 3.
c             if this is the case, the name must be declared external in
c             the user's calling program.  the routine users is called
c             by sdriv3 when certain linear systems must be solved.  the
c             user may choose any method to form, store and solve these
c             systems in order to obtain the solution result that is
c             returned to sdriv3.  in particular, this allows sparse
c             matrix methods to be used.  the call sequence for this
c             routine is:
c
c                subroutine users (y, yh, ywt, save1, save2, t, h, el,
c               8                  impl, n, nde, iflag)
c                real y(*), yh(*), ywt(*), save1(*),
c               8     save2(*), t, h, el
c
c             the input variable iflag indicates what action is to be
c             taken.  subroutine users should perform the following
c             operations, depending on the value of iflag and impl.
c
c               iflag = 0
c                 impl = 0.  users is not called.
c                 impl = 1, 2 or 3.  solve the system a*x = save2,
c                   returning the result in save2.  the array save1 can
c                   be used as a work array.  for impl = 1, there are n
c                   components to the system, and for impl = 2 or 3,
c                   there are nde components to the system.
c
c               iflag = 1
c                 impl = 0.  compute, decompose and store the matrix
c                   (i - h*el*j), where i is the identity matrix and j
c                   is the jacobian matrix of the right hand side.  the
c                   array save1 can be used as a work array.
c                 impl = 1, 2 or 3. compute, decompose and store the
c                   matrix (a - h*el*j).  the array save1 can be used as
c                   a work array.
c
c               iflag = 2
c                 impl = 0.   solve the system
c                     (i - h*el*j)*x = h*save2 - yh - save1,
c                   returning the result in save2.
c                 impl = 1, 2 or 3.  solve the system
c                   (a - h*el*j)*x = h*save2 - a*(yh + save1)
c                   returning the result in save2.
c                 the array save1 should not be altered.
c             if iflag is 0 and impl is 1 or 2 and the matrix a is
c             singular, or if iflag is 1 and one of the matrices
c             (i - h*el*j), (a - h*el*j) is singular, the integer
c             variable iflag is to be set to -1 before returning.
c             normally a return from users passes control back to
c             sdriv3.  however, if the user would like to abort the
c             calculation, i.e., return control to the program which
c             calls sdriv3, he should set n to zero.  sdriv3 will signal
c             this by returning a value of nstate equal to +10(-10).
c             altering the value of n in users has no effect on the
c             value of n in the call sequence of sdriv3.
c
c    ierflg = an error flag.  the error number associated with a
c             diagnostic message (see section iii-a below) is the same
c             as the corresponding value of ierflg.  the meaning of
c             ierflg:
c               0  the routine completed successfully. (no message is
c                  issued.)
c               3  (warning) the number of steps required to reach tout
c                  exceeds mxstep.
c               4  (warning) the value of eps is too small.
c              11  (warning) for ntask = 2 or 3, t is beyond tout.
c                  the solution was obtained by interpolation.
c              15  (warning) the integration step size is below the
c                  roundoff level of t.  (the program issues this
c                  message as a warning but does not return control to
c                  the user.)
c              22  (recoverable) n is not positive.
c              23  (recoverable) mint is less than 1 or greater than 3 .
c              24  (recoverable) miter is less than 0 or greater than
c                  5 .
c              25  (recoverable) impl is less than 0 or greater than 3 .
c              26  (recoverable) the value of nstate is less than 1 or
c                  greater than 12 .
c              27  (recoverable) eps is less than zero.
c              28  (recoverable) mxord is not positive.
c              29  (recoverable) for mint = 3, either miter = 0 or 3, or
c                  impl = 0 .
c              30  (recoverable) for miter = 0, impl is not 0 .
c              31  (recoverable) for mint = 1, impl is 2 or 3 .
c              32  (recoverable) insufficient storage has been allocated
c                  for the work array.
c              33  (recoverable) insufficient storage has been allocated
c                  for the iwork array.
c              41  (recoverable) the integration step size has gone
c                  to zero.
c              42  (recoverable) the integration step size has been
c                  reduced about 50 times without advancing the
c                  solution.  the problem setup may not be correct.
c              43  (recoverable)  for impl greater than 0, the matrix a
c                  is singular.
c             999  (fatal) the value of nstate is 12 .
c
c  iii.  other communication to the user  ..............................
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
c  iv.  remarks  .......................................................
c
c    a. other routines used:
c         sdntp, sdzro, sdstp, sdntl, sdpst, sdcor, sdcst,
c         sdpsc, and sdscl;
c         sgefa, sgesl, sgbfa, sgbsl, and snrm2 (from linpack)
c         r1mach (from the bell laboratories machine constants package)
c         xermsg (from the slatec common math library)
c       the last seven routines above, not having been written by the
c       present authors, are not explicitly part of this package.
c
c    b. on any return from sdriv3 all information necessary to continue
c       the calculation is contained in the call sequence parameters,
c       including the work arrays.  thus it is possible to suspend one
c       problem, integrate another, and then return to the first.
c
c    c. if this package is to be used in an overlay situation, the user
c       must declare in the primary overlay the variables in the call
c       sequence to sdriv3.
c
c    d. changing parameters during an integration.
c       the value of nroot, eps, ewt, ierror, mint, miter, or hmax may
c       be altered by the user between calls to sdriv3.  for example, if
c       too much accuracy has been requested (the program returns with
c       nstate = 4 and an increased value of eps) the user may wish to
c       increase eps further.  in general, prudence is necessary when
c       making changes in parameters since such changes are not
c       implemented until the next integration step, which is not
c       necessarily the next call to sdriv3.  this can happen if the
c       program has already integrated to a point which is beyond the
c       new point tout.
c
c    e. as the price for complete control of matrix algebra, the sdriv3
c       users option puts all responsibility for jacobian matrix
c       evaluation on the user.  it is often useful to approximate
c       numerically all or part of the jacobian matrix.  however this
c       must be done carefully.  the fortran sequence below illustrates
c       the method we recommend.  it can be inserted directly into
c       subroutine users to approximate jacobian elements in rows i1
c       to i2 and columns j1 to j2.
c              real dfdy(n,n), epsj, h, r, r1mach,
c             8     save1(n), save2(n), t, uround, y(n), yj, ywt(n)
c              uround = r1mach(4)
c              epsj = sqrt(uround)
c              do 30 j = j1,j2
c                r = epsj*max(abs(ywt(j)), abs(y(j)))
c                if (r .eq. 0.e0) r = ywt(j)
c                yj = y(j)
c                y(j) = y(j) + r
c                call f (n, t, y, save1)
c                if (n .eq. 0) return
c                y(j) = yj
c                do 20 i = i1,i2
c         20       dfdy(i,j) = (save1(i) - save2(i))/r
c         30     continue
c       many problems give rise to structured sparse jacobians, e.g.,
c       block banded.  it is possible to approximate them with fewer
c       function evaluations than the above procedure uses; see curtis,
c       powell and reid, j. inst. maths applics, (1974), vol. 13,
c       pp. 117-119.
c
c    f. when any of the routines jacobn, fa, g, or users, is not
c       required, difficulties associated with unsatisfied externals can
c       be avoided by using the name of the routine which calculates the
c       right hand side of the differential equations in place of the
c       corresponding name in the call sequence of sdriv3.
c
c***references  c. w. gear, numerical initial value problems in
c                 ordinary differential equations, prentice-hall, 1971.
c***routines called  r1mach, sdntp, sdstp, sdzro, sgbfa, sgbsl, sgefa,
c                    sgesl, snrm2, xermsg
c***revision history  (yymmdd)
c   790601  date written
c   900329  initial submission to slatec.
c***end prologue  sdriv3
      external f, jacobn, fa, g, users
      real ae, big, eps, ewt(*), g, glast, gnow, h, hmax,
     8     hsign, hused, nround, re, r1mach, size, snrm2, sum, t, tlast,
     8     tout, troot, uround, work(*), y(*)
      integer i, ia, iavgh, iavgrd, icnvrg, idfdy, iel, ierflg, ierror,
     8        ifac, iflag, ignow, ih, ihmax, ihold, ihsign, ihused,
     8        ijroot, ijstpl, ijtask, imnt, imntld, impl, imtr, imtrld,
     8        imtrsv, imxerr, imxord, imxrds, indmxr, indprt, indpvt,
     8        indtrt, infe, info, inje, inq, inquse, inroot, inrtld,
     8        instep, inwait, irc, irmax, iroot, imach1, imach4, isave1,
     8        isave2, it, itout, itq, itrend, itroot, iwork(*), iyh,
     8        iywt, j, jstate, jtroot, lenchk, leniw, lenw, liwchk,
     8        matdim, maxord, mint, miter, ml, mu, mxord, mxstep, n,
     8        nde, ndecom, npar, nroot, nstate, nstepl, ntask
      logical convrg
      character intgr1*8, intgr2*8, rl1*16, rl2*16
      parameter(nround = 20.e0)
      parameter(iavgh = 1, ihused = 2, iavgrd = 3,
     8          iel = 4, ih = 160, ihmax = 161, ihold = 162,
     8          ihsign = 163, irc = 164, irmax = 165, it = 166,
     8          itout = 167, itq = 168, itrend = 204, imach1 = 205,
     8          imach4 = 206, iyh = 251,
     8          indmxr = 1, inquse = 2, instep = 3, infe = 4, inje = 5,
     8          inroot = 6, icnvrg = 7, ijroot = 8, ijtask = 9,
     8          imntld = 10, imtrld = 11, inq = 12, inrtld = 13,
     8          indtrt = 14, inwait = 15, imnt = 16, imtrsv = 17,
     8          imtr = 18, imxrds = 19, imxord = 20, indprt = 21,
     8          ijstpl = 22, indpvt = 51)
c***first executable statement  sdriv3
      if (nstate .eq. 12) then
        ierflg = 999
        call xermsg('slatec', 'sdriv3',
     8  'illegal input.  the value of nstate is 12 .', ierflg, 2)
        return
      else if (nstate .lt. 1 .or. nstate .gt. 12) then
        write(intgr1, '(i8)') nstate
        ierflg = 26
        call xermsg('slatec', 'sdriv3',
     8  'illegal input.  improper value for nstate(= '//intgr1//').',
     8  ierflg, 1)
        nstate = 12
        return
      end if
      npar = n
      if (eps .lt. 0.e0) then
        write(rl1, '(e16.8)') eps
        ierflg = 27
        call xermsg('slatec', 'sdriv3',
     8  'illegal input.  eps, '//rl1//', is negative.', ierflg, 1)
        nstate = 12
        return
      end if
      if (n .le. 0) then
        write(intgr1, '(i8)') n
        ierflg = 22
        call xermsg('slatec', 'sdriv3',
     8  'illegal input.  number of equations, '//intgr1//
     8  ', is not positive.', ierflg, 1)
        nstate = 12
        return
      end if
      if (mxord .le. 0) then
        write(intgr1, '(i8)') mxord
        ierflg = 28
        call xermsg('slatec', 'sdriv3',
     8  'illegal input.  maximum order, '//intgr1//
     8  ', is not positive.', ierflg, 1)
        nstate = 12
        return
      end if
      if (mint .lt. 1 .or. mint .gt. 3) then
        write(intgr1, '(i8)') mint
        ierflg = 23
        call xermsg('slatec', 'sdriv3',
     8  'illegal input.  improper value for the integration method '//
     8  'flag, '//intgr1//' .', ierflg, 1)
        nstate = 12
        return
      else if (miter .lt. 0 .or. miter .gt. 5) then
        write(intgr1, '(i8)') miter
        ierflg = 24
        call xermsg('slatec', 'sdriv3',
     8  'illegal input.  improper value for miter(= '//intgr1//').',
     8  ierflg, 1)
        nstate = 12
        return
      else if (impl .lt. 0 .or. impl .gt. 3) then
        write(intgr1, '(i8)') impl
        ierflg = 25
        call xermsg('slatec', 'sdriv3',
     8  'illegal input.  improper value for impl(= '//intgr1//').',
     8  ierflg, 1)
        nstate = 12
        return
      else if (mint .eq. 3 .and.
     8  (miter .eq. 0 .or. miter .eq. 3 .or. impl .ne. 0)) then
        write(intgr1, '(i8)') miter
        write(intgr2, '(i8)') impl
        ierflg = 29
        call xermsg('slatec', 'sdriv3',
     8  'illegal input.  for mint = 3, the value of miter, '//intgr1//
     8  ', and/or impl, '//intgr2//', is not allowed.', ierflg, 1)
        nstate = 12
        return
      else if ((impl .ge. 1 .and. impl .le. 3) .and. miter .eq. 0) then
        write(intgr1, '(i8)') impl
        ierflg = 30
        call xermsg('slatec', 'sdriv3',
     8  'illegal input.  for miter = 0, the value of impl, '//intgr1//
     8  ', is not allowed.', ierflg, 1)
        nstate = 12
        return
      else if ((impl .eq. 2 .or. impl .eq. 3) .and. mint .eq. 1) then
        write(intgr1, '(i8)') impl
        ierflg = 31
        call xermsg('slatec', 'sdriv3',
     8  'illegal input.  for mint = 1, the value of impl, '//intgr1//
     8  ', is not allowed.', ierflg, 1)
        nstate = 12
        return
      end if
      if (miter .eq. 0 .or. miter .eq. 3) then
        liwchk = indpvt - 1
      else if (miter .eq. 1 .or. miter .eq. 2 .or. miter .eq. 4 .or.
     8  miter .eq. 5) then
        liwchk = indpvt + n - 1
      end if
      if (leniw .lt. liwchk) then
        write(intgr1, '(i8)') liwchk
        ierflg = 33
        call xermsg('slatec', 'sdriv3',
     8  'illegal input.  insufficient storage allocated for the '//
     8  'iwork array.  based on the value of the input parameters '//
     8  'involved, the required storage is '//intgr1//' .', ierflg, 1)
        nstate = 12
        return
      end if
c                                                allocate the work array
c                                         iyh is the index of yh in work
      if (mint .eq. 1 .or. mint .eq. 3) then
        maxord = min(mxord, 12)
      else if (mint .eq. 2) then
        maxord = min(mxord, 5)
      end if
      idfdy = iyh + (maxord + 1)*n
c                                             idfdy is the index of dfdy
c
      if (miter .eq. 0 .or. miter .eq. 3) then
        iywt = idfdy
      else if (miter .eq. 1 .or. miter .eq. 2) then
        iywt = idfdy + n*n
      else if (miter .eq. 4 .or. miter .eq. 5) then
        iywt = idfdy + (2*ml + mu + 1)*n
      end if
c                                               iywt is the index of ywt
      isave1 = iywt + n
c                                           isave1 is the index of save1
      isave2 = isave1 + n
c                                           isave2 is the index of save2
      ignow = isave2 + n
c                                             ignow is the index of gnow
      itroot = ignow + nroot
c                                           itroot is the index of troot
      ifac = itroot + nroot
c                                               ifac is the index of fac
      if (miter .eq. 2 .or. miter .eq. 5 .or. mint .eq. 3) then
        ia = ifac + n
      else
        ia = ifac
      end if
c                                                   ia is the index of a
      if (impl .eq. 0 .or. miter .eq. 3) then
        lenchk = ia - 1
      else if (impl .eq. 1 .and. (miter .eq. 1 .or. miter .eq. 2)) then
        lenchk = ia - 1 + n*n
      else if (impl .eq. 1 .and. (miter .eq. 4 .or. miter .eq. 5)) then
        lenchk = ia - 1 + (2*ml + mu + 1)*n
      else if (impl .eq. 2 .and. miter .ne. 3) then
        lenchk = ia - 1 + n
      else if (impl .eq. 3 .and. (miter .eq. 1 .or. miter .eq. 2)) then
        lenchk = ia - 1 + n*nde
      else if (impl .eq. 3 .and. (miter .eq. 4 .or. miter .eq. 5)) then
        lenchk = ia - 1 + (2*ml + mu + 1)*nde
      end if
      if (lenw .lt. lenchk) then
        write(intgr1, '(i8)') lenchk
        ierflg = 32
        call xermsg('slatec', 'sdriv3',
     8  'illegal input.  insufficient storage allocated for the '//
     8  'work array.  based on the value of the input parameters '//
     8  'involved, the required storage is '//intgr1//' .', ierflg, 1)
        nstate = 12
        return
      end if
      if (miter .eq. 0 .or. miter .eq. 3) then
        matdim = 1
      else if (miter .eq. 1 .or. miter .eq. 2) then
        matdim = n
      else if (miter .eq. 4 .or. miter .eq. 5) then
        matdim = 2*ml + mu + 1
      end if
      if (impl .eq. 0 .or. impl .eq. 1) then
        ndecom = n
      else if (impl .eq. 2 .or. impl .eq. 3) then
        ndecom = nde
      end if
      if (nstate .eq. 1) then
c                                                  initialize parameters
        if (mint .eq. 1 .or. mint .eq. 3) then
          iwork(imxord) = min(mxord, 12)
        else if (mint .eq. 2) then
          iwork(imxord) = min(mxord, 5)
        end if
        iwork(imxrds) = mxord
        if (mint .eq. 1 .or. mint .eq. 2) then
          iwork(imnt) = mint
          iwork(imtr) = miter
          iwork(imntld) = mint
          iwork(imtrld) = miter
        else if (mint .eq. 3) then
          iwork(imnt) = 1
          iwork(imtr) = 0
          iwork(imntld) = iwork(imnt)
          iwork(imtrld) = iwork(imtr)
          iwork(imtrsv) = miter
        end if
        work(ihmax) = hmax
        uround = r1mach (4)
        work(imach4) = uround
        work(imach1) = r1mach (1)
        if (nroot .ne. 0) then
          re = uround
          ae = work(imach1)
        end if
        h = (tout - t)*(1.e0 - 4.e0*uround)
        h = sign(min(abs(h), hmax), h)
        work(ih) = h
        hsign = sign(1.e0, h)
        work(ihsign) = hsign
        iwork(ijtask) = 0
        work(iavgh) = 0.e0
        work(ihused) = 0.e0
        work(iavgrd) = 0.e0
        iwork(indmxr) = 0
        iwork(inquse) = 0
        iwork(instep) = 0
        iwork(ijstpl) = 0
        iwork(infe) = 0
        iwork(inje) = 0
        iwork(inroot) = 0
        work(it) = t
        iwork(icnvrg) = 0
        iwork(indprt) = 0
c                                                 set initial conditions
        do 30 i = 1,n
 30       work(i+iyh-1) = y(i)
        if (t .eq. tout) return
        go to 180
      else
        uround = work(imach4)
        if (nroot .ne. 0) then
          re = uround
          ae = work(imach1)
        end if
      end if
c                                             on a continuation, check
c                                             that output points have
c                                             been or will be overtaken.
      if (iwork(icnvrg) .eq. 1) then
        convrg = .true.
      else
        convrg = .false.
      end if
      t = work(it)
      h = work(ih)
      hsign = work(ihsign)
      if (iwork(ijtask) .eq. 0) go to 180
c
c                                   iwork(ijroot) flags unreported
c                                   roots, and is set to the value of
c                                   ntask when a root was last selected.
c                                   it is set to zero when all roots
c                                   have been reported.  iwork(inroot)
c                                   contains the index and work(itout)
c                                   contains the value of the root last
c                                   selected to be reported.
c                                   iwork(inrtld) contains the value of
c                                   nroot and iwork(indtrt) contains
c                                   the value of itroot when the array
c                                   of roots was last calculated.
      if (nroot .ne. 0) then
        if (iwork(ijroot) .gt. 0) then
c                                      tout has just been reported.
c                                      if troot .le. tout, report troot.
          if (nstate .ne. 5) then
            if (tout*hsign .ge. work(itout)*hsign) then
              troot = work(itout)
              call sdntp (h, 0, n, iwork(inq), t, troot, work(iyh),  y)
              t = troot
              nstate = 5
              ierflg = 0
              go to 580
            end if
c                                         a root has just been reported.
c                                         select the next root.
          else
            troot = t
            iroot = 0
            do 50 i = 1,iwork(inrtld)
              jtroot = i + iwork(indtrt) - 1
              if (work(jtroot)*hsign .le. troot*hsign) then
c
c                                              check for multiple roots.
c
                if (work(jtroot) .eq. work(itout) .and.
     8          i .gt. iwork(inroot)) then
                  iroot = i
                  troot = work(jtroot)
                  go to 60
                end if
                if (work(jtroot)*hsign .gt. work(itout)*hsign) then
                  iroot = i
                  troot = work(jtroot)
                end if
              end if
 50           continue
 60         iwork(inroot) = iroot
            work(itout) = troot
            iwork(ijroot) = ntask
            if (ntask .eq. 1) then
              if (iroot .eq. 0) then
                iwork(ijroot) = 0
              else
                if (tout*hsign .ge. troot*hsign) then
                  call sdntp (h, 0, n, iwork(inq), t, troot, work(iyh),
     8                        y)
                  nstate = 5
                  t = troot
                  ierflg = 0
                  go to 580
                end if
              end if
            else if (ntask .eq. 2 .or. ntask .eq. 3) then
c
c                                     if there are no more roots, or the
c                                     user has altered tout to be less
c                                     than a root, set ijroot to zero.
c
              if (iroot .eq. 0 .or. (tout*hsign .lt. troot*hsign)) then
                iwork(ijroot) = 0
              else
                call sdntp (h, 0, n, iwork(inq), t, troot, work(iyh),
     8                      y)
                nstate = 5
                ierflg = 0
                t = troot
                go to 580
              end if
            end if
          end if
        end if
      end if
c
      if (ntask .eq. 1) then
        nstate = 2
        if (t*hsign .ge. tout*hsign) then
          call sdntp (h, 0, n, iwork(inq), t, tout, work(iyh),  y)
          t = tout
          ierflg = 0
          go to 580
        end if
      else if (ntask .eq. 2) then
c                                                      check if tout has
c                                                      been reset .lt. t
        if (t*hsign .gt. tout*hsign) then
          write(rl1, '(e16.8)') t
          write(rl2, '(e16.8)') tout
          ierflg = 11
          call xermsg('slatec', 'sdriv3',
     8    'while integrating exactly to tout, t, '//rl1//
     8    ', was beyond tout, '//rl2//' .  solution obtained by '//
     8    'interpolation.', ierflg, 0)
          nstate = 11
          call sdntp (h, 0, n, iwork(inq), t, tout, work(iyh),  y)
          t = tout
          go to 580
        end if
c                                   determine if tout has been overtaken
c
        if (abs(tout - t).le.nround*uround*max(abs(t), abs(tout))) then
          t = tout
          nstate = 2
          ierflg = 0
          go to 560
        end if
c                                             if there are no more roots
c                                             to report, report t.
        if (nstate .eq. 5) then
          nstate = 2
          ierflg = 0
          go to 560
        end if
        nstate = 2
c                                                       see if tout will
c                                                       be overtaken.
        if ((t + h)*hsign .gt. tout*hsign) then
          h = tout - t
          if ((t + h)*hsign .gt. tout*hsign) h = h*(1.e0 - 4.e0*uround)
          work(ih) = h
          if (h .eq. 0.e0) go to 670
          iwork(ijtask) = -1
        end if
      else if (ntask .eq. 3) then
        nstate = 2
        if (t*hsign .gt. tout*hsign) then
          write(rl1, '(e16.8)') t
          write(rl2, '(e16.8)') tout
          ierflg = 11
          call xermsg('slatec', 'sdriv3',
     8    'while integrating exactly to tout, t, '//rl1//
     8    ', was beyond tout, '//rl2//' .  solution obtained by '//
     8    'interpolation.', ierflg, 0)
          nstate = 11
          call sdntp (h, 0, n, iwork(inq), t, tout, work(iyh),  y)
          t = tout
          go to 580
        end if
        if (abs(tout - t).le.nround*uround*max(abs(t), abs(tout))) then
          t = tout
          ierflg = 0
          go to 560
        end if
        if ((t + h)*hsign .gt. tout*hsign) then
          h = tout - t
          if ((t + h)*hsign .gt. tout*hsign) h = h*(1.e0 - 4.e0*uround)
          work(ih) = h
          if (h .eq. 0.e0) go to 670
          iwork(ijtask) = -1
        end if
      end if
c                         implement changes in mint, miter, and/or hmax.
c
      if ((mint .ne. iwork(imntld) .or. miter .ne. iwork(imtrld)) .and.
     8  mint .ne. 3 .and. iwork(imntld) .ne. 3) iwork(ijtask) = -1
      if (hmax .ne. work(ihmax)) then
        h = sign(min(abs(h), hmax), h)
        if (h .ne. work(ih)) then
          iwork(ijtask) = -1
          work(ih) = h
        end if
        work(ihmax) = hmax
      end if
c
 180  nstepl = iwork(instep)
      do 190 i = 1,n
 190    y(i) = work(i+iyh-1)
      if (nroot .ne. 0) then
        do 200 i = 1,nroot
          work(i+ignow-1) = g (npar, t, y, i)
          if (npar .eq. 0) then
            iwork(inroot) = i
            nstate = 7
            return
          end if
 200     continue
      end if
      if (ierror .eq. 1) then
        do 230 i = 1,n
 230      work(i+iywt-1) = 1.e0
        go to 410
      else if (ierror .eq. 5) then
        do 250 i = 1,n
 250      work(i+iywt-1) = ewt(i)
        go to 410
      end if
c                                       reset ywt array.  looping point.
 260  if (ierror .eq. 2) then
        do 280 i = 1,n
          if (y(i) .eq. 0.e0) go to 290
 280      work(i+iywt-1) = abs(y(i))
        go to 410
 290    if (iwork(ijtask) .eq. 0) then
          call f (npar, t, y, work(isave2))
          if (npar .eq. 0) then
            nstate = 6
            return
          end if
          iwork(infe) = iwork(infe) + 1
          if (miter .eq. 3 .and. impl .ne. 0) then
            iflag = 0
            call users (y, work(iyh), work(iywt), work(isave1),
     8                  work(isave2), t, h, work(iel), impl, npar,
     8                  ndecom, iflag)
            if (iflag .eq. -1) go to 690
            if (npar .eq. 0) then
              nstate = 10
              return
            end if
          else if (impl .eq. 1) then
            if (miter .eq. 1 .or. miter .eq. 2) then
              call fa (npar, t, y, work(ia), matdim, ml, mu, ndecom)
              if (npar .eq. 0) then
                nstate = 9
                return
              end if
              call sgefa (work(ia), matdim, n, iwork(indpvt), info)
              if (info .ne. 0) go to 690
              call sgesl (work(ia), matdim, n, iwork(indpvt),
     8                    work(isave2), 0)
            else if (miter .eq. 4 .or. miter .eq. 5) then
              call fa (npar, t, y, work(ia+ml), matdim, ml, mu, ndecom)
              if (npar .eq. 0) then
                nstate = 9
                return
              end if
              call sgbfa (work(ia), matdim, n, ml, mu, iwork(indpvt),
     8                    info)
              if (info .ne. 0) go to 690
              call sgbsl (work(ia), matdim, n, ml, mu, iwork(indpvt),
     8                    work(isave2), 0)
            end if
          else if (impl .eq. 2) then
            call fa (npar, t, y, work(ia), matdim, ml, mu, ndecom)
            if (npar .eq. 0) then
              nstate = 9
              return
            end if
            do 340 i = 1,ndecom
              if (work(i+ia-1) .eq. 0.e0) go to 690
 340          work(i+isave2-1) = work(i+isave2-1)/work(i+ia-1)
          else if (impl .eq. 3) then
            if (miter .eq. 1 .or. miter .eq. 2) then
              call fa (npar, t, y, work(ia), matdim, ml, mu, ndecom)
              if (npar .eq. 0) then
                nstate = 9
                return
              end if
              call sgefa (work(ia), matdim, nde, iwork(indpvt), info)
              if (info .ne. 0) go to 690
              call sgesl (work(ia), matdim, nde, iwork(indpvt),
     8                    work(isave2), 0)
            else if (miter .eq. 4 .or. miter .eq. 5) then
              call fa (npar, t, y, work(ia+ml), matdim, ml, mu, ndecom)
              if (npar .eq. 0) then
                nstate = 9
                return
              end if
              call sgbfa (work(ia), matdim, nde, ml, mu, iwork(indpvt),
     8                    info)
              if (info .ne. 0) go to 690
              call sgbsl (work(ia), matdim, nde, ml, mu, iwork(indpvt),
     8                    work(isave2), 0)
            end if
          end if
        end if
        do 360 j = i,n
          if (y(j) .ne. 0.e0) then
            work(j+iywt-1) = abs(y(j))
          else
            if (iwork(ijtask) .eq. 0) then
              work(j+iywt-1) = abs(h*work(j+isave2-1))
            else
              work(j+iywt-1) = abs(work(j+iyh+n-1))
            end if
          end if
          if (work(j+iywt-1) .eq. 0.e0) work(j+iywt-1) = uround
 360      continue
      else if (ierror .eq. 3) then
        do 380 i = 1,n
 380      work(i+iywt-1) = max(ewt(1), abs(y(i)))
      else if (ierror .eq. 4) then
        do 400 i = 1,n
 400      work(i+iywt-1) = max(ewt(i), abs(y(i)))
      end if
c
 410  do 420 i = 1,n
 420    work(i+isave2-1) = y(i)/work(i+iywt-1)
      sum = snrm2(n, work(isave2), 1)/sqrt(real(n))
      sum = max(1.e0, sum)
      if (eps .lt. sum*uround) then
        eps = sum*uround*(1.e0 + 10.e0*uround)
        write(rl1, '(e16.8)') t
        write(rl2, '(e16.8)') eps
        ierflg = 4
        call xermsg('slatec', 'sdriv3',
     8  'at t, '//rl1//', the requested accuracy, eps, was not '//
     8  'obtainable with the machine precision.  eps has been '//
     8  'increased to '//rl2//' .', ierflg, 0)
        nstate = 4
        go to 560
      end if
      if (abs(h) .ge. uround*abs(t)) then
        iwork(indprt) = 0
      else if (iwork(indprt) .eq. 0) then
        write(rl1, '(e16.8)') t
        write(rl2, '(e16.8)') h
        ierflg = 15
        call xermsg('slatec', 'sdriv3',
     8  'at t, '//rl1//', the step size, '//rl2//', is smaller '//
     8  'than the roundoff level of t.  this may occur if there is '//
     8  'an abrupt change in the right hand side of the '//
     8  'differential equations.', ierflg, 0)
        iwork(indprt) = 1
      end if
      if (ntask.ne.2) then
        if ((iwork(instep)-nstepl) .eq. mxstep) then
          write(rl1, '(e16.8)') t
          write(intgr1, '(i8)') mxstep
          write(rl2, '(e16.8)') tout
          ierflg = 3
          call xermsg('slatec', 'sdriv3',
     8    'at t, '//rl1//', '//intgr1//' steps have been taken '//
     8    'without reaching tout, '//rl2//' .', ierflg, 0)
          nstate = 3
          go to 560
        end if
      end if
c
c     call sdstp (eps, f, fa, hmax, impl, ierror, jacobn, matdim,
c    8            maxord, mint, miter, ml, mu, n, nde, ywt, uround,
c    8            users,  avgh, avgord, h, hused, jtask, mntold, mtrold,
c    8            nfe, nje, nqused, nstep, t, y, yh,  a, convrg,
c    8            dfdy, el, fac, hold, ipvt, jstate, jstepl, nq, nwait,
c    8            rc, rmax, save1, save2, tq, trend, iswflg, mtrsv,
c    8            mxrdsv)
c
      call sdstp (eps, f, fa, work(ihmax), impl, ierror, jacobn,
     8            matdim, iwork(imxord), iwork(imnt), iwork(imtr), ml,
     8            mu, npar, ndecom, work(iywt), uround, users,
     8            work(iavgh), work(iavgrd), work(ih), hused,
     8            iwork(ijtask), iwork(imntld), iwork(imtrld),
     8            iwork(infe), iwork(inje), iwork(inquse),
     8            iwork(instep), work(it), y, work(iyh), work(ia),
     8            convrg, work(idfdy), work(iel), work(ifac),
     8            work(ihold), iwork(indpvt), jstate, iwork(ijstpl),
     8            iwork(inq), iwork(inwait), work(irc), work(irmax),
     8            work(isave1), work(isave2), work(itq), work(itrend),
     8            mint, iwork(imtrsv), iwork(imxrds))
      t = work(it)
      h = work(ih)
      if (convrg) then
        iwork(icnvrg) = 1
      else
        iwork(icnvrg) = 0
      end if
      go to (470, 670, 680, 690, 690, 660, 660, 660, 660, 660), jstate
 470  iwork(ijtask) = 1
c                                 determine if a root has been overtaken
      if (nroot .ne. 0) then
        iroot = 0
        do 500 i = 1,nroot
          glast = work(i+ignow-1)
          gnow = g (npar, t, y, i)
          if (npar .eq. 0) then
            iwork(inroot) = i
            nstate = 7
            return
          end if
          work(i+ignow-1) = gnow
          if (glast*gnow .gt. 0.e0) then
            work(i+itroot-1) = t + h
          else
            if (gnow .eq. 0.e0) then
              work(i+itroot-1) = t
              iroot = i
            else
              if (glast .eq. 0.e0) then
                work(i+itroot-1) = t + h
              else
                if (abs(hused) .ge. uround*abs(t)) then
                  tlast = t - hused
                  iroot = i
                  troot = t
                  call sdzro (ae, g, h, npar, iwork(inq), iroot, re, t,
     8                        work(iyh), uround,  troot, tlast,
     8                        gnow, glast,  y)
                  do 480 j = 1,n
 480                y(j) = work(iyh+j-1)
                  if (npar .eq. 0) then
                    iwork(inroot) = i
                    nstate = 7
                    return
                  end if
                  work(i+itroot-1) = troot
                else
                  work(i+itroot-1) = t
                  iroot = i
                end if
              end if
            end if
          end if
 500      continue
        if (iroot .eq. 0) then
          iwork(ijroot) = 0
c                                                  select the first root
        else
          iwork(ijroot) = ntask
          iwork(inrtld) = nroot
          iwork(indtrt) = itroot
          troot = t + h
          do 510 i = 1,nroot
            if (work(i+itroot-1)*hsign .lt. troot*hsign) then
              troot = work(i+itroot-1)
              iroot = i
            end if
 510        continue
          iwork(inroot) = iroot
          work(itout) = troot
          if (troot*hsign .le. tout*hsign) then
            call sdntp (h, 0, n, iwork(inq), t, troot, work(iyh),  y)
            nstate = 5
            t = troot
            ierflg = 0
            go to 580
          end if
        end if
      end if
c                               test for ntask condition to be satisfied
      nstate = 2
      if (ntask .eq. 1) then
        if (t*hsign .lt. tout*hsign) go to 260
        call sdntp (h, 0, n, iwork(inq), t, tout, work(iyh),  y)
        t = tout
        ierflg = 0
        go to 580
c                               tout is assumed to have been attained
c                               exactly if t is within twenty roundoff
c                               units of tout, relative to max(tout, t).
c
      else if (ntask .eq. 2) then
        if (abs(tout - t).le.nround*uround*max(abs(t), abs(tout))) then
          t = tout
        else
          if ((t + h)*hsign .gt. tout*hsign) then
            h = tout - t
            if ((t + h)*hsign.gt.tout*hsign) h = h*(1.e0 - 4.e0*uround)
            work(ih) = h
            if (h .eq. 0.e0) go to 670
            iwork(ijtask) = -1
          end if
        end if
      else if (ntask .eq. 3) then
        if (abs(tout - t).le.nround*uround*max(abs(t), abs(tout))) then
          t = tout
        else
          if ((t + h)*hsign .gt. tout*hsign) then
            h = tout - t
            if ((t + h)*hsign.gt.tout*hsign) h = h*(1.e0 - 4.e0*uround)
            work(ih) = h
            if (h .eq. 0.e0) go to 670
            iwork(ijtask) = -1
          end if
          go to 260
        end if
      end if
      ierflg = 0
c                                      all returns are made through this
c                                      section.  imxerr is determined.
 560  do 570 i = 1,n
 570    y(i) = work(i+iyh-1)
 580  if (iwork(ijtask) .eq. 0) return
      big = 0.e0
      imxerr = 1
      do  590 i = 1,n
c                                            size = abs(error(i)/ywt(i))
        size = abs(work(i+isave1-1)/work(i+iywt-1))
        if (big .lt. size) then
          big = size
          imxerr = i
        end if
 590    continue
      iwork(indmxr) = imxerr
      work(ihused) = hused
      return
c
 660  nstate = jstate
      return
c                                        fatal errors are processed here
c
 670  write(rl1, '(e16.8)') t
      ierflg = 41
      call xermsg('slatec', 'sdriv3',
     8  'at t, '//rl1//', the attempted step size has gone to '//
     8  'zero.  often this occurs if the problem setup is incorrect.',
     8  ierflg, 1)
      nstate = 12
      return
c
 680  write(rl1, '(e16.8)') t
      ierflg = 42
      call xermsg('slatec', 'sdriv3',
     8  'at t, '//rl1//', the step size has been reduced about 50 '//
     8  'times without advancing the solution.  often this occurs '//
     8  'if the problem setup is incorrect.', ierflg, 1)
      nstate = 12
      return
c
 690  write(rl1, '(e16.8)') t
      ierflg = 43
      call xermsg('slatec', 'sdriv3',
     8  'at t, '//rl1//', while solving a*ydot = f, a is singular.',
     8  ierflg, 1)
      nstate = 12
      return
      end
