*deck debdf
      subroutine debdf (f, neq, t, y, tout, info, rtol, atol, idid,
     +   rwork, lrw, iwork, liw, rpar, ipar, jac)
c***begin prologue  debdf
c***purpose  solve an initial value problem in ordinary differential
c            equations using backward differentiation formulas.  it is
c            intended primarily for stiff problems.
c***library   slatec (depac)
c***category  i1a2
c***type      single precision (debdf-s, ddebdf-d)
c***keywords  backward differentiation formulas, depac,
c             initial value problems, ode,
c             ordinary differential equations, stiff
c***author  shampine, l. f., (snla)
c           watts, h. a., (snla)
c***description
c
c   this is the backward differentiation code in the package of
c   differential equation solvers depac, consisting of the codes
c   derkf, deabm, and debdf.  design of the package was by
c   l. f. shampine and h. a. watts.  it is documented in
c        sand-79-2374 , depac - design of a user oriented package of ode
c                              solvers.
c   debdf is a driver for a modification of the code lsode written by
c             a. c. hindmarsh
c             lawrence livermore laboratory
c             livermore, california 94550
c
c **********************************************************************
c **             depac package overview           **
c **********************************************************************
c
c        you have a choice of three differential equation solvers from
c        depac.  the following brief descriptions are meant to aid you
c        in choosing the most appropriate code for your problem.
c
c        derkf is a fifth order runge-kutta code.  it is the simplest of
c        the three choices, both algorithmically and in the use of the
c        code.  derkf is primarily designed to solve non-stiff and mild-
c        ly stiff differential equations when derivative evaluations are
c        not expensive.  it should generally not be used to get high
c        accuracy results nor answers at a great many specific points.
c        because derkf has very low overhead costs, it will usually
c        result in the least expensive integration when solving
c        problems requiring a modest amount of accuracy and having
c        equations that are not costly to evaluate.  derkf attempts to
c        discover when it is not suitable for the task posed.
c
c        deabm is a variable order (one through twelve) adams code.
c        its complexity lies somewhere between that of derkf and debdf.
c        deabm is primarily designed to solve non-stiff and mildly
c        stiff differential equations when derivative evaluations are
c        expensive, high accuracy results are needed or answers at
c        many specific points are required.  deabm attempts to discover
c        when it is not suitable for the task posed.
c
c        debdf is a variable order (one through five) backward
c        differentiation formula code.  it is the most complicated of
c        the three choices.  debdf is primarily designed to solve stiff
c        differential equations at crude to moderate tolerances.
c        if the problem is very stiff at all, derkf and deabm will be
c        quite inefficient compared to debdf.  however, debdf will be
c        inefficient compared to derkf and deabm on non-stiff problems
c        because it uses much more storage, has a much larger overhead,
c        and the low order formulas will not give high accuracies
c        efficiently.
c
c        the concept of stiffness cannot be described in a few words.
c        if you do not know the problem to be stiff, try either derkf
c        or deabm.  both of these codes will inform you of stiffness
c        when the cost of solving such problems becomes important.
c
c **********************************************************************
c ** abstract **
c **********************************************************************
c
c   subroutine debdf uses the backward differentiation formulas of
c   orders one through five to integrate a system of neq first order
c   ordinary differential equations of the form
c                         du/dx = f(x,u)
c   when the vector y(*) of initial values for u(*) at x=t is given. the
c   subroutine integrates from t to tout.  it is easy to continue the
c   integration to get results at additional tout.  this is the interval
c   mode of operation.  it is also easy for the routine to return with
c   the solution at each intermediate step on the way to tout.  this is
c   the intermediate-output mode of operation.
c
c **********************************************************************
c ** description of the arguments to debdf (an overview) **
c **********************************************************************
c
c   the parameters are:
c
c      f -- this is the name of a subroutine which you provide to
c             define the differential equations.
c
c      neq -- this is the number of (first order) differential
c             equations to be integrated.
c
c      t -- this is a value of the independent variable.
c
c      y(*) -- this array contains the solution components at t.
c
c      tout -- this is a point at which a solution is desired.
c
c      info(*) -- the basic task of the code is to integrate the
c             differential equations from t to tout and return an
c             answer at tout.  info(*) is an integer array which is used
c             to communicate exactly how you want this task to be
c             carried out.
c
c      rtol, atol -- these quantities
c             represent relative and absolute error tolerances which you
c             provide to indicate how accurately you wish the solution
c             to be computed.  you may choose them to be both scalars
c             or else both vectors.
c
c      idid -- this scalar quantity is an indicator reporting what
c             the code did.  you must monitor this integer variable to
c             decide what action to take next.
c
c      rwork(*), lrw -- rwork(*) is a real work array of
c             length lrw which provides the code with needed storage
c             space.
c
c      iwork(*), liw -- iwork(*) is an integer work array of length liw
c             which provides the code with needed storage space and an
c             across call flag.
c
c      rpar, ipar -- these are real and integer parameter
c             arrays which you can use for communication between your
c             calling program and the f subroutine (and the jac
c             subroutine).
c
c      jac -- this is the name of a subroutine which you may choose to
c             provide for defining the jacobian matrix of partial
c             derivatives df/du.
c
c  quantities which are used as input items are
c             neq, t, y(*), tout, info(*),
c             rtol, atol, rwork(1), lrw,
c             iwork(1), iwork(2), and liw.
c
c  quantities which may be altered by the code are
c             t, y(*), info(1), rtol, atol,
c             idid, rwork(*) and iwork(*).
c
c **********************************************************************
c * input -- what to do on the first call to debdf *
c **********************************************************************
c
c   the first call of the code is defined to be the start of each new
c   problem.  read through the descriptions of all the following items,
c   provide sufficient storage space for designated arrays, set
c   appropriate variables for the initialization of the problem, and
c   give information about how you want the problem to be solved.
c
c
c      f -- provide a subroutine of the form
c                               f(x,u,uprime,rpar,ipar)
c             to define the system of first order differential equations
c             which is to be solved. for the given values of x and the
c             vector  u(*)=(u(1),u(2),...,u(neq)) , the subroutine must
c             evaluate the neq components of the system of differential
c             equations  du/dx=f(x,u)  and store the derivatives in the
c             array uprime(*), that is,  uprime(i) = * du(i)/dx *  for
c             equations i=1,...,neq.
c
c             subroutine f must not alter x or u(*).  you must declare
c             the name f in an external statement in your program that
c             calls debdf.  you must dimension u and uprime in f.
c
c             rpar and ipar are real and integer parameter arrays which
c             you can use for communication between your calling program
c             and subroutine f.  they are not used or altered by debdf.
c             if you do not need rpar or ipar, ignore these parameters
c             by treating them as dummy arguments.  if you do choose to
c             use them, dimension them in your calling program and in f
c             as arrays of appropriate length.
c
c      neq -- set it to the number of differential equations.
c             (neq .ge. 1)
c
c      t -- set it to the initial point of the integration.
c             you must use a program variable for t because the code
c             changes its value.
c
c      y(*) -- set this vector to the initial values of the neq solution
c             components at the initial point.  you must dimension y at
c             least neq in your calling program.
c
c      tout -- set it to the first point at which a solution is desired.
c             you can take tout = t, in which case the code
c             will evaluate the derivative of the solution at t and
c             return.  integration either forward in t  (tout .gt. t)
c             or backward in t  (tout .lt. t)  is permitted.
c
c             the code advances the solution from t to tout using
c             step sizes which are automatically selected so as to
c             achieve the desired accuracy.  if you wish, the code will
c             return with the solution and its derivative following
c             each intermediate step (intermediate-output mode) so that
c             you can monitor them, but you still must provide tout in
c             accord with the basic aim of the code.
c
c             the first step taken by the code is a critical one
c             because it must reflect how fast the solution changes near
c             the initial point.  the code automatically selects an
c             initial step size which is practically always suitable for
c             the problem.  by using the fact that the code will not
c             step past tout in the first step, you could, if necessary,
c             restrict the length of the initial step size.
c
c             for some problems it may not be permissible to integrate
c             past a point tstop because a discontinuity occurs there
c             or the solution or its derivative is not defined beyond
c             tstop.  when you have declared a tstop point (see info(4)
c             and rwork(1)), you have told the code not to integrate
c             past tstop.  in this case any tout beyond tstop is invalid
c             input.
c
c      info(*) -- use the info array to give the code more details about
c             how you want your problem solved.  this array should be
c             dimensioned of length 15 to accommodate other members of
c             depac or possible future extensions, though debdf uses
c             only the first six entries.  you must respond to all of
c             the following items which are arranged as questions.  the
c             simplest use of the code corresponds to answering all
c             questions as yes ,i.e. setting all entries of info to 0.
c
c        info(1) -- this parameter enables the code to initialize
c               itself.  you must set it to indicate the start of every
c               new problem.
c
c            **** is this the first call for this problem ...
c                  yes -- set info(1) = 0
c                   no -- not applicable here.
c                         see below for continuation calls.  ****
c
c        info(2) -- how much accuracy you want of your solution
c               is specified by the error tolerances rtol and atol.
c               the simplest use is to take them both to be scalars.
c               to obtain more flexibility, they can both be vectors.
c               the code must be told your choice.
c
c            **** are both error tolerances rtol, atol scalars ...
c                  yes -- set info(2) = 0
c                         and input scalars for both rtol and atol
c                   no -- set info(2) = 1
c                         and input arrays for both rtol and atol ****
c
c        info(3) -- the code integrates from t in the direction
c               of tout by steps.  if you wish, it will return the
c               computed solution and derivative at the next
c               intermediate step (the intermediate-output mode) or
c               tout, whichever comes first.  this is a good way to
c               proceed if you want to see the behavior of the solution.
c               if you must have solutions at a great many specific
c               tout points, this code will compute them efficiently.
c
c            **** do you want the solution only at
c                 tout (and not at the next intermediate step) ...
c                  yes -- set info(3) = 0
c                   no -- set info(3) = 1 ****
c
c        info(4) -- to handle solutions at a great many specific
c               values tout efficiently, this code may integrate past
c               tout and interpolate to obtain the result at tout.
c               sometimes it is not possible to integrate beyond some
c               point tstop because the equation changes there or it is
c               not defined past tstop.  then you must tell the code
c               not to go past.
c
c            **** can the integration be carried out without any
c                 restrictions on the independent variable t ...
c                  yes -- set info(4)=0
c                   no -- set info(4)=1
c                         and define the stopping point tstop by
c                         setting rwork(1)=tstop ****
c
c        info(5) -- to solve stiff problems it is necessary to use the
c               jacobian matrix of partial derivatives of the system
c               of differential equations.  if you do not provide a
c               subroutine to evaluate it analytically (see the
c               description of the item jac in the call list), it will
c               be approximated by numerical differencing in this code.
c               although it is less trouble for you to have the code
c               compute partial derivatives by numerical differencing,
c               the solution will be more reliable if you provide the
c               derivatives via jac.  sometimes numerical differencing
c               is cheaper than evaluating derivatives in jac and
c               sometimes it is not - this depends on your problem.
c
c               if your problem is linear, i.e. has the form
c               du/dx = f(x,u) = j(x)*u + g(x)   for some matrix j(x)
c               and vector g(x), the jacobian matrix  df/du = j(x).
c               since you must provide a subroutine to evaluate f(x,u)
c               analytically, it is little extra trouble to provide
c               subroutine jac for evaluating j(x) analytically.
c               furthermore, in such cases, numerical differencing is
c               much more expensive than analytic evaluation.
c
c            **** do you want the code to evaluate the partial
c                 derivatives automatically by numerical differences ...
c                  yes -- set info(5)=0
c                   no -- set info(5)=1
c                         and provide subroutine jac for evaluating the
c                         jacobian matrix ****
c
c        info(6) -- debdf will perform much better if the jacobian
c               matrix is banded and the code is told this.  in this
c               case, the storage needed will be greatly reduced,
c               numerical differencing will be performed more cheaply,
c               and a number of important algorithms will execute much
c               faster.  the differential equation is said to have
c               half-bandwidths ml (lower) and mu (upper) if equation i
c               involves only unknowns y(j) with
c                              i-ml .le. j .le. i+mu
c               for all i=1,2,...,neq.  thus, ml and mu are the widths
c               of the lower and upper parts of the band, respectively,
c               with the main diagonal being excluded.  if you do not
c               indicate that the equation has a banded jacobian,
c               the code works with a full matrix of neq**2 elements
c               (stored in the conventional way).  computations with
c               banded matrices cost less time and storage than with
c               full matrices if  2*ml+mu .lt. neq.  if you tell the
c               code that the jacobian matrix has a banded structure and
c               you want to provide subroutine jac to compute the
c               partial derivatives, then you must be careful to store
c               the elements of the jacobian matrix in the special form
c               indicated in the description of jac.
c
c            **** do you want to solve the problem using a full
c                 (dense) jacobian matrix (and not a special banded
c                 structure) ...
c                  yes -- set info(6)=0
c                   no -- set info(6)=1
c                         and provide the lower (ml) and upper (mu)
c                         bandwidths by setting
c                         iwork(1)=ml
c                         iwork(2)=mu ****
c
c      rtol, atol -- you must assign relative (rtol) and absolute (atol)
c             error tolerances to tell the code how accurately you want
c             the solution to be computed.  they must be defined as
c             program variables because the code may change them.  you
c             have two choices --
c                  both rtol and atol are scalars. (info(2)=0)
c                  both rtol and atol are vectors. (info(2)=1)
c             in either case all components must be non-negative.
c
c             the tolerances are used by the code in a local error test
c             at each step which requires roughly that
c                     abs(local error) .le. rtol*abs(y)+atol
c             for each vector component.
c             (more specifically, a root-mean-square norm is used to
c             measure the size of vectors, and the error test uses the
c             magnitude of the solution at the beginning of the step.)
c
c             the true (global) error is the difference between the true
c             solution of the initial value problem and the computed
c             approximation.  practically all present day codes,
c             including this one, control the local error at each step
c             and do not even attempt to control the global error
c             directly.  roughly speaking, they produce a solution y(t)
c             which satisfies the differential equations with a
c             residual r(t),    dy(t)/dt = f(t,y(t)) + r(t)   ,
c             and, almost always, r(t) is bounded by the error
c             tolerances.  usually, but not always, the true accuracy of
c             the computed y is comparable to the error tolerances. this
c             code will usually, but not always, deliver a more accurate
c             solution if you reduce the tolerances and integrate again.
c             by comparing two such solutions you can get a fairly
c             reliable idea of the true error in the solution at the
c             bigger tolerances.
c
c             setting atol=0. results in a pure relative error test on
c             that component.  setting rtol=0. results in a pure abso-
c             lute error test on that component.  a mixed test with non-
c             zero rtol and atol corresponds roughly to a relative error
c             test when the solution component is much bigger than atol
c             and to an absolute error test when the solution component
c             is smaller than the threshold atol.
c
c             proper selection of the absolute error control parameters
c             atol  requires you to have some idea of the scale of the
c             solution components.  to acquire this information may mean
c             that you will have to solve the problem more than once. in
c             the absence of scale information, you should ask for some
c             relative accuracy in all the components (by setting  rtol
c             values non-zero) and perhaps impose extremely small
c             absolute error tolerances to protect against the danger of
c             a solution component becoming zero.
c
c             the code will not attempt to compute a solution at an
c             accuracy unreasonable for the machine being used.  it will
c             advise you if you ask for too much accuracy and inform
c             you as to the maximum accuracy it believes possible.
c
c      rwork(*) -- dimension this real work array of length lrw in your
c             calling program.
c
c      rwork(1) -- if you have set info(4)=0, you can ignore this
c             optional input parameter.  otherwise you must define a
c             stopping point tstop by setting   rwork(1) = tstop.
c             (for some problems it may not be permissible to integrate
c             past a point tstop because a discontinuity occurs there
c             or the solution or its derivative is not defined beyond
c             tstop.)
c
c      lrw -- set it to the declared length of the rwork array.
c             you must have
c                  lrw .ge. 250+10*neq+neq**2
c             for the full (dense) jacobian case (when info(6)=0),  or
c                  lrw .ge. 250+10*neq+(2*ml+mu+1)*neq
c             for the banded jacobian case (when info(6)=1).
c
c      iwork(*) -- dimension this integer work array of length liw in
c             your calling program.
c
c      iwork(1), iwork(2) -- if you have set info(6)=0, you can ignore
c             these optional input parameters. otherwise you must define
c             the half-bandwidths ml (lower) and mu (upper) of the
c             jacobian matrix by setting    iwork(1) = ml   and
c             iwork(2) = mu.  (the code will work with a full matrix
c             of neq**2 elements unless it is told that the problem has
c             a banded jacobian, in which case the code will work with
c             a matrix containing at most  (2*ml+mu+1)*neq  elements.)
c
c      liw -- set it to the declared length of the iwork array.
c             you must have liw .ge. 56+neq.
c
c      rpar, ipar -- these are parameter arrays, of real and integer
c             type, respectively.  you can use them for communication
c             between your program that calls debdf and the  f
c             subroutine (and the jac subroutine).  they are not used or
c             altered by debdf.  if you do not need rpar or ipar, ignore
c             these parameters by treating them as dummy arguments.  if
c             you do choose to use them, dimension them in your calling
c             program and in f (and in jac) as arrays of appropriate
c             length.
c
c      jac -- if you have set info(5)=0, you can ignore this parameter
c             by treating it as a dummy argument. (for some compilers
c             you may have to write a dummy subroutine named  jac  in
c             order to avoid problems associated with missing external
c             routine names.)  otherwise, you must provide a subroutine
c             of the form
c                          jac(x,u,pd,nrowpd,rpar,ipar)
c             to define the jacobian matrix of partial derivatives df/du
c             of the system of differential equations   du/dx = f(x,u).
c             for the given values of x and the vector
c             u(*)=(u(1),u(2),...,u(neq)), the subroutine must evaluate
c             the non-zero partial derivatives  df(i)/du(j)  for each
c             differential equation i=1,...,neq and each solution
c             component j=1,...,neq , and store these values in the
c             matrix pd.  the elements of pd are set to zero before each
c             call to jac so only non-zero elements need to be defined.
c
c             subroutine jac must not alter x, u(*), or nrowpd.  you
c             must declare the name jac in an external statement in your
c             program that calls debdf.  nrowpd is the row dimension of
c             the pd matrix and is assigned by the code.  therefore you
c             must dimension pd in jac according to
c                              dimension pd(nrowpd,1)
c             you must also dimension u in jac.
c
c             the way you must store the elements into the pd matrix
c             depends on the structure of the jacobian which you
c             indicated by info(6).
c             *** info(6)=0 -- full (dense) jacobian ***
c                 when you evaluate the (non-zero) partial derivative
c                 of equation i with respect to variable j, you must
c                 store it in pd according to
c                                pd(i,j) = * df(i)/du(j) *
c             *** info(6)=1 -- banded jacobian with ml lower and mu
c                 upper diagonal bands (refer to info(6) description of
c                 ml and mu) ***
c                 when you evaluate the (non-zero) partial derivative
c                 of equation i with respect to variable j, you must
c                 store it in pd according to
c                                irow = i - j + ml + mu + 1
c                                pd(irow,j) = * df(i)/du(j) *
c
c             rpar and ipar are real and integer parameter
c             arrays which you can use for communication between your
c             calling program and your jacobian subroutine jac.  they
c             are not altered by debdf.  if you do not need rpar or
c             ipar, ignore these parameters by treating them as dummy
c             arguments.  if you do choose to use them, dimension them
c             in your calling program and in jac as arrays of
c             appropriate length.
c
c **********************************************************************
c * output -- after any return from ddebdf *
c **********************************************************************
c
c   the principal aim of the code is to return a computed solution at
c   tout, although it is also possible to obtain intermediate results
c   along the way.  to find out whether the code achieved its goal
c   or if the integration process was interrupted before the task was
c   completed, you must check the idid parameter.
c
c
c      t -- the solution was successfully advanced to the
c             output value of t.
c
c      y(*) -- contains the computed solution approximation at t.
c             you may also be interested in the approximate derivative
c             of the solution at t.  it is contained in
c             rwork(21),...,rwork(20+neq).
c
c      idid -- reports what the code did
c
c                         *** task completed ***
c                   reported by positive values of idid
c
c             idid = 1 -- a step was successfully taken in the
c                       intermediate-output mode.  the code has not
c                       yet reached tout.
c
c             idid = 2 -- the integration to tout was successfully
c                       completed (t=tout) by stepping exactly to tout.
c
c             idid = 3 -- the integration to tout was successfully
c                       completed (t=tout) by stepping past tout.
c                       y(*) is obtained by interpolation.
c
c                         *** task interrupted ***
c                   reported by negative values of idid
c
c             idid = -1 -- a large amount of work has been expended.
c                       (500 steps attempted)
c
c             idid = -2 -- the error tolerances are too stringent.
c
c             idid = -3 -- the local error test cannot be satisfied
c                       because you specified a zero component in atol
c                       and the corresponding computed solution
c                       component is zero.  thus, a pure relative error
c                       test is impossible for this component.
c
c             idid = -4,-5  -- not applicable for this code but used
c                       by other members of depac.
c
c             idid = -6 -- debdf had repeated convergence test failures
c                       on the last attempted step.
c
c             idid = -7 -- debdf had repeated error test failures on
c                       the last attempted step.
c
c             idid = -8,..,-32  -- not applicable for this code but
c                       used by other members of depac or possible
c                       future extensions.
c
c                         *** task terminated ***
c                   reported by the value of idid=-33
c
c             idid = -33 -- the code has encountered trouble from which
c                       it cannot recover.  a message is printed
c                       explaining the trouble and control is returned
c                       to the calling program.  for example, this
c                       occurs when invalid input is detected.
c
c      rtol, atol -- these quantities remain unchanged except when
c             idid = -2.  in this case, the error tolerances have been
c             increased by the code to values which are estimated to be
c             appropriate for continuing the integration.  however, the
c             reported solution at t was obtained using the input values
c             of rtol and atol.
c
c      rwork, iwork -- contain information which is usually of no
c             interest to the user but necessary for subsequent calls.
c             however, you may find use for
c
c             rwork(11)--which contains the step size h to be
c                        attempted on the next step.
c
c             rwork(12)--if the tolerances have been increased by the
c                        code (idid = -2) , they were multiplied by the
c                        value in rwork(12).
c
c             rwork(13)--which contains the current value of the
c                        independent variable, i.e. the farthest point
c                        integration has reached.  this will be
c                        different from t only when interpolation has
c                        been performed (idid=3).
c
c             rwork(20+i)--which contains the approximate derivative
c                        of the solution component y(i).  in debdf, it
c                        is never obtained by calling subroutine f to
c                        evaluate the differential equation using t and
c                        y(*), except at the initial point of
c                        integration.
c
c **********************************************************************
c ** input -- what to do to continue the integration **
c **             (calls after the first)             **
c **********************************************************************
c
c        this code is organized so that subsequent calls to continue the
c        integration involve little (if any) additional effort on your
c        part. you must monitor the idid parameter in order to determine
c        what to do next.
c
c        recalling that the principal task of the code is to integrate
c        from t to tout (the interval mode), usually all you will need
c        to do is specify a new tout upon reaching the current tout.
c
c        do not alter any quantity not specifically permitted below,
c        in particular do not alter neq, t, y(*), rwork(*), iwork(*) or
c        the differential equation in subroutine f.  any such alteration
c        constitutes a new problem and must be treated as such, i.e.
c        you must start afresh.
c
c        you cannot change from vector to scalar error control or vice
c        versa (info(2)) but you can change the size of the entries of
c        rtol, atol.  increasing a tolerance makes the equation easier
c        to integrate.  decreasing a tolerance will make the equation
c        harder to integrate and should generally be avoided.
c
c        you can switch from the intermediate-output mode to the
c        interval mode (info(3)) or vice versa at any time.
c
c        if it has been necessary to prevent the integration from going
c        past a point tstop (info(4), rwork(1)), keep in mind that the
c        code will not integrate to any tout beyond the currently
c        specified tstop.  once tstop has been reached you must change
c        the value of tstop or set info(4)=0.  you may change info(4)
c        or tstop at any time but you must supply the value of tstop in
c        rwork(1) whenever you set info(4)=1.
c
c        do not change info(5), info(6), iwork(1), or iwork(2)
c        unless you are going to restart the code.
c
c        the parameter info(1) is used by the code to indicate the
c        beginning of a new problem and to indicate whether integration
c        is to be continued.  you must input the value  info(1) = 0
c        when starting a new problem.  you must input the value
c        info(1) = 1  if you wish to continue after an interrupted task.
c        do not set  info(1) = 0  on a continuation call unless you
c        want the code to restart at the current t.
c
c                         *** following a completed task ***
c         if
c             idid = 1, call the code again to continue the integration
c                     another step in the direction of tout.
c
c             idid = 2 or 3, define a new tout and call the code again.
c                     tout must be different from t.  you cannot change
c                     the direction of integration without restarting.
c
c                         *** following an interrupted task ***
c                     to show the code that you realize the task was
c                     interrupted and that you want to continue, you
c                     must take appropriate action and reset info(1) = 1
c         if
c             idid = -1, the code has attempted 500 steps.
c                     if you want to continue, set info(1) = 1 and
c                     call the code again.  an additional 500 steps
c                     will be allowed.
c
c             idid = -2, the error tolerances rtol, atol have been
c                     increased to values the code estimates appropriate
c                     for continuing.  you may want to change them
c                     yourself.  if you are sure you want to continue
c                     with relaxed error tolerances, set info(1)=1 and
c                     call the code again.
c
c             idid = -3, a solution component is zero and you set the
c                     corresponding component of atol to zero.  if you
c                     are sure you want to continue, you must first
c                     alter the error criterion to use positive values
c                     for those components of atol corresponding to zero
c                     solution components, then set info(1)=1 and call
c                     the code again.
c
c             idid = -4,-5  --- cannot occur with this code but used
c                     by other members of depac.
c
c             idid = -6, repeated convergence test failures occurred
c                     on the last attempted step in debdf.  an inaccu-
c                     rate jacobian may be the problem.  if you are
c                     absolutely certain you want to continue, restart
c                     the integration at the current t by setting
c                     info(1)=0 and call the code again.
c
c             idid = -7, repeated error test failures occurred on the
c                     last attempted step in debdf.  a singularity in
c                     the solution may be present.  you should re-
c                     examine the problem being solved.  if you are
c                     absolutely certain you want to continue, restart
c                     the integration at the current t by setting
c                     info(1)=0 and call the code again.
c
c             idid = -8,..,-32  --- cannot occur with this code but
c                     used by other members of depac or possible future
c                     extensions.
c
c                         *** following a terminated task ***
c         if
c             idid = -33, you cannot continue the solution of this
c                     problem.  an attempt to do so will result in your
c                     run being terminated.
c
c **********************************************************************
c
c         ***** warning *****
c
c     if debdf is to be used in an overlay situation, you must save and
c     restore certain items used internally by debdf  (values in the
c     common block debdf1).  this can be accomplished as follows.
c
c     to save the necessary values upon return from debdf, simply call
c        svco(rwork(22+neq),iwork(21+neq)).
c
c     to restore the necessary values before the next call to debdf,
c     simply call    rsco(rwork(22+neq),iwork(21+neq)).
c
c***references  l. f. shampine and h. a. watts, depac - design of a user
c                 oriented package of ode solvers, report sand79-2374,
c                 sandia laboratories, 1979.
c***routines called  lsod, xermsg
c***common blocks    debdf1
c***revision history  (yymmdd)
c   800901  date written
c   890831  modified array declarations.  (wrb)
c   891024  changed references from vnorm to hvnrm.  (wrb)
c   891024  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900326  removed duplicate information from description section.
c           (wrb)
c   900510  convert xerrwv calls to xermsg calls, change prologue
c           comments to agree with ddebdf.  (rwc)
c   920501  reformatted the references section.  (wrb)
c***end prologue  debdf
c
c
      logical intout
      character*8 xern1, xern2
      character*16 xern3
c
      dimension y(*),info(15),rtol(*),atol(*),rwork(*),iwork(*),
     1          rpar(*),ipar(*)
c
      common /debdf1/ told, rowns(210),
     1   el0, h, hmin, hmxi, hu, tn, uround,
     2   iquit, init, iyh, iewt, iacor, isavf, iwm, ksteps,
     3   ibegin, itol, iinteg, itstop, ijac, iband, iowns(6),
     4   ier, jstart, kflag, l, meth, miter, maxord, n, nq, nst, nfe,
     5   nje, nqu
c
      external f, jac
c
c        check for an apparent infinite loop
c
c***first executable statement  debdf
      if (info(1) .eq. 0) iwork(liw) = 0
c
      if (iwork(liw).ge. 5) then
         if (t .eq. rwork(21+neq)) then
            write (xern3, '(1pe15.6)') t
            call xermsg ('slatec', 'debdf',
     *         'an apparent infinite loop has been detected.$$' //
     *         'you have made repeated calls at t = ' // xern3 //
     *         ' and the integration has not advanced.  check the ' //
     *         'way you have set parameters for the call to the ' //
     *         'code particularly info(1).', 13, 2)
            return
         endif
      endif
c
      idid = 0
c
c        check validity of info parameters
c
      if (info(1) .ne. 0 .and. info(1) .ne. 1) then
         write (xern1, '(i8)') info(1)
         call xermsg ('slatec', 'debdf', 'info(1) must be set to 0 ' //
     *      'for the  start of a new problem, and must be set to 1 ' //
     *      'following an interrupted task.  you are attempting to ' //
     *      'continue the integration illegally by calling the ' //
     *      'code with  info(1) = ' // xern1, 3, 1)
         idid = -33
      endif
c
      if (info(2) .ne. 0 .and. info(2) .ne. 1) then
         write (xern1, '(i8)') info(2)
         call xermsg ('slatec', 'debdf', 'info(2) must be 0 or 1 ' //
     *      'indicating scalar and vector error tolerances, ' //
     *      'respectively.  you have called the code with info(2) = ' //
     *      xern1, 4, 1)
         idid = -33
      endif
c
      if (info(3) .ne. 0 .and. info(3) .ne. 1) then
         write (xern1, '(i8)') info(3)
         call xermsg ('slatec', 'debdf', 'info(3) must be 0 or 1 ' //
     *      'indicating the interval or intermediate-output mode of ' //
     *      'integration, respectively.  you have called the code ' //
     *      'with  info(3) = ' // xern1, 5, 1)
         idid = -33
      endif
c
      if (info(4) .ne. 0 .and. info(4) .ne. 1) then
         write (xern1, '(i8)') info(4)
         call xermsg ('slatec', 'debdf', 'info(4) must be 0 or 1 ' //
     *      'indicating whether or not the integration interval is ' //
     *      'to be restricted by a point tstop.  you have called ' //
     *      'the code  with info(4) = ' // xern1, 14, 1)
         idid = -33
      endif
c
      if (info(5) .ne. 0 .and. info(5) .ne. 1) then
         write (xern1, '(i8)') info(5)
         call xermsg ('slatec',  'debdf', 'info(5) must be 0 or 1 ' //
     *      'indicating whether the code is told to form the ' //
     *      'jacobian matrix by numerical differencing or you ' //
     *      'provide a subroutine to evaluate it analytically.  ' //
     *      'you have called the code with info(5) = ' // xern1, 15, 1)
         idid = -33
      endif
c
      if (info(6) .ne. 0 .and. info(6) .ne. 1) then
         write (xern1, '(i8)') info(6)
         call xermsg ('slatec', 'debdf', 'info(6) must be 0 or 1 ' //
     *      'indicating whether the code is told to treat the ' //
     *      'jacobian as a full (dense) matrix or as having a ' //
     *      'special banded structure.  you have called the code ' //
     *      'with info(6) = ' // xern1, 16, 1)
         idid = -33
      endif
c
      ilrw = neq
      if (info(6) .ne. 0) then
c
c        check bandwidth parameters
c
         ml = iwork(1)
         mu = iwork(2)
         ilrw = 2*ml + mu + 1
c
         if (ml.lt.0 .or. ml.ge.neq .or. mu.lt.0 .or. mu.ge.neq) then
            write (xern1, '(i8)') ml
            write (xern2, '(i8)') mu
            call xermsg ('slatec', 'debdf', 'you have set info(6) ' //
     *         '= 1, telling the code that the jacobian matrix has ' //
     *         'a special banded structure.  however, the lower ' //
     *         '(upper) bandwidths  ml (mu) violate the constraints ' //
     *         'ml,mu .ge. 0 and  ml,mu .lt. neq.  you have called ' //
     *         'the code with ml = ' // xern1 // ' and mu = ' // xern2,
     *         17, 1)
            idid = -33
         endif
      endif
c
c        check lrw and liw for sufficient storage allocation
c
      if (lrw .lt. 250 + (10 + ilrw)*neq) then
         write (xern1, '(i8)') lrw
         if (info(6) .eq. 0) then
            call xermsg ('slatec', 'debdf', 'length of array rwork ' //
     *         'must be at least 250 + 10*neq + neq*neq.$$' //
     *         'you have called the code with  lrw = ' // xern1, 1, 1)
         else
            call xermsg ('slatec', 'debdf', 'length of array rwork ' //
     *         'must be at least 250 + 10*neq + (2*ml+mu+1)*neq.$$' //
     *         'you have called the code with  lrw = ' // xern1, 18, 1)
         endif
         idid = -33
      endif
c
      if (liw .lt. 56 + neq) then
         write (xern1, '(i8)') liw
         call xermsg ('slatec', 'debdf', 'length of array iwork ' //
     *      'be at least  56 + neq.  you have called the code with ' //
     *      'liw = ' // xern1, 2, 1)
         idid = -33
      endif
c
c        compute the indices for the arrays to be stored in the work
c        array and restore common block data
c
      icomi = 21 + neq
      iinout = icomi + 33
c
      iypout = 21
      itstar = 21 + neq
      icomr = 22 + neq
c
      if (info(1) .ne. 0) intout = iwork(iinout) .ne. (-1)
c     call rsco(rwork(icomr),iwork(icomi))
c
      iyh = icomr + 218
      iewt = iyh + 6*neq
      isavf = iewt + neq
      iacor = isavf + neq
      iwm = iacor + neq
      idelsn = iwm + 2 + ilrw*neq
c
      ibegin = info(1)
      itol = info(2)
      iinteg = info(3)
      itstop = info(4)
      ijac = info(5)
      iband = info(6)
      rwork(itstar) = t
c
      call lsod(f,neq,t,y,tout,rtol,atol,idid,rwork(iypout),
     1          rwork(iyh),rwork(iyh),rwork(iewt),rwork(isavf),
     2          rwork(iacor),rwork(iwm),iwork(1),jac,intout,
     3          rwork(1),rwork(12),rwork(idelsn),rpar,ipar)
c
      iwork(iinout) = -1
      if (intout) iwork(iinout) = 1
c
      if (idid .ne. (-2)) iwork(liw) = iwork(liw) + 1
      if (t .ne. rwork(itstar)) iwork(liw) = 0
c     call svco(rwork(icomr),iwork(icomi))
      rwork(11) = h
      rwork(13) = tn
      info(1) = ibegin
c
      return
      end
