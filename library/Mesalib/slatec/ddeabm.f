*deck ddeabm
      subroutine ddeabm (df, neq, t, y, tout, info, rtol, atol, idid,
     +   rwork, lrw, iwork, liw, rpar, ipar)
c***begin prologue  ddeabm
c***purpose  solve an initial value problem in ordinary differential
c            equations using an adams-bashforth method.
c***library   slatec (depac)
c***category  i1a1b
c***type      double precision (deabm-s, ddeabm-d)
c***keywords  adams-bashforth method, depac, initial value problems,
c             ode, ordinary differential equations, predictor-corrector
c***author  shampine, l. f., (snla)
c           watts, h. a., (snla)
c***description
c
c   this is the adams code in the package of differential equation
c   solvers depac, consisting of the codes dderkf, ddeabm, and ddebdf.
c   design of the package was by l. f. shampine and h. a. watts.
c   it is documented in
c        sand79-2374 , depac - design of a user oriented package of ode
c                              solvers.
c   ddeabm is a driver for a modification of the code ode written by
c             l. f. shampine and m. k. gordon
c             sandia laboratories
c             albuquerque, new mexico 87185
c
c **********************************************************************
c * abstract *
c ************
c
c   subroutine ddeabm uses the adams-bashforth-moulton
c   predictor-corrector formulas of orders one through twelve to
c   integrate a system of neq first order ordinary differential
c   equations of the form
c                         du/dx = df(x,u)
c   when the vector y(*) of initial values for u(*) at x=t is given.
c   the subroutine integrates from t to tout. it is easy to continue the
c   integration to get results at additional tout.  this is the interval
c   mode of operation.  it is also easy for the routine to return with
c   the solution at each intermediate step on the way to tout.  this is
c   the intermediate-output mode of operation.
c
c   ddeabm uses subprograms ddes, dsteps, dintp, dhstrt, dhvnrm,
c   d1mach, and the error handling routine xermsg.  the only machine
c   dependent parameters to be assigned appear in d1mach.
c
c **********************************************************************
c * description of the arguments to ddeabm (an overview) *
c **********************************************************************
c
c   the parameters are
c
c      df -- this is the name of a subroutine which you provide to
c             define the differential equations.
c
c      neq -- this is the number of (first order) differential
c             equations to be integrated.
c
c      t -- this is a double precision value of the independent
c           variable.
c
c      y(*) -- this double precision array contains the solution
c              components at t.
c
c      tout -- this is a double precision point at which a solution is
c              desired.
c
c      info(*) -- the basic task of the code is to integrate the
c             differential equations from t to tout and return an
c             answer at tout.  info(*) is an integer array which is used
c             to communicate exactly how you want this task to be
c             carried out.
c
c      rtol, atol -- these double precision quantities represent
c                    relative and absolute error tolerances which you
c                    provide to indicate how accurately you wish the
c                    solution to be computed.  you may choose them to be
c                    both scalars or else both vectors.
c
c      idid -- this scalar quantity is an indicator reporting what
c             the code did.  you must monitor this integer variable to
c             decide what action to take next.
c
c      rwork(*), lrw -- rwork(*) is a double precision work array of
c             length lrw which provides the code with needed storage
c             space.
c
c      iwork(*), liw -- iwork(*) is an integer work array of length liw
c             which provides the code with needed storage space and an
c             across call flag.
c
c      rpar, ipar -- these are double precision and integer parameter
c             arrays which you can use for communication between your
c             calling program and the df subroutine.
c
c  quantities which are used as input items are
c             neq, t, y(*), tout, info(*),
c             rtol, atol, rwork(1), lrw and liw.
c
c  quantities which may be altered by the code are
c             t, y(*), info(1), rtol, atol,
c             idid, rwork(*) and iwork(*).
c
c **********************************************************************
c * input -- what to do on the first call to ddeabm *
c **********************************************************************
c
c   the first call of the code is defined to be the start of each new
c   problem.  read through the descriptions of all the following items,
c   provide sufficient storage space for designated arrays, set
c   appropriate variables for the initialization of the problem, and
c   give information about how you want the problem to be solved.
c
c
c      df -- provide a subroutine of the form
c                               df(x,u,uprime,rpar,ipar)
c             to define the system of first order differential equations
c             which is to be solved.  for the given values of x and the
c             vector  u(*)=(u(1),u(2),...,u(neq)) , the subroutine must
c             evaluate the neq components of the system of differential
c             equations  du/dx=df(x,u)  and store the derivatives in the
c             array uprime(*), that is,  uprime(i) = * du(i)/dx *  for
c             equations i=1,...,neq.
c
c             subroutine df must not alter x or u(*).  you must declare
c             the name df in an external statement in your program that
c             calls ddeabm.  you must dimension u and uprime in df.
c
c             rpar and ipar are double precision and integer parameter
c             arrays which you can use for communication between your
c             calling program and subroutine df. they are not used or
c             altered by ddeabm.  if you do not need rpar or ipar,
c             ignore these parameters by treating them as dummy
c             arguments. if you do choose to use them, dimension them in
c             your calling program and in df as arrays of appropriate
c             length.
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
c      tout -- set it to the first point at which a solution
c             is desired.  you can take tout = t, in which case the code
c             will evaluate the derivative of the solution at t and
c             return. integration either forward in t  (tout .gt. t)  or
c             backward in t  (tout .lt. t)  is permitted.
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
c             the problem. by using the fact that the code will not step
c             past tout in the first step, you could, if necessary,
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
c             depac or possible future extensions, though ddeabm uses
c             only the first four entries.  you must respond to all of
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
c             (more specifically, a euclidean norm is used to measure
c             the size of vectors, and the error test uses the magnitude
c             of the solution at the beginning of the step.)
c
c             the true (global) error is the difference between the true
c             solution of the initial value problem and the computed
c             approximation.  practically all present day codes,
c             including this one, control the local error at each step
c             and do not even attempt to control the global error
c             directly.  roughly speaking, they produce a solution y(t)
c             which satisfies the differential equations with a
c             residual r(t),    dy(t)/dt = df(t,y(t)) + r(t)   ,
c             and, almost always, r(t) is bounded by the error
c             tolerances.  usually, but not always, the true accuracy of
c             the computed y is comparable to the error tolerances. this
c             code will usually, but not always, deliver a more accurate
c             solution if you reduce the tolerances and integrate again.
c             by comparing two such solutions you can get a fairly
c             reliable idea of the true error in the solution at the
c             bigger tolerances.
c
c             setting atol=0.d0 results in a pure relative error test on
c             that component. setting rtol=0. results in a pure absolute
c             error test on that component.  a mixed test with non-zero
c             rtol and atol corresponds roughly to a relative error
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
c      rwork(*) -- dimension this double precision work array of length
c             lrw in your calling program.
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
c             you must have  lrw .ge. 130+21*neq
c
c      iwork(*) -- dimension this integer work array of length liw in
c             your calling program.
c
c      liw -- set it to the declared length of the iwork array.
c             you must have  liw .ge. 51
c
c      rpar, ipar -- these are parameter arrays, of double precision and
c             integer type, respectively.  you can use them for
c             communication between your program that calls ddeabm and
c             the  df subroutine.  they are not used or altered by
c             ddeabm.  if you do not need rpar or ipar, ignore these
c             parameters by treating them as dummy arguments.  if you do
c             choose to use them, dimension them in your calling program
c             and in df as arrays of appropriate length.
c
c **********************************************************************
c * output -- after any return from ddeabm *
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
c             idid = -4 -- the problem appears to be stiff.
c
c             idid = -5,-6,-7,..,-32  -- not applicable for this code
c                       but used by other members of depac or possible
c                       future extensions.
c
c                         *** task terminated ***
c                   reported by the value of idid=-33
c
c             idid = -33 -- the code has encountered trouble from which
c                       it cannot recover.  a message is printed
c                       explaining the trouble and control is returned
c                       to the calling program. for example, this occurs
c                       when invalid input is detected.
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
c                        integration has reached. this will be different
c                        from t only when interpolation has been
c                        performed (idid=3).
c
c             rwork(20+i)--which contains the approximate derivative
c                        of the solution component y(i).  in ddeabm, it
c                        is obtained by calling subroutine df to
c                        evaluate the differential equation using t and
c                        y(*) when idid=1 or 2, and by interpolation
c                        when idid=3.
c
c **********************************************************************
c * input -- what to do to continue the integration *
c *             (calls after the first)             *
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
c        the differential equation in subroutine df. any such alteration
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
c                     tout must be different from t. you cannot change
c                     the direction of integration without restarting.
c
c                         *** following an interrupted task ***
c                     to show the code that you realize the task was
c                     interrupted and that you want to continue, you
c                     must take appropriate action and reset info(1) = 1
c         if
c             idid = -1, the code has attempted 500 steps.
c                     if you want to continue, set info(1) = 1 and
c                     call the code again. an additional 500 steps
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
c             idid = -4, the problem appears to be stiff.  it is very
c                     inefficient to solve such problems with ddeabm.
c                     the code ddebdf in depac handles this task
c                     efficiently.  if you are absolutely sure you want
c                     to continue with ddeabm, set info(1)=1 and call
c                     the code again.
c
c             idid = -5,-6,-7,..,-32  --- cannot occur with this code
c                     but used by other members of depac or possible
c                     future extensions.
c
c                         *** following a terminated task ***
c         if
c             idid = -33, you cannot continue the solution of this
c                     problem.  an attempt to do so will result in your
c                     run being terminated.
c
c **********************************************************************
c *long description:
c
c **********************************************************************
c *             depac package overview           *
c **********************************************************************
c
c ....   you have a choice of three differential equation solvers from
c ....   depac. the following brief descriptions are meant to aid you in
c ....   choosing the most appropriate code for your problem.
c
c ....   dderkf is a fifth order runge-kutta code. it is the simplest of
c ....   the three choices, both algorithmically and in the use of the
c ....   code. dderkf is primarily designed to solve non-stiff and
c ....   mildly stiff differential equations when derivative evaluations
c ....   are not expensive. it should generally not be used to get high
c ....   accuracy results nor answers at a great many specific points.
c ....   because dderkf has very low overhead costs, it will usually
c ....   result in the least expensive integration when solving
c ....   problems requiring a modest amount of accuracy and having
c ....   equations that are not costly to evaluate. dderkf attempts to
c ....   discover when it is not suitable for the task posed.
c
c ....   ddeabm is a variable order (one through twelve) adams code.
c ....   its complexity lies somewhere between that of dderkf and
c ....   ddebdf.  ddeabm is primarily designed to solve non-stiff and
c ....   mildly stiff differential equations when derivative evaluations
c ....   are expensive, high accuracy results are needed or answers at
c ....   many specific points are required. ddeabm attempts to discover
c ....   when it is not suitable for the task posed.
c
c ....   ddebdf is a variable order (one through five) backward
c ....   differentiation formula code. it is the most complicated of
c ....   the three choices. ddebdf is primarily designed to solve stiff
c ....   differential equations at crude to moderate tolerances.
c ....   if the problem is very stiff at all, dderkf and ddeabm will be
c ....   quite inefficient compared to ddebdf. however, ddebdf will be
c ....   inefficient compared to dderkf and ddeabm on non-stiff problems
c ....   because it uses much more storage, has a much larger overhead,
c ....   and the low order formulas will not give high accuracies
c ....   efficiently.
c
c ....   the concept of stiffness cannot be described in a few words.
c ....   if you do not know the problem to be stiff, try either dderkf
c ....   or ddeabm. both of these codes will inform you of stiffness
c ....   when the cost of solving such problems becomes important.
c
c *********************************************************************
c
c***references  l. f. shampine and h. a. watts, depac - design of a user
c                 oriented package of ode solvers, report sand79-2374,
c                 sandia laboratories, 1979.
c***routines called  ddes, xermsg
c***revision history  (yymmdd)
c   820301  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   891006  cosmetic changes to prologue.  (wrb)
c   891024  changed references from dvnorm to dhvnrm.  (wrb)
c   891024  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900510  convert xerrwv calls to xermsg calls.  (rwc)
c   920501  reformatted the references section.  (wrb)
c***end prologue  ddeabm
c
      integer ialpha, ibeta, idelsn, idid, ifouru, ig, ihold,
     1      info, ip, ipar, iphi, ipsi, isig, itold, itstar, itwou,
     2      iv, iw, iwork, iwt, iyp, iypout, iyy, liw, lrw, neq
      double precision atol, rpar, rtol, rwork, t, tout, y
      logical start,phase1,nornd,stiff,intout
c
      dimension y(*),info(15),rtol(*),atol(*),rwork(*),iwork(*),
     1          rpar(*),ipar(*)
c
      character*8 xern1
      character*16 xern3
c
      external df
c
c     check for an apparent infinite loop
c
c***first executable statement  ddeabm
      if ( info(1) .eq. 0 ) iwork(liw) = 0
      if (iwork(liw) .ge. 5) then
         if (t .eq. rwork(21 + neq)) then
            write (xern3, '(1pe15.6)') t
            call xermsg ('slatec', 'ddeabm',
     *         'an apparent infinite loop has been detected.$$' //
     *         'you have made repeated calls at t = ' // xern3 //
     *         ' and the integration has not advanced.  check the ' //
     *         'way you have set parameters for the call to the ' //
     *         'code, particularly info(1).', 13, 2)
            return
         endif
      endif
c
c     check lrw and liw for sufficient storage allocation
c
      idid=0
      if (lrw .lt. 130+21*neq) then
         write (xern1, '(i8)') lrw
         call xermsg ('slatec', 'ddeabm', 'the length of the rwork ' //
     *      'array must be at least 130 + 21*neq.$$' //
     *      'you have called the code with lrw = ' // xern1, 1, 1)
         idid=-33
      endif
c
      if (liw .lt. 51) then
         write (xern1, '(i8)') liw
         call xermsg ('slatec', 'ddeabm', 'the length of the iwork ' //
     *      'array must be at least 51.$$you have called the code ' //
     *      'with liw = ' // xern1, 2, 1)
         idid=-33
      endif
c
c     compute the indices for the arrays to be stored in the work array
c
      iypout = 21
      itstar = neq + 21
      iyp = 1 + itstar
      iyy = neq + iyp
      iwt = neq + iyy
      ip = neq + iwt
      iphi = neq + ip
      ialpha = (neq*16) + iphi
      ibeta = 12 + ialpha
      ipsi = 12 + ibeta
      iv = 12 + ipsi
      iw = 12 + iv
      isig = 12 + iw
      ig = 13 + isig
      igi = 13 + ig
      ixold = 11 + igi
      ihold = 1 + ixold
      itold = 1 + ihold
      idelsn = 1 + itold
      itwou = 1 + idelsn
      ifouru = 1 + itwou
c
      rwork(itstar) = t
      if (info(1) .eq. 0) go to 50
      start = iwork(21) .ne. (-1)
      phase1 = iwork(22) .ne. (-1)
      nornd = iwork(23) .ne. (-1)
      stiff = iwork(24) .ne. (-1)
      intout = iwork(25) .ne. (-1)
c
 50   call ddes(df,neq,t,y,tout,info,rtol,atol,idid,rwork(iypout),
     1         rwork(iyp),rwork(iyy),rwork(iwt),rwork(ip),rwork(iphi),
     2         rwork(ialpha),rwork(ibeta),rwork(ipsi),rwork(iv),
     3         rwork(iw),rwork(isig),rwork(ig),rwork(igi),rwork(11),
     4         rwork(12),rwork(13),rwork(ixold),rwork(ihold),
     5         rwork(itold),rwork(idelsn),rwork(1),rwork(itwou),
     5         rwork(ifouru),start,phase1,nornd,stiff,intout,iwork(26),
     6         iwork(27),iwork(28),iwork(29),iwork(30),iwork(31),
     7         iwork(32),iwork(33),iwork(34),iwork(35),iwork(45),
     8         rpar,ipar)
c
      iwork(21) = -1
      if (start) iwork(21) = 1
      iwork(22) = -1
      if (phase1) iwork(22) = 1
      iwork(23) = -1
      if (nornd) iwork(23) = 1
      iwork(24) = -1
      if (stiff) iwork(24) = 1
      iwork(25) = -1
      if (intout) iwork(25) = 1
c
      if (idid .ne. (-2)) iwork(liw) = iwork(liw) + 1
      if (t .ne. rwork(itstar)) iwork(liw) = 0
c
      return
      end
