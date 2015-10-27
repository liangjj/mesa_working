*deck sdassl
      subroutine sdassl (res, neq, t, y, yprime, tout, info, rtol, atol,
     *   idid, rwork, lrw, iwork, liw, rpar, ipar, jac)
c***begin prologue  sdassl
c***purpose  this code solves a system of differential/algebraic
c            equations of the form g(t,y,yprime) = 0.
c***library   slatec (dassl)
c***category  i1a2
c***type      single precision (sdassl-s, ddassl-d)
c***keywords  backward differentiation formulas, dassl,
c             differential/algebraic, implicit differential systems
c***author  petzold, linda r., (llnl)
c             computing and mathematics research division
c             lawrence livermore national laboratory
c             l - 316, p.o. box 808,
c             livermore, ca.    94550
c***description
c
c *usage:
c
c      external res, jac
c      integer neq, info(n), idid, lrw, liw, iwork(liw), ipar
c      real t, y(neq), yprime(neq), tout, rtol, atol,
c     *   rwork(lrw), rpar
c
c      call sdassl (res, neq, t, y, yprime, tout, info, rtol, atol,
c     *   idid, rwork, lrw, iwork, liw, rpar, ipar, jac)
c
c
c *arguments:
c
c  res:ext     this is a subroutine which you provide to define the
c              differential/algebraic system.
c
c  neq:in      this is the number of equations to be solved.
c
c  t:inout     this is the current value of the independent variable.
c
c  y(*):inout  this array contains the solution components at t.
c
c  yprime(*):inout  this array contains the derivatives of the solution
c              components at t.
c
c  tout:in     this is a point at which a solution is desired.
c
c  info(n):in  the basic task of the code is to solve the system from t
c              to tout and return an answer at tout.  info is an integer
c              array which is used to communicate exactly how you want
c              this task to be carried out.  (see below for details.)
c              n must be greater than or equal to 15.
c
c  rtol,atol:inout  these quantities represent relative and absolute
c              error tolerances which you provide to indicate how
c              accurately you wish the solution to be computed.  you
c              may choose them to be both scalars or else both vectors.
c              caution:  in fortran 77, a scalar is not the same as an
c                        array of length 1.  some compilers may object
c                        to using scalars for rtol,atol.
c
c  idid:out    this scalar quantity is an indicator reporting what the
c              code did.  you must monitor this integer variable to
c              decide  what action to take next.
c
c  rwork:work  a real work array of length lrw which provides the
c              code with needed storage space.
c
c  lrw:in      the length of rwork.  (see below for required length.)
c
c  iwork:work  an integer work array of length liw which provides the
c              code with needed storage space.
c
c  liw:in      the length of iwork.  (see below for required length.)
c
c  rpar,ipar:in  these are real and integer parameter arrays which
c              you can use for communication between your calling
c              program and the res subroutine (and the jac subroutine)
c
c  jac:ext     this is the name of a subroutine which you may choose
c              to provide for defining a matrix of partial derivatives
c              described below.
c
c  quantities which may be altered by sdassl are:
c     t, y(*), yprime(*), info(1), rtol, atol,
c     idid, rwork(*) and iwork(*)
c
c *description
c
c  subroutine sdassl uses the backward differentiation formulas of
c  orders one through five to solve a system of the above form for y and
c  yprime.  values for y and yprime at the initial time must be given as
c  input.  these values must be consistent, (that is, if t,y,yprime are
c  the given initial values, they must satisfy g(t,y,yprime) = 0.).  the
c  subroutine solves the system from t to tout.  it is easy to continue
c  the solution to get results at additional tout.  this is the interval
c  mode of operation.  intermediate results can also be obtained easily
c  by using the intermediate-output capability.
c
c  the following detailed description is divided into subsections:
c    1. input required for the first call to sdassl.
c    2. output after any return from sdassl.
c    3. what to do to continue the integration.
c    4. error messages.
c
c
c  -------- input -- what to do on the first call to sdassl ------------
c
c  the first call of the code is defined to be the start of each new
c  problem. read through the descriptions of all the following items,
c  provide sufficient storage space for designated arrays, set
c  appropriate variables for the initialization of the problem, and
c  give information about how you want the problem to be solved.
c
c
c  res -- provide a subroutine of the form
c             subroutine res(t,y,yprime,delta,ires,rpar,ipar)
c         to define the system of differential/algebraic
c         equations which is to be solved. for the given values
c         of t,y and yprime, the subroutine should
c         return the residual of the differential/algebraic
c         system
c             delta = g(t,y,yprime)
c         (delta(*) is a vector of length neq which is
c         output for res.)
c
c         subroutine res must not alter t,y or yprime.
c         you must declare the name res in an external
c         statement in your program that calls sdassl.
c         you must dimension y,yprime and delta in res.
c
c         ires is an integer flag which is always equal to
c         zero on input. subroutine res should alter ires
c         only if it encounters an illegal value of y or
c         a stop condition. set ires = -1 if an input value
c         is illegal, and sdassl will try to solve the problem
c         without getting ires = -1. if ires = -2, sdassl
c         will return control to the calling program
c         with idid = -11.
c
c         rpar and ipar are real and integer parameter arrays which
c         you can use for communication between your calling program
c         and subroutine res. they are not altered by sdassl. if you
c         do not need rpar or ipar, ignore these parameters by treat-
c         ing them as dummy arguments. if you do choose to use them,
c         dimension them in your calling program and in res as arrays
c         of appropriate length.
c
c  neq -- set it to the number of differential equations.
c         (neq .ge. 1)
c
c  t -- set it to the initial point of the integration.
c         t must be defined as a variable.
c
c  y(*) -- set this vector to the initial values of the neq solution
c         components at the initial point. you must dimension y of
c         length at least neq in your calling program.
c
c  yprime(*) -- set this vector to the initial values of the neq
c         first derivatives of the solution components at the initial
c         point.  you must dimension yprime at least neq in your
c         calling program. if you do not know initial values of some
c         of the solution components, see the explanation of info(11).
c
c  tout -- set it to the first point at which a solution
c         is desired. you can not take tout = t.
c         integration either forward in t (tout .gt. t) or
c         backward in t (tout .lt. t) is permitted.
c
c         the code advances the solution from t to tout using
c         step sizes which are automatically selected so as to
c         achieve the desired accuracy. if you wish, the code will
c         return with the solution and its derivative at
c         intermediate steps (intermediate-output mode) so that
c         you can monitor them, but you still must provide tout in
c         accord with the basic aim of the code.
c
c         the first step taken by the code is a critical one
c         because it must reflect how fast the solution changes near
c         the initial point. the code automatically selects an
c         initial step size which is practically always suitable for
c         the problem. by using the fact that the code will not step
c         past tout in the first step, you could, if necessary,
c         restrict the length of the initial step size.
c
c         for some problems it may not be permissible to integrate
c         past a point tstop because a discontinuity occurs there
c         or the solution or its derivative is not defined beyond
c         tstop. when you have declared a tstop point (see info(4)
c         and rwork(1)), you have told the code not to integrate
c         past tstop. in this case any tout beyond tstop is invalid
c         input.
c
c  info(*) -- use the info array to give the code more details about
c         how you want your problem solved.  this array should be
c         dimensioned of length 15, though sdassl uses only the first
c         eleven entries.  you must respond to all of the following
c         items, which are arranged as questions.  the simplest use
c         of the code corresponds to answering all questions as yes,
c         i.e. setting all entries of info to 0.
c
c       info(1) - this parameter enables the code to initialize
c              itself. you must set it to indicate the start of every
c              new problem.
c
c          **** is this the first call for this problem ...
c                yes - set info(1) = 0
c                 no - not applicable here.
c                      see below for continuation calls.  ****
c
c       info(2) - how much accuracy you want of your solution
c              is specified by the error tolerances rtol and atol.
c              the simplest use is to take them both to be scalars.
c              to obtain more flexibility, they can both be vectors.
c              the code must be told your choice.
c
c          **** are both error tolerances rtol, atol scalars ...
c                yes - set info(2) = 0
c                      and input scalars for both rtol and atol
c                 no - set info(2) = 1
c                      and input arrays for both rtol and atol ****
c
c       info(3) - the code integrates from t in the direction
c              of tout by steps. if you wish, it will return the
c              computed solution and derivative at the next
c              intermediate step (the intermediate-output mode) or
c              tout, whichever comes first. this is a good way to
c              proceed if you want to see the behavior of the solution.
c              if you must have solutions at a great many specific
c              tout points, this code will compute them efficiently.
c
c          **** do you want the solution only at
c                tout (and not at the next intermediate step) ...
c                 yes - set info(3) = 0
c                  no - set info(3) = 1 ****
c
c       info(4) - to handle solutions at a great many specific
c              values tout efficiently, this code may integrate past
c              tout and interpolate to obtain the result at tout.
c              sometimes it is not possible to integrate beyond some
c              point tstop because the equation changes there or it is
c              not defined past tstop. then you must tell the code
c              not to go past.
c
c           **** can the integration be carried out without any
c                restrictions on the independent variable t ...
c                 yes - set info(4)=0
c                  no - set info(4)=1
c                       and define the stopping point tstop by
c                       setting rwork(1)=tstop ****
c
c       info(5) - to solve differential/algebraic problems it is
c              necessary to use a matrix of partial derivatives of the
c              system of differential equations. if you do not
c              provide a subroutine to evaluate it analytically (see
c              description of the item jac in the call list), it will
c              be approximated by numerical differencing in this code.
c              although it is less trouble for you to have the code
c              compute partial derivatives by numerical differencing,
c              the solution will be more reliable if you provide the
c              derivatives via jac. sometimes numerical differencing
c              is cheaper than evaluating derivatives in jac and
c              sometimes it is not - this depends on your problem.
c
c           **** do you want the code to evaluate the partial
c                derivatives automatically by numerical differences ...
c                   yes - set info(5)=0
c                    no - set info(5)=1
c                  and provide subroutine jac for evaluating the
c                  matrix of partial derivatives ****
c
c       info(6) - sdassl will perform much better if the matrix of
c              partial derivatives, dg/dy + cj*dg/dyprime,
c              (here cj is a scalar determined by sdassl)
c              is banded and the code is told this. in this
c              case, the storage needed will be greatly reduced,
c              numerical differencing will be performed much cheaper,
c              and a number of important algorithms will execute much
c              faster. the differential equation is said to have
c              half-bandwidths ml (lower) and mu (upper) if equation i
c              involves only unknowns y(j) with
c                             i-ml .le. j .le. i+mu
c              for all i=1,2,...,neq. thus, ml and mu are the widths
c              of the lower and upper parts of the band, respectively,
c              with the main diagonal being excluded. if you do not
c              indicate that the equation has a banded matrix of partial
c              derivatives, the code works with a full matrix of neq**2
c              elements (stored in the conventional way). computations
c              with banded matrices cost less time and storage than with
c              full matrices if 2*ml+mu .lt. neq. if you tell the
c              code that the matrix of partial derivatives has a banded
c              structure and you want to provide subroutine jac to
c              compute the partial derivatives, then you must be careful
c              to store the elements of the matrix in the special form
c              indicated in the description of jac.
c
c          **** do you want to solve the problem using a full
c               (dense) matrix (and not a special banded
c               structure) ...
c                yes - set info(6)=0
c                 no - set info(6)=1
c                       and provide the lower (ml) and upper (mu)
c                       bandwidths by setting
c                       iwork(1)=ml
c                       iwork(2)=mu ****
c
c
c        info(7) -- you can specify a maximum (absolute value of)
c              stepsize, so that the code
c              will avoid passing over very
c              large regions.
c
c          ****  do you want the code to decide
c                on its own maximum stepsize?
c                yes - set info(7)=0
c                 no - set info(7)=1
c                      and define hmax by setting
c                      rwork(2)=hmax ****
c
c        info(8) -- differential/algebraic problems
c              may occasionally suffer from
c              severe scaling difficulties on the
c              first step. if you know a great deal
c              about the scaling of your problem, you can
c              help to alleviate this problem by
c              specifying an initial stepsize ho.
c
c          ****  do you want the code to define
c                its own initial stepsize?
c                yes - set info(8)=0
c                 no - set info(8)=1
c                      and define ho by setting
c                      rwork(3)=ho ****
c
c        info(9) -- if storage is a severe problem,
c              you can save some locations by
c              restricting the maximum order maxord.
c              the default value is 5. for each
c              order decrease below 5, the code
c              requires neq fewer locations, however
c              it is likely to be slower. in any
c              case, you must have 1 .le. maxord .le. 5
c          ****  do you want the maximum order to
c                default to 5?
c                yes - set info(9)=0
c                 no - set info(9)=1
c                      and define maxord by setting
c                      iwork(3)=maxord ****
c
c        info(10) --if you know that the solutions to your equations
c               will always be nonnegative, it may help to set this
c               parameter. however, it is probably best to
c               try the code without using this option first,
c               and only to use this option if that doesn't
c               work very well.
c           ****  do you want the code to solve the problem without
c                 invoking any special nonnegativity constraints?
c                  yes - set info(10)=0
c                   no - set info(10)=1
c
c        info(11) --sdassl normally requires the initial t,
c               y, and yprime to be consistent. that is,
c               you must have g(t,y,yprime) = 0 at the initial
c               time. if you do not know the initial
c               derivative precisely, you can let sdassl try
c               to compute it.
c          ****   are the initial t, y, yprime consistent?
c                 yes - set info(11) = 0
c                  no - set info(11) = 1,
c                       and set yprime to an initial approximation
c                       to yprime.  (if you have no idea what
c                       yprime should be, set it to zero. note
c                       that the initial y should be such
c                       that there must exist a yprime so that
c                       g(t,y,yprime) = 0.)
c
c  rtol, atol -- you must assign relative (rtol) and absolute (atol
c         error tolerances to tell the code how accurately you
c         want the solution to be computed.  they must be defined
c         as variables because the code may change them.  you
c         have two choices --
c               both rtol and atol are scalars. (info(2)=0)
c               both rtol and atol are vectors. (info(2)=1)
c         in either case all components must be non-negative.
c
c         the tolerances are used by the code in a local error
c         test at each step which requires roughly that
c               abs(local error) .le. rtol*abs(y)+atol
c         for each vector component.
c         (more specifically, a root-mean-square norm is used to
c         measure the size of vectors, and the error test uses the
c         magnitude of the solution at the beginning of the step.)
c
c         the true (global) error is the difference between the
c         true solution of the initial value problem and the
c         computed approximation.  practically all present day
c         codes, including this one, control the local error at
c         each step and do not even attempt to control the global
c         error directly.
c         usually, but not always, the true accuracy of the
c         computed y is comparable to the error tolerances. this
c         code will usually, but not always, deliver a more
c         accurate solution if you reduce the tolerances and
c         integrate again.  by comparing two such solutions you
c         can get a fairly reliable idea of the true error in the
c         solution at the bigger tolerances.
c
c         setting atol=0. results in a pure relative error test on
c         that component.  setting rtol=0. results in a pure
c         absolute error test on that component.  a mixed test
c         with non-zero rtol and atol corresponds roughly to a
c         relative error test when the solution component is much
c         bigger than atol and to an absolute error test when the
c         solution component is smaller than the threshhold atol.
c
c         the code will not attempt to compute a solution at an
c         accuracy unreasonable for the machine being used.  it will
c         advise you if you ask for too much accuracy and inform
c         you as to the maximum accuracy it believes possible.
c
c  rwork(*) --  dimension this real work array of length lrw in your
c         calling program.
c
c  lrw -- set it to the declared length of the rwork array.
c               you must have
c                    lrw .ge. 40+(maxord+4)*neq+neq**2
c               for the full (dense) jacobian case (when info(6)=0), or
c                    lrw .ge. 40+(maxord+4)*neq+(2*ml+mu+1)*neq
c               for the banded user-defined jacobian case
c               (when info(5)=1 and info(6)=1), or
c                     lrw .ge. 40+(maxord+4)*neq+(2*ml+mu+1)*neq
c                           +2*(neq/(ml+mu+1)+1)
c               for the banded finite-difference-generated jacobian case
c               (when info(5)=0 and info(6)=1)
c
c  iwork(*) --  dimension this integer work array of length liw in
c         your calling program.
c
c  liw -- set it to the declared length of the iwork array.
c               you must have liw .ge. 20+neq
c
c  rpar, ipar -- these are parameter arrays, of real and integer
c         type, respectively.  you can use them for communication
c         between your program that calls sdassl and the
c         res subroutine (and the jac subroutine).  they are not
c         altered by sdassl.  if you do not need rpar or ipar,
c         ignore these parameters by treating them as dummy
c         arguments.  if you do choose to use them, dimension
c         them in your calling program and in res (and in jac)
c         as arrays of appropriate length.
c
c  jac -- if you have set info(5)=0, you can ignore this parameter
c         by treating it as a dummy argument.  otherwise, you must
c         provide a subroutine of the form
c             subroutine jac(t,y,yprime,pd,cj,rpar,ipar)
c         to define the matrix of partial derivatives
c             pd=dg/dy+cj*dg/dyprime
c         cj is a scalar which is input to jac.
c         for the given values of t,y,yprime, the
c         subroutine must evaluate the non-zero partial
c         derivatives for each equation and each solution
c         component, and store these values in the
c         matrix pd.  the elements of pd are set to zero
c         before each call to jac so only non-zero elements
c         need to be defined.
c
c         subroutine jac must not alter t,y,(*),yprime(*), or cj.
c         you must declare the name jac in an external statement in
c         your program that calls sdassl.  you must dimension y,
c         yprime and pd in jac.
c
c         the way you must store the elements into the pd matrix
c         depends on the structure of the matrix which you
c         indicated by info(6).
c               *** info(6)=0 -- full (dense) matrix ***
c                   give pd a first dimension of neq.
c                   when you evaluate the (non-zero) partial derivative
c                   of equation i with respect to variable j, you must
c                   store it in pd according to
c                   pd(i,j) = "dg(i)/dy(j)+cj*dg(i)/dyprime(j)"
c               *** info(6)=1 -- banded jacobian with ml lower and mu
c                   upper diagonal bands (refer to info(6) description
c                   of ml and mu) ***
c                   give pd a first dimension of 2*ml+mu+1.
c                   when you evaluate the (non-zero) partial derivative
c                   of equation i with respect to variable j, you must
c                   store it in pd according to
c                   irow = i - j + ml + mu + 1
c                   pd(irow,j) = "dg(i)/dy(j)+cj*dg(i)/dyprime(j)"
c
c         rpar and ipar are real and integer parameter arrays
c         which you can use for communication between your calling
c         program and your jacobian subroutine jac. they are not
c         altered by sdassl. if you do not need rpar or ipar,
c         ignore these parameters by treating them as dummy
c         arguments. if you do choose to use them, dimension
c         them in your calling program and in jac as arrays of
c         appropriate length.
c
c
c  optionally replaceable norm routine:
c
c     sdassl uses a weighted norm sdanrm to measure the size
c     of vectors such as the estimated error in each step.
c     a function subprogram
c       real function sdanrm(neq,v,wt,rpar,ipar)
c       dimension v(neq),wt(neq)
c     is used to define this norm. here, v is the vector
c     whose norm is to be computed, and wt is a vector of
c     weights.  a sdanrm routine has been included with sdassl
c     which computes the weighted root-mean-square norm
c     given by
c       sdanrm=sqrt((1/neq)*sum(v(i)/wt(i))**2)
c     this norm is suitable for most problems. in some
c     special cases, it may be more convenient and/or
c     efficient to define your own norm by writing a function
c     subprogram to be called instead of sdanrm. this should,
c     however, be attempted only after careful thought and
c     consideration.
c
c
c  -------- output -- after any return from sdassl ---------------------
c
c  the principal aim of the code is to return a computed solution at
c  tout, although it is also possible to obtain intermediate results
c  along the way. to find out whether the code achieved its goal
c  or if the integration process was interrupted before the task was
c  completed, you must check the idid parameter.
c
c
c  t -- the solution was successfully advanced to the
c               output value of t.
c
c  y(*) -- contains the computed solution approximation at t.
c
c  yprime(*) -- contains the computed derivative
c               approximation at t.
c
c  idid -- reports what the code did.
c
c                     *** task completed ***
c                reported by positive values of idid
c
c           idid = 1 -- a step was successfully taken in the
c                   intermediate-output mode. the code has not
c                   yet reached tout.
c
c           idid = 2 -- the integration to tstop was successfully
c                   completed (t=tstop) by stepping exactly to tstop.
c
c           idid = 3 -- the integration to tout was successfully
c                   completed (t=tout) by stepping past tout.
c                   y(*) is obtained by interpolation.
c                   yprime(*) is obtained by interpolation.
c
c                    *** task interrupted ***
c                reported by negative values of idid
c
c           idid = -1 -- a large amount of work has been expended.
c                   (about 500 steps)
c
c           idid = -2 -- the error tolerances are too stringent.
c
c           idid = -3 -- the local error test cannot be satisfied
c                   because you specified a zero component in atol
c                   and the corresponding computed solution
c                   component is zero. thus, a pure relative error
c                   test is impossible for this component.
c
c           idid = -6 -- sdassl had repeated error test
c                   failures on the last attempted step.
c
c           idid = -7 -- the corrector could not converge.
c
c           idid = -8 -- the matrix of partial derivatives
c                   is singular.
c
c           idid = -9 -- the corrector could not converge.
c                   there were repeated error test failures
c                   in this step.
c
c           idid =-10 -- the corrector could not converge
c                   because ires was equal to minus one.
c
c           idid =-11 -- ires equal to -2 was encountered
c                   and control is being returned to the
c                   calling program.
c
c           idid =-12 -- sdassl failed to compute the initial
c                   yprime.
c
c
c
c           idid = -13,..,-32 -- not applicable for this code
c
c                    *** task terminated ***
c                reported by the value of idid=-33
c
c           idid = -33 -- the code has encountered trouble from which
c                   it cannot recover. a message is printed
c                   explaining the trouble and control is returned
c                   to the calling program. for example, this occurs
c                   when invalid input is detected.
c
c  rtol, atol -- these quantities remain unchanged except when
c               idid = -2. in this case, the error tolerances have been
c               increased by the code to values which are estimated to
c               be appropriate for continuing the integration. however,
c               the reported solution at t was obtained using the input
c               values of rtol and atol.
c
c  rwork, iwork -- contain information which is usually of no
c               interest to the user but necessary for subsequent calls.
c               however, you may find use for
c
c               rwork(3)--which contains the step size h to be
c                       attempted on the next step.
c
c               rwork(4)--which contains the current value of the
c                       independent variable, i.e., the farthest point
c                       integration has reached. this will be different
c                       from t only when interpolation has been
c                       performed (idid=3).
c
c               rwork(7)--which contains the stepsize used
c                       on the last successful step.
c
c               iwork(7)--which contains the order of the method to
c                       be attempted on the next step.
c
c               iwork(8)--which contains the order of the method used
c                       on the last step.
c
c               iwork(11)--which contains the number of steps taken so
c                        far.
c
c               iwork(12)--which contains the number of calls to res
c                        so far.
c
c               iwork(13)--which contains the number of evaluations of
c                        the matrix of partial derivatives needed so
c                        far.
c
c               iwork(14)--which contains the total number
c                        of error test failures so far.
c
c               iwork(15)--which contains the total number
c                        of convergence test failures so far.
c                        (includes singular iteration matrix
c                        failures.)
c
c
c  -------- input -- what to do to continue the integration ------------
c                    (calls after the first)
c
c  this code is organized so that subsequent calls to continue the
c  integration involve little (if any) additional effort on your
c  part. you must monitor the idid parameter in order to determine
c  what to do next.
c
c  recalling that the principal task of the code is to integrate
c  from t to tout (the interval mode), usually all you will need
c  to do is specify a new tout upon reaching the current tout.
c
c  do not alter any quantity not specifically permitted below,
c  in particular do not alter neq,t,y(*),yprime(*),rwork(*),iwork(*)
c  or the differential equation in subroutine res. any such
c  alteration constitutes a new problem and must be treated as such,
c  i.e., you must start afresh.
c
c  you cannot change from vector to scalar error control or vice
c  versa (info(2)), but you can change the size of the entries of
c  rtol, atol. increasing a tolerance makes the equation easier
c  to integrate. decreasing a tolerance will make the equation
c  harder to integrate and should generally be avoided.
c
c  you can switch from the intermediate-output mode to the
c  interval mode (info(3)) or vice versa at any time.
c
c  if it has been necessary to prevent the integration from going
c  past a point tstop (info(4), rwork(1)), keep in mind that the
c  code will not integrate to any tout beyond the currently
c  specified tstop. once tstop has been reached you must change
c  the value of tstop or set info(4)=0. you may change info(4)
c  or tstop at any time but you must supply the value of tstop in
c  rwork(1) whenever you set info(4)=1.
c
c  do not change info(5), info(6), iwork(1), or iwork(2)
c  unless you are going to restart the code.
c
c                 *** following a completed task ***
c  if
c     idid = 1, call the code again to continue the integration
c                  another step in the direction of tout.
c
c     idid = 2 or 3, define a new tout and call the code again.
c                  tout must be different from t. you cannot change
c                  the direction of integration without restarting.
c
c                 *** following an interrupted task ***
c               to show the code that you realize the task was
c               interrupted and that you want to continue, you
c               must take appropriate action and set info(1) = 1
c  if
c    idid = -1, the code has taken about 500 steps.
c                  if you want to continue, set info(1) = 1 and
c                  call the code again. an additional 500 steps
c                  will be allowed.
c
c    idid = -2, the error tolerances rtol, atol have been
c                  increased to values the code estimates appropriate
c                  for continuing. you may want to change them
c                  yourself. if you are sure you want to continue
c                  with relaxed error tolerances, set info(1)=1 and
c                  call the code again.
c
c    idid = -3, a solution component is zero and you set the
c                  corresponding component of atol to zero. if you
c                  are sure you want to continue, you must first
c                  alter the error criterion to use positive values
c                  for those components of atol corresponding to zero
c                  solution components, then set info(1)=1 and call
c                  the code again.
c
c    idid = -4,-5  --- cannot occur with this code.
c
c    idid = -6, repeated error test failures occurred on the
c                  last attempted step in sdassl. a singularity in the
c                  solution may be present. if you are absolutely
c                  certain you want to continue, you should restart
c                  the integration. (provide initial values of y and
c                  yprime which are consistent)
c
c    idid = -7, repeated convergence test failures occurred
c                  on the last attempted step in sdassl. an inaccurate
c                  or ill-conditioned jacobian may be the problem. if
c                  you are absolutely certain you want to continue, you
c                  should restart the integration.
c
c    idid = -8, the matrix of partial derivatives is singular.
c                  some of your equations may be redundant.
c                  sdassl cannot solve the problem as stated.
c                  it is possible that the redundant equations
c                  could be removed, and then sdassl could
c                  solve the problem. it is also possible
c                  that a solution to your problem either
c                  does not exist or is not unique.
c
c    idid = -9, sdassl had multiple convergence test
c                  failures, preceded by multiple error
c                  test failures, on the last attempted step.
c                  it is possible that your problem
c                  is ill-posed, and cannot be solved
c                  using this code. or, there may be a
c                  discontinuity or a singularity in the
c                  solution. if you are absolutely certain
c                  you want to continue, you should restart
c                  the integration.
c
c    idid =-10, sdassl had multiple convergence test failures
c                  because ires was equal to minus one.
c                  if you are absolutely certain you want
c                  to continue, you should restart the
c                  integration.
c
c    idid =-11, ires=-2 was encountered, and control is being
c                  returned to the calling program.
c
c    idid =-12, sdassl failed to compute the initial yprime.
c                  this could happen because the initial
c                  approximation to yprime was not very good, or
c                  if a yprime consistent with the initial y
c                  does not exist. the problem could also be caused
c                  by an inaccurate or singular iteration matrix.
c
c    idid = -13,..,-32  --- cannot occur with this code.
c
c
c                 *** following a terminated task ***
c
c  if idid= -33, you cannot continue the solution of this problem.
c                  an attempt to do so will result in your
c                  run being terminated.
c
c
c  -------- error messages ---------------------------------------------
c
c      the slatec error print routine xermsg is called in the event of
c   unsuccessful completion of a task.  most of these are treated as
c   "recoverable errors", which means that (unless the user has directed
c   otherwise) control will be returned to the calling program for
c   possible action after the message has been printed.
c
c   in the event of a negative value of idid other than -33, an appro-
c   priate message is printed and the "error number" printed by xermsg
c   is the value of idid.  there are quite a number of illegal input
c   errors that can lead to a returned value idid=-33.  the conditions
c   and their printed "error numbers" are as follows:
c
c   error number       condition
c
c        1       some element of info vector is not zero or one.
c        2       neq .le. 0
c        3       maxord not in range.
c        4       lrw is less than the required length for rwork.
c        5       liw is less than the required length for iwork.
c        6       some element of rtol is .lt. 0
c        7       some element of atol is .lt. 0
c        8       all elements of rtol and atol are zero.
c        9       info(4)=1 and tstop is behind tout.
c       10       hmax .lt. 0.0
c       11       tout is behind t.
c       12       info(8)=1 and h0=0.0
c       13       some element of wt is .le. 0.0
c       14       tout is too close to t to start integration.
c       15       info(4)=1 and tstop is behind t.
c       16       --( not used in this version )--
c       17       ml illegal.  either .lt. 0 or .gt. neq
c       18       mu illegal.  either .lt. 0 or .gt. neq
c       19       tout = t.
c
c   if sdassl is called again without any action taken to remove the
c   cause of an unsuccessful return, xermsg will be called with a fatal
c   error flag, which will cause unconditional termination of the
c   program.  there are two such fatal errors:
c
c   error number -998:  the last step was terminated with a negative
c       value of idid other than -33, and no appropriate action was
c       taken.
c
c   error number -999:  the previous call was terminated because of
c       illegal input (idid=-33) and there is illegal input in the
c       present call, as well.  (suspect infinite loop.)
c
c  ---------------------------------------------------------------------
c
c***references  a description of dassl: a differential/algebraic
c                 system solver, l. r. petzold, sand82-8637,
c                 sandia national laboratories, september 1982.
c***routines called  r1mach, sdaini, sdanrm, sdastp, sdatrp, sdawts,
c                    xermsg
c***revision history  (yymmdd)
c   830315  date written
c   880387  code changes made.  all common statements have been
c           replaced by a data statement, which defines pointers into
c           rwork, and parameter statements which define pointers
c           into iwork.  as well the documentation has gone through
c           grammatical changes.
c   881005  the prologue has been changed to mixed case.
c           the subordinate routines had revision dates changed to
c           this date, although the documentation for these routines
c           is all upper case.  no code changes.
c   890511  code changes made.  the data statement in the declaration
c           section of sdassl was replaced with a parameter
c           statement.  also the statement s = 100.e0 was removed
c           from the top of the newton iteration in sdastp.
c           the subordinate routines had revision dates changed to
c           this date.
c   890517  the revision date syntax was replaced with the revision
c           history syntax.  also the "deck" comment was added to
c           the top of all subroutines.  these changes are consistent
c           with new slatec guidelines.
c           the subordinate routines had revision dates changed to
c           this date.  no code changes.
c   891013  code changes made.
c           removed all occurrences of float.  all operations
c           are now performed with "mixed-mode" arithmetic.
c           also, specific function names were replaced with generic
c           function names to be consistent with new slatec guidelines.
c           in particular:
c              replaced amin1 with min everywhere.
c              replaced min0 with min everywhere.
c              replaced amax1 with max everywhere.
c              replaced max0 with max everywhere.
c           also replaced revision date with revision history in all
c           subordinate routines.
c   901004  miscellaneous changes to prologue to complete conversion
c           to slatec 4.0 format.  no code changes.  (f.n.fritsch)
c   901009  corrected gams classification code and converted subsidiary
c           routines to 4.0 format.  no code changes.  (f.n.fritsch)
c   901010  converted xerrwv calls to xermsg calls.  (r.clemens, afwl)
c   901019  code changes made.
c           merged slatec 4.0 changes with previous changes made
c           by c. ulrich.  below is a history of the changes made by
c           c. ulrich. (changes in subsidiary routines are implied
c           by this history)
c           891228  bug was found and repaired inside the sdassl
c                   and sdaini routines.  sdaini was incorrectly
c                   returning the initial t with y and yprime
c                   computed at t+h.  the routine now returns t+h
c                   rather than the initial t.
c                   cosmetic changes made to sdastp.
c           900904  three modifications were made to fix a bug (inside
c                   sdassl) re interpolation for continuation calls and
c                   cases where tn is very close to tstop:
c
c                   1) in testing for whether h is too large, just
c                      compare h to (tstop - tn), rather than
c                      (tstop - tn) * (1-4*uround), and set h to
c                      tstop - tn.  this will force sdastp to step
c                      exactly to tstop under certain situations
c                      (i.e. when h returned from sdastp would otherwise
c                      take tn beyond tstop).
c
c                   2) inside the sdastp loop, interpolate exactly to
c                      tstop if tn is very close to tstop (rather than
c                      interpolating to within roundoff of tstop).
c
c                   3) modified idid description for idid = 2 to say
c                      that the solution is returned by stepping exactly
c                      to tstop, rather than tout.  (in some cases the
c                      solution is actually obtained by extrapolating
c                      over a distance near unit roundoff to tstop,
c                      but this small distance is deemed acceptable in
c                      these circumstances.)
c   901026  added explicit declarations for all variables and minor
c           cosmetic changes to prologue, removed unreferenced labels,
c           and improved xermsg calls.  (fnf)
c   901030  added error messages section and reworked other sections to
c           be of more uniform format.  (fnf)
c   910624  fixed minor bug related to hmax (six lines after label
c           525).  (lrp)
c***end prologue  sdassl
c
c**end
c
c     declare arguments.
c
      integer  neq, info(15), idid, lrw, iwork(*), liw, ipar(*)
      real  t, y(*), yprime(*), tout, rtol(*), atol(*), rwork(*),
     *   rpar(*)
      external  res, jac
c
c     declare externals.
c
      external  r1mach, sdaini, sdanrm, sdastp, sdatrp, sdawts, xermsg
      real  r1mach, sdanrm
c
c     declare local variables.
c
      integer  i, itemp, lalpha, lbeta, lcj, lcjold, lctf, ldelta,
     *   leniw, lenpd, lenrw, le, letf, lgamma, lh, lhmax, lhold, lipvt,
     *   ljcalc, lk, lkold, liwm, lml, lmtype, lmu, lmxord, lnje, lnpd,
     *   lnre, lns, lnst, lnstl, lpd, lphase, lphi, lpsi, lround, ls,
     *   lsigma, ltn, ltstop, lwm, lwt, mband, msave, mxord, npd, ntemp,
     *   nzflg
      real  atoli, h, hmax, hmin, ho, r, rh, rtoli, tdist, tn, tnext,
     *   tstop, uround, ypnorm
      logical  done
c       auxiliary variables for conversion of values to be included in
c       error messages.
      character*8  xern1, xern2
      character*16 xern3, xern4
c
c     set pointers into iwork
      parameter (lml=1, lmu=2, lmxord=3, lmtype=4, lnst=11,
     *  lnre=12, lnje=13, letf=14, lctf=15, lnpd=16,
     *  lipvt=21, ljcalc=5, lphase=6, lk=7, lkold=8,
     *  lns=9, lnstl=10, liwm=1)
c
c     set relative offset into rwork
      parameter (npd=1)
c
c     set pointers into rwork
      parameter (ltstop=1, lhmax=2, lh=3, ltn=4,
     *  lcj=5, lcjold=6, lhold=7, ls=8, lround=9,
     *  lalpha=11, lbeta=17, lgamma=23,
     *  lpsi=29, lsigma=35, ldelta=41)
c
c***first executable statement  sdassl
      if(info(1).ne.0)go to 100
c
c-----------------------------------------------------------------------
c     this block is executed for the initial call only.
c     it contains checking of inputs and initializations.
c-----------------------------------------------------------------------
c
c     first check info array to make sure all elements of info
c     are either zero or one.
      do 10 i=2,11
         if(info(i).ne.0.and.info(i).ne.1)go to 701
10       continue
c
      if(neq.le.0)go to 702
c
c     check and compute maximum order
      mxord=5
      if(info(9).eq.0)go to 20
         mxord=iwork(lmxord)
         if(mxord.lt.1.or.mxord.gt.5)go to 703
20       iwork(lmxord)=mxord
c
c     compute mtype,lenpd,lenrw.check ml and mu.
      if(info(6).ne.0)go to 40
         lenpd=neq**2
         lenrw=40+(iwork(lmxord)+4)*neq+lenpd
         if(info(5).ne.0)go to 30
            iwork(lmtype)=2
            go to 60
30          iwork(lmtype)=1
            go to 60
40    if(iwork(lml).lt.0.or.iwork(lml).ge.neq)go to 717
      if(iwork(lmu).lt.0.or.iwork(lmu).ge.neq)go to 718
      lenpd=(2*iwork(lml)+iwork(lmu)+1)*neq
      if(info(5).ne.0)go to 50
         iwork(lmtype)=5
         mband=iwork(lml)+iwork(lmu)+1
         msave=(neq/mband)+1
         lenrw=40+(iwork(lmxord)+4)*neq+lenpd+2*msave
         go to 60
50       iwork(lmtype)=4
         lenrw=40+(iwork(lmxord)+4)*neq+lenpd
c
c     check lengths of rwork and iwork
60    leniw=20+neq
      iwork(lnpd)=lenpd
      if(lrw.lt.lenrw)go to 704
      if(liw.lt.leniw)go to 705
c
c     check to see that tout is different from t
      if(tout .eq. t)go to 719
c
c     check hmax
      if(info(7).eq.0)go to 70
         hmax=rwork(lhmax)
         if(hmax.le.0.0e0)go to 710
70    continue
c
c     initialize counters
      iwork(lnst)=0
      iwork(lnre)=0
      iwork(lnje)=0
c
      iwork(lnstl)=0
      idid=1
      go to 200
c
c-----------------------------------------------------------------------
c     this block is for continuation calls
c     only. here we check info(1), and if the
c     last step was interrupted we check whether
c     appropriate action was taken.
c-----------------------------------------------------------------------
c
100   continue
      if(info(1).eq.1)go to 110
      if(info(1).ne.-1)go to 701
c
c     if we are here, the last step was interrupted
c     by an error condition from sdastp, and
c     appropriate action was not taken. this
c     is a fatal error.
      write (xern1, '(i8)') idid
      call xermsg ('slatec', 'sdassl',
     *   'the last step terminated with a negative value of idid = ' //
     *   xern1 // ' and no appropriate action was taken.  ' //
     *   'run terminated', -998, 2)
      return
110   continue
      iwork(lnstl)=iwork(lnst)
c
c-----------------------------------------------------------------------
c     this block is executed on all calls.
c     the error tolerance parameters are
c     checked, and the work array pointers
c     are set.
c-----------------------------------------------------------------------
c
200   continue
c     check rtol,atol
      nzflg=0
      rtoli=rtol(1)
      atoli=atol(1)
      do 210 i=1,neq
         if(info(2).eq.1)rtoli=rtol(i)
         if(info(2).eq.1)atoli=atol(i)
         if(rtoli.gt.0.0e0.or.atoli.gt.0.0e0)nzflg=1
         if(rtoli.lt.0.0e0)go to 706
         if(atoli.lt.0.0e0)go to 707
210      continue
      if(nzflg.eq.0)go to 708
c
c     set up rwork storage.iwork storage is fixed
c     in data statement.
      le=ldelta+neq
      lwt=le+neq
      lphi=lwt+neq
      lpd=lphi+(iwork(lmxord)+1)*neq
      lwm=lpd
      ntemp=npd+iwork(lnpd)
      if(info(1).eq.1)go to 400
c
c-----------------------------------------------------------------------
c     this block is executed on the initial call
c     only. set the initial step size, and
c     the error weight vector, and phi.
c     compute initial yprime, if necessary.
c-----------------------------------------------------------------------
c
      tn=t
      idid=1
c
c     set error weight vector wt
      call sdawts(neq,info(2),rtol,atol,y,rwork(lwt),rpar,ipar)
      do 305 i = 1,neq
         if(rwork(lwt+i-1).le.0.0e0) go to 713
305      continue
c
c     compute unit roundoff and hmin
      uround = r1mach(4)
      rwork(lround) = uround
      hmin = 4.0e0*uround*max(abs(t),abs(tout))
c
c     check initial interval to see that it is long enough
      tdist = abs(tout - t)
      if(tdist .lt. hmin) go to 714
c
c     check ho, if this was input
      if (info(8) .eq. 0) go to 310
         ho = rwork(lh)
         if ((tout - t)*ho .lt. 0.0e0) go to 711
         if (ho .eq. 0.0e0) go to 712
         go to 320
310    continue
c
c     compute initial stepsize, to be used by either
c     sdastp or sdaini, depending on info(11)
      ho = 0.001e0*tdist
      ypnorm = sdanrm(neq,yprime,rwork(lwt),rpar,ipar)
      if (ypnorm .gt. 0.5e0/ho) ho = 0.5e0/ypnorm
      ho = sign(ho,tout-t)
c     adjust ho if necessary to meet hmax bound
320   if (info(7) .eq. 0) go to 330
         rh = abs(ho)/rwork(lhmax)
         if (rh .gt. 1.0e0) ho = ho/rh
c     compute tstop, if applicable
330   if (info(4) .eq. 0) go to 340
         tstop = rwork(ltstop)
         if ((tstop - t)*ho .lt. 0.0e0) go to 715
         if ((t + ho - tstop)*ho .gt. 0.0e0) ho = tstop - t
         if ((tstop - tout)*ho .lt. 0.0e0) go to 709
c
c     compute initial derivative, updating tn and y, if applicable
340   if (info(11) .eq. 0) go to 350
      call sdaini(tn,y,yprime,neq,
     *  res,jac,ho,rwork(lwt),idid,rpar,ipar,
     *  rwork(lphi),rwork(ldelta),rwork(le),
     *  rwork(lwm),iwork(liwm),hmin,rwork(lround),
     *  info(10),ntemp)
      if (idid .lt. 0) go to 390
c
c     load h with ho.  store h in rwork(lh)
350   h = ho
      rwork(lh) = h
c
c     load y and h*yprime into phi(*,1) and phi(*,2)
      itemp = lphi + neq
      do 370 i = 1,neq
         rwork(lphi + i - 1) = y(i)
370      rwork(itemp + i - 1) = h*yprime(i)
c
390   go to 500
c
c-------------------------------------------------------
c     this block is for continuation calls only. its
c     purpose is to check stop conditions before
c     taking a step.
c     adjust h if necessary to meet hmax bound
c-------------------------------------------------------
c
400   continue
      uround=rwork(lround)
      done = .false.
      tn=rwork(ltn)
      h=rwork(lh)
      if(info(7) .eq. 0) go to 410
         rh = abs(h)/rwork(lhmax)
         if(rh .gt. 1.0e0) h = h/rh
410   continue
      if(t .eq. tout) go to 719
      if((t - tout)*h .gt. 0.0e0) go to 711
      if(info(4) .eq. 1) go to 430
      if(info(3) .eq. 1) go to 420
      if((tn-tout)*h.lt.0.0e0)go to 490
      call sdatrp(tn,tout,y,yprime,neq,iwork(lkold),
     *  rwork(lphi),rwork(lpsi))
      t=tout
      idid = 3
      done = .true.
      go to 490
420   if((tn-t)*h .le. 0.0e0) go to 490
      if((tn - tout)*h .gt. 0.0e0) go to 425
      call sdatrp(tn,tn,y,yprime,neq,iwork(lkold),
     *  rwork(lphi),rwork(lpsi))
      t = tn
      idid = 1
      done = .true.
      go to 490
425   continue
      call sdatrp(tn,tout,y,yprime,neq,iwork(lkold),
     *  rwork(lphi),rwork(lpsi))
      t = tout
      idid = 3
      done = .true.
      go to 490
430   if(info(3) .eq. 1) go to 440
      tstop=rwork(ltstop)
      if((tn-tstop)*h.gt.0.0e0) go to 715
      if((tstop-tout)*h.lt.0.0e0)go to 709
      if((tn-tout)*h.lt.0.0e0)go to 450
      call sdatrp(tn,tout,y,yprime,neq,iwork(lkold),
     *   rwork(lphi),rwork(lpsi))
      t=tout
      idid = 3
      done = .true.
      go to 490
440   tstop = rwork(ltstop)
      if((tn-tstop)*h .gt. 0.0e0) go to 715
      if((tstop-tout)*h .lt. 0.0e0) go to 709
      if((tn-t)*h .le. 0.0e0) go to 450
      if((tn - tout)*h .gt. 0.0e0) go to 445
      call sdatrp(tn,tn,y,yprime,neq,iwork(lkold),
     *  rwork(lphi),rwork(lpsi))
      t = tn
      idid = 1
      done = .true.
      go to 490
445   continue
      call sdatrp(tn,tout,y,yprime,neq,iwork(lkold),
     *  rwork(lphi),rwork(lpsi))
      t = tout
      idid = 3
      done = .true.
      go to 490
450   continue
c     check whether we are within roundoff of tstop
      if(abs(tn-tstop).gt.100.0e0*uround*
     *   (abs(tn)+abs(h)))go to 460
      call sdatrp(tn,tstop,y,yprime,neq,iwork(lkold),
     *  rwork(lphi),rwork(lpsi))
      idid=2
      t=tstop
      done = .true.
      go to 490
460   tnext=tn+h
      if((tnext-tstop)*h.le.0.0e0)go to 490
      h=tstop-tn
      rwork(lh)=h
c
490   if (done) go to 580
c
c-------------------------------------------------------
c     the next block contains the call to the
c     one-step integrator sdastp.
c     this is a looping point for the integration steps.
c     check for too many steps.
c     update wt.
c     check for too much accuracy requested.
c     compute minimum stepsize.
c-------------------------------------------------------
c
500   continue
c     check for failure to compute initial yprime
      if (idid .eq. -12) go to 527
c
c     check for too many steps
      if((iwork(lnst)-iwork(lnstl)).lt.500)
     *   go to 510
           idid=-1
           go to 527
c
c     update wt
510   call sdawts(neq,info(2),rtol,atol,rwork(lphi),
     *  rwork(lwt),rpar,ipar)
      do 520 i=1,neq
         if(rwork(i+lwt-1).gt.0.0e0)go to 520
           idid=-3
           go to 527
520   continue
c
c     test for too much accuracy requested.
      r=sdanrm(neq,rwork(lphi),rwork(lwt),rpar,ipar)*
     *   100.0e0*uround
      if(r.le.1.0e0)go to 525
c     multiply rtol and atol by r and return
      if(info(2).eq.1)go to 523
           rtol(1)=r*rtol(1)
           atol(1)=r*atol(1)
           idid=-2
           go to 527
523   do 524 i=1,neq
           rtol(i)=r*rtol(i)
524        atol(i)=r*atol(i)
      idid=-2
      go to 527
525   continue
c
c     compute minimum stepsize
      hmin=4.0e0*uround*max(abs(tn),abs(tout))
c
c     test h vs. hmax
      if (info(7) .ne. 0) then
         rh = abs(h)/rwork(lhmax)
         if (rh .gt. 1.0e0) h = h/rh
      endif
c
      call sdastp(tn,y,yprime,neq,
     *   res,jac,h,rwork(lwt),info(1),idid,rpar,ipar,
     *   rwork(lphi),rwork(ldelta),rwork(le),
     *   rwork(lwm),iwork(liwm),
     *   rwork(lalpha),rwork(lbeta),rwork(lgamma),
     *   rwork(lpsi),rwork(lsigma),
     *   rwork(lcj),rwork(lcjold),rwork(lhold),
     *   rwork(ls),hmin,rwork(lround),
     *   iwork(lphase),iwork(ljcalc),iwork(lk),
     *   iwork(lkold),iwork(lns),info(10),ntemp)
527   if(idid.lt.0)go to 600
c
c--------------------------------------------------------
c     this block handles the case of a successful return
c     from sdastp (idid=1).  test for stop conditions.
c--------------------------------------------------------
c
      if(info(4).ne.0)go to 540
           if(info(3).ne.0)go to 530
             if((tn-tout)*h.lt.0.0e0)go to 500
             call sdatrp(tn,tout,y,yprime,neq,
     *         iwork(lkold),rwork(lphi),rwork(lpsi))
             idid=3
             t=tout
             go to 580
530          if((tn-tout)*h.ge.0.0e0)go to 535
             t=tn
             idid=1
             go to 580
535          call sdatrp(tn,tout,y,yprime,neq,
     *         iwork(lkold),rwork(lphi),rwork(lpsi))
             idid=3
             t=tout
             go to 580
540   if(info(3).ne.0)go to 550
      if((tn-tout)*h.lt.0.0e0)go to 542
         call sdatrp(tn,tout,y,yprime,neq,
     *     iwork(lkold),rwork(lphi),rwork(lpsi))
         t=tout
         idid=3
         go to 580
542   if(abs(tn-tstop).le.100.0e0*uround*
     *   (abs(tn)+abs(h)))go to 545
      tnext=tn+h
      if((tnext-tstop)*h.le.0.0e0)go to 500
      h=tstop-tn
      go to 500
545   call sdatrp(tn,tstop,y,yprime,neq,
     *  iwork(lkold),rwork(lphi),rwork(lpsi))
      idid=2
      t=tstop
      go to 580
550   if((tn-tout)*h.ge.0.0e0)go to 555
      if(abs(tn-tstop).le.100.0e0*uround*(abs(tn)+abs(h)))go to 552
      t=tn
      idid=1
      go to 580
552   call sdatrp(tn,tstop,y,yprime,neq,
     *  iwork(lkold),rwork(lphi),rwork(lpsi))
      idid=2
      t=tstop
      go to 580
555   call sdatrp(tn,tout,y,yprime,neq,
     *   iwork(lkold),rwork(lphi),rwork(lpsi))
      t=tout
      idid=3
      go to 580
c
c--------------------------------------------------------
c     all successful returns from sdassl are made from
c     this block.
c--------------------------------------------------------
c
580   continue
      rwork(ltn)=tn
      rwork(lh)=h
      return
c
c-----------------------------------------------------------------------
c     this block handles all unsuccessful
c     returns other than for illegal input.
c-----------------------------------------------------------------------
c
600   continue
      itemp=-idid
      go to (610,620,630,690,690,640,650,660,670,675,
     *  680,685), itemp
c
c     the maximum number of steps was taken before
c     reaching tout
610   write (xern3, '(1p,e15.6)') tn
      call xermsg ('slatec', 'sdassl',
     *   'at current t = ' // xern3 // ' 500 steps taken on this ' //
     *   'call before reaching tout', idid, 1)
      go to 690
c
c     too much accuracy for machine precision
620   write (xern3, '(1p,e15.6)') tn
      call xermsg ('slatec', 'sdassl',
     *   'at t = ' // xern3 // ' too much accuracy requested for ' //
     *   'precision of machine. rtol and atol were increased to ' //
     *   'appropriate values', idid, 1)
      go to 690
c
c     wt(i) .le. 0.0 for some i (not at start of problem)
630   write (xern3, '(1p,e15.6)') tn
      call xermsg ('slatec', 'sdassl',
     *   'at t = ' // xern3 // ' some element of wt has become .le. ' //
     *   '0.0', idid, 1)
      go to 690
c
c     error test failed repeatedly or with h=hmin
640   write (xern3, '(1p,e15.6)') tn
      write (xern4, '(1p,e15.6)') h
      call xermsg ('slatec', 'sdassl',
     *   'at t = ' // xern3 // ' and stepsize h = ' // xern4 //
     *   ' the error test failed repeatedly or with abs(h)=hmin',
     *   idid, 1)
      go to 690
c
c     corrector convergence failed repeatedly or with h=hmin
650   write (xern3, '(1p,e15.6)') tn
      write (xern4, '(1p,e15.6)') h
      call xermsg ('slatec', 'sdassl',
     *   'at t = ' // xern3 // ' and stepsize h = ' // xern4 //
     *   ' the corrector failed to converge repeatedly or with ' //
     *   'abs(h)=hmin', idid, 1)
      go to 690
c
c     the iteration matrix is singular
660   write (xern3, '(1p,e15.6)') tn
      write (xern4, '(1p,e15.6)') h
      call xermsg ('slatec', 'sdassl',
     *   'at t = ' // xern3 // ' and stepsize h = ' // xern4 //
     *   ' the iteration matrix is singular', idid, 1)
      go to 690
c
c     corrector failure preceded by error test failures.
670   write (xern3, '(1p,e15.6)') tn
      write (xern4, '(1p,e15.6)') h
      call xermsg ('slatec', 'sdassl',
     *   'at t = ' // xern3 // ' and stepsize h = ' // xern4 //
     *   ' the corrector could not converge.  also, the error test ' //
     *   'failed repeatedly.', idid, 1)
      go to 690
c
c     corrector failure because ires = -1
675   write (xern3, '(1p,e15.6)') tn
      write (xern4, '(1p,e15.6)') h
      call xermsg ('slatec', 'sdassl',
     *   'at t = ' // xern3 // ' and stepsize h = ' // xern4 //
     *   ' the corrector could not converge because ires was equal ' //
     *   'to minus one', idid, 1)
      go to 690
c
c     failure because ires = -2
680   write (xern3, '(1p,e15.6)') tn
      write (xern4, '(1p,e15.6)') h
      call xermsg ('slatec', 'sdassl',
     *   'at t = ' // xern3 // ' and stepsize h = ' // xern4 //
     *   ' ires was equal to minus two', idid, 1)
      go to 690
c
c     failed to compute initial yprime
685   write (xern3, '(1p,e15.6)') tn
      write (xern4, '(1p,e15.6)') ho
      call xermsg ('slatec', 'sdassl',
     *   'at t = ' // xern3 // ' and stepsize h = ' // xern4 //
     *   ' the initial yprime could not be computed', idid, 1)
      go to 690
c
690   continue
      info(1)=-1
      t=tn
      rwork(ltn)=tn
      rwork(lh)=h
      return
c
c-----------------------------------------------------------------------
c     this block handles all error returns due
c     to illegal input, as detected before calling
c     sdastp. first the error message routine is
c     called. if this happens twice in
c     succession, execution is terminated
c
c-----------------------------------------------------------------------
701   call xermsg ('slatec', 'sdassl',
     *   'some element of info vector is not zero or one', 1, 1)
      go to 750
c
702   write (xern1, '(i8)') neq
      call xermsg ('slatec', 'sdassl',
     *   'neq = ' // xern1 // ' .le. 0', 2, 1)
      go to 750
c
703   write (xern1, '(i8)') mxord
      call xermsg ('slatec', 'sdassl',
     *   'maxord = ' // xern1 // ' not in range', 3, 1)
      go to 750
c
704   write (xern1, '(i8)') lenrw
      write (xern2, '(i8)') lrw
      call xermsg ('slatec', 'sdassl',
     *   'rwork length needed, lenrw = ' // xern1 //
     *   ', exceeds lrw = ' // xern2, 4, 1)
      go to 750
c
705   write (xern1, '(i8)') leniw
      write (xern2, '(i8)') liw
      call xermsg ('slatec', 'sdassl',
     *   'iwork length needed, leniw = ' // xern1 //
     *   ', exceeds liw = ' // xern2, 5, 1)
      go to 750
c
706   call xermsg ('slatec', 'sdassl',
     *   'some element of rtol is .lt. 0', 6, 1)
      go to 750
c
707   call xermsg ('slatec', 'sdassl',
     *   'some element of atol is .lt. 0', 7, 1)
      go to 750
c
708   call xermsg ('slatec', 'sdassl',
     *   'all elements of rtol and atol are zero', 8, 1)
      go to 750
c
709   write (xern3, '(1p,e15.6)') tstop
      write (xern4, '(1p,e15.6)') tout
      call xermsg ('slatec', 'sdassl',
     *   'info(4) = 1 and tstop = ' // xern3 // ' behind tout = ' //
     *   xern4, 9, 1)
      go to 750
c
710   write (xern3, '(1p,e15.6)') hmax
      call xermsg ('slatec', 'sdassl',
     *   'hmax = ' // xern3 // ' .lt. 0.0', 10, 1)
      go to 750
c
711   write (xern3, '(1p,e15.6)') tout
      write (xern4, '(1p,e15.6)') t
      call xermsg ('slatec', 'sdassl',
     *   'tout = ' // xern3 // ' behind t = ' // xern4, 11, 1)
      go to 750
c
712   call xermsg ('slatec', 'sdassl',
     *   'info(8)=1 and h0=0.0', 12, 1)
      go to 750
c
713   call xermsg ('slatec', 'sdassl',
     *   'some element of wt is .le. 0.0', 13, 1)
      go to 750
c
714   write (xern3, '(1p,e15.6)') tout
      write (xern4, '(1p,e15.6)') t
      call xermsg ('slatec', 'sdassl',
     *   'tout = ' // xern3 // ' too close to t = ' // xern4 //
     *   ' to start integration', 14, 1)
      go to 750
c
715   write (xern3, '(1p,e15.6)') tstop
      write (xern4, '(1p,e15.6)') t
      call xermsg ('slatec', 'sdassl',
     *   'info(4)=1 and tstop = ' // xern3 // ' behind t = ' // xern4,
     *   15, 1)
      go to 750
c
717   write (xern1, '(i8)') iwork(lml)
      call xermsg ('slatec', 'sdassl',
     *   'ml = ' // xern1 // ' illegal.  either .lt. 0 or .gt. neq',
     *   17, 1)
      go to 750
c
718   write (xern1, '(i8)') iwork(lmu)
      call xermsg ('slatec', 'sdassl',
     *   'mu = ' // xern1 // ' illegal.  either .lt. 0 or .gt. neq',
     *   18, 1)
      go to 750
c
719   write (xern3, '(1p,e15.6)') tout
      call xermsg ('slatec', 'sdassl',
     *  'tout = t = ' // xern3, 19, 1)
      go to 750
c
750   idid=-33
      if(info(1).eq.-1) then
         call xermsg ('slatec', 'sdassl',
     *      'repeated occurrences of illegal input$$' //
     *      'run terminated. apparent infinite loop', -999, 2)
      endif
c
      info(1)=-1
      return
c-----------end of subroutine sdassl------------------------------------
      end
