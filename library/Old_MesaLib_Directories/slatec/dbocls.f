*deck dbocls
      subroutine dbocls (w, mdw, mcon, mrows, ncols, bl, bu, ind, iopt,
     +   x, rnormc, rnorm, mode, rw, iw)
c***begin prologue  dbocls
c***purpose  solve the bounded and constrained least squares
c            problem consisting of solving the equation
c                      e*x = f  (in the least squares sense)
c             subject to the linear constraints
c                            c*x = y.
c***library   slatec
c***category  k1a2a, g2e, g2h1, g2h2
c***type      double precision (sbocls-s, dbocls-d)
c***keywords  bounds, constraints, inequality, least squares, linear
c***author  hanson, r. j., (snla)
c***description
c
c   **** all input and output real variables are double precision ****
c
c     this subprogram solves the bounded and constrained least squares
c     problem. the problem statement is:
c
c     solve e*x = f (least squares sense), subject to constraints
c     c*x=y.
c
c     in this formulation both x and y are unknowns, and both may
c     have bounds on any of their components.  this formulation
c     of the problem allows the user to have equality and inequality
c     constraints as well as simple bounds on the solution components.
c
c     this constrained linear least squares subprogram solves e*x=f
c     subject to c*x=y, where e is mrows by ncols, c is mcon by ncols.
c
c      the user must have dimension statements of the form
c
c      dimension w(mdw,ncols+mcon+1), bl(ncols+mcon), bu(ncols+mcon),
c     * x(2*(ncols+mcon)+2+nx), rw(6*ncols+5*mcon)
c       integer ind(ncols+mcon), iopt(17+ni), iw(2*(ncols+mcon))
c
c     (here nx=number of extra locations required for the options; nx=0
c     if no options are in use. also ni=number of extra locations
c     for options 1-9.)
c
c    input
c    -----
c
c    -------------------------
c    w(mdw,*),mcon,mrows,ncols
c    -------------------------
c     the array w contains the (possibly null) matrix [c:*] followed by
c     [e:f].  this must be placed in w as follows:
c          [c  :  *]
c     w  = [       ]
c          [e  :  f]
c     the (*) after c indicates that this data can be undefined. the
c     matrix [e:f] has mrows rows and ncols+1 columns. the matrix c is
c     placed in the first mcon rows of w(*,*) while [e:f]
c     follows in rows mcon+1 through mcon+mrows of w(*,*). the vector f
c     is placed in rows mcon+1 through mcon+mrows, column ncols+1. the
c     values of mdw and ncols must be positive; the value of mcon must
c     be nonnegative. an exception to this occurs when using option 1
c     for accumulation of blocks of equations. in that case mrows is an
c     output variable only, and the matrix data for [e:f] is placed in
c     w(*,*), one block of rows at a time. see iopt(*) contents, option
c     number 1, for further details. the row dimension, mdw, of the
c     array w(*,*) must satisfy the inequality:
c
c     if using option 1,
c                     mdw .ge. mcon + max(max. number of
c                     rows accumulated, ncols) + 1.
c     if using option 8,
c                     mdw .ge. mcon + mrows.
c     else
c                     mdw .ge. mcon + max(mrows, ncols).
c
c     other values are errors, but this is checked only when using
c     option=2.  the value of mrows is an output parameter when
c     using option number 1 for accumulating large blocks of least
c     squares equations before solving the problem.
c     see iopt(*) contents for details about option 1.
c
c    ------------------
c    bl(*),bu(*),ind(*)
c    ------------------
c     these arrays contain the information about the bounds that the
c     solution values are to satisfy. the value of ind(j) tells the
c     type of bound and bl(j) and bu(j) give the explicit values for
c     the respective upper and lower bounds on the unknowns x and y.
c     the first nvars entries of ind(*), bl(*) and bu(*) specify
c     bounds on x; the next mcon entries specify bounds on y.
c
c    1.    for ind(j)=1, require x(j) .ge. bl(j);
c          if j.gt.ncols,        y(j-ncols) .ge. bl(j).
c          (the value of bu(j) is not used.)
c    2.    for ind(j)=2, require x(j) .le. bu(j);
c          if j.gt.ncols,        y(j-ncols) .le. bu(j).
c          (the value of bl(j) is not used.)
c    3.    for ind(j)=3, require x(j) .ge. bl(j) and
c                                x(j) .le. bu(j);
c          if j.gt.ncols,        y(j-ncols) .ge. bl(j) and
c                                y(j-ncols) .le. bu(j).
c          (to impose equality constraints have bl(j)=bu(j)=
c          constraining value.)
c    4.    for ind(j)=4, no bounds on x(j) or y(j-ncols) are required.
c          (the values of bl(j) and bu(j) are not used.)
c
c     values other than 1,2,3 or 4 for ind(j) are errors. in the case
c     ind(j)=3 (upper and lower bounds) the condition bl(j) .gt. bu(j)
c     is  an  error.   the values bl(j), bu(j), j .gt. ncols, will be
c     changed.  significant changes mean that the constraints are
c     infeasible.  (users must make this decision themselves.)
c     the new values for bl(j), bu(j), j .gt. ncols, define a
c     region such that the perturbed problem is feasible.  if users
c     know that their problem is feasible, this step can be skipped
c     by using option number 8 described below.
c     see iopt(*) description.
c
c
c    -------
c    iopt(*)
c    -------
c     this is the array where the user can specify nonstandard options
c     for dbocls( ). most of the time this feature can be ignored by
c     setting the input value iopt(1)=99. occasionally users may have
c     needs that require use of the following subprogram options. for
c     details about how to use the options see below: iopt(*) contents.
c
c     option number   brief statement of purpose
c     ------ ------   ----- --------- -- -------
c           1         return to user for accumulation of blocks
c                     of least squares equations.  the values
c                     of iopt(*) are changed with this option.
c                     the changes are updates to pointers for
c                     placing the rows of equations into position
c                     for processing.
c           2         check lengths of all arrays used in the
c                     subprogram.
c           3         column scaling of the data matrix, [c].
c                                                        [e]
c           4         user provides column scaling for matrix [c].
c                                                             [e]
c           5         provide option array to the low-level
c                     subprogram sbols( ).
c           6         provide option array to the low-level
c                     subprogram sbolsm( ).
c           7         move the iopt(*) processing pointer.
c           8         do not preprocess the constraints to
c                     resolve infeasibilities.
c           9         do not pretriangularize the least squares matrix.
c          99         no more options to change.
c
c    ----
c    x(*)
c    ----
c     this array is used to pass data associated with options 4,5 and
c     6. ignore this parameter (on input) if no options are used.
c     otherwise see below: iopt(*) contents.
c
c
c    output
c    ------
c
c    -----------------
c    x(*),rnormc,rnorm
c    -----------------
c     the array x(*) contains a solution (if mode .ge.0 or .eq.-22) for
c     the constrained least squares problem. the value rnormc is the
c     minimum residual vector length for the constraints c*x - y = 0.
c     the value rnorm is the minimum residual vector length for the
c     least squares equations. normally rnormc=0, but in the case of
c     inconsistent constraints this value will be nonzero.
c     the values of x are returned in the first nvars entries of x(*).
c     the values of y are returned in the last mcon entries of x(*).
c
c    ----
c    mode
c    ----
c     the sign of mode determines whether the subprogram has completed
c     normally, or encountered an error condition or abnormal status. a
c     value of mode .ge. 0 signifies that the subprogram has completed
c     normally. the value of mode (.ge. 0) is the number of variables
c     in an active status: not at a bound nor at the value zero, for
c     the case of free variables. a negative value of mode will be one
c     of the cases (-57)-(-41), (-37)-(-22), (-19)-(-2). values .lt. -1
c     correspond to an abnormal completion of the subprogram. these
c     error messages are in groups for the subprograms dbocls(),
c     sbolsm(), and sbols().  an approximate solution will be returned
c     to the user only when max. iterations is reached, mode=-22.
c
c    -----------
c    rw(*),iw(*)
c    -----------
c     these are working arrays.  (normally the user can ignore the
c     contents of these arrays.)
c
c    iopt(*) contents
c    ------- --------
c     the option array allows a user to modify some internal variables
c     in the subprogram without recompiling the source code. a central
c     goal of the initial software design was to do a good job for most
c     people. thus the use of options will be restricted to a select
c     group of users. the processing of the option array proceeds as
c     follows: a pointer, here called lp, is initially set to the value
c     1. at the pointer position the option number is extracted and
c     used for locating other information that allows for options to be
c     changed. the portion of the array iopt(*) that is used for each
c     option is fixed; the user and the subprogram both know how many
c     locations are needed for each option. the value of lp is updated
c     for each option based on the amount of storage in iopt(*) that is
c     required. a great deal of error checking is done by the
c     subprogram on the contents of the option array. nevertheless it
c     is still possible to give the subprogram optional input that is
c     meaningless. for example option 4 uses the locations
c     x(ncols+ioff),...,x(ncols+ioff+ncols-1) for passing scaling data.
c     the user must manage the allocation of these locations.
c
c   1
c   -
c     this option allows the user to solve problems with a large number
c     of rows compared to the number of variables. the idea is that the
c     subprogram returns to the user (perhaps many times) and receives
c     new least squares equations from the calling program unit.
c     eventually the user signals "that's all" and a solution is then
c     computed. the value of mrows is an output variable when this
c     option is used. its value is always in the range 0 .le. mrows
c     .le. ncols+1. it is the number of rows after the
c     triangularization of the entire set of equations. if lp is the
c     processing pointer for iopt(*), the usage for the sequential
c     processing of blocks of equations is
c
c
c        iopt(lp)=1
c         move block of equations to w(*,*) starting at
c         the first row of w(*,*).
c        iopt(lp+3)=# of rows in the block; user defined
c
c     the user now calls dbocls( ) in a loop. the value of iopt(lp+1)
c     directs the user's action. the value of iopt(lp+2) points to
c     where the subsequent rows are to be placed in w(*,*). both of
c     these values are first defined in the subprogram. the user
c     changes the value of iopt(lp+1) (to 2) as a signal that all of
c     the rows have been processed.
c
c
c      .<loop
c      . call dbocls( )
c      . if(iopt(lp+1) .eq. 1) then
c      .    iopt(lp+3)=# of rows in the new block; user defined
c      .    place new block of iopt(lp+3) rows in
c      .    w(*,*) starting at row mcon + iopt(lp+2).
c      .
c      .    if( this is the last block of equations ) then
c      .       iopt(lp+1)=2
c      .<------cycle loop
c      .    else if (iopt(lp+1) .eq. 2) then
c      <-------exit loop solution computed if mode .ge. 0
c      . else
c      . error condition; should not happen.
c      .<end loop
c
c     use of this option adds 4 to the required length of iopt(*).
c
c   2
c   -
c     this option is useful for checking the lengths of all arrays used
c     by dbocls( ) against their actual requirements for this problem.
c     the idea is simple: the user's program unit passes the declared
c     dimension information of the arrays. these values are compared
c     against the problem-dependent needs within the subprogram. if any
c     of the dimensions are too small an error message is printed and a
c     negative value of mode is returned, -41 to -47. the printed error
c     message tells how long the dimension should be. if lp is the
c     processing pointer for iopt(*),
c
c        iopt(lp)=2
c        iopt(lp+1)=row dimension of w(*,*)
c        iopt(lp+2)=col. dimension of w(*,*)
c        iopt(lp+3)=dimensions of bl(*),bu(*),ind(*)
c        iopt(lp+4)=dimension of x(*)
c        iopt(lp+5)=dimension of rw(*)
c        iopt(lp+6)=dimension of iw(*)
c        iopt(lp+7)=dimension of iopt(*)
c         .
c        call dbocls( )
c
c     use of this option adds 8 to the required length of iopt(*).
c
c   3
c   -
c     this option can change the type of scaling for the data matrix.
c     nominally each nonzero column of the matrix is scaled so that the
c     magnitude of its largest entry is equal to the value one. if lp
c     is the processing pointer for iopt(*),
c
c        iopt(lp)=3
c        iopt(lp+1)=1,2 or 3
c            1= nominal scaling as noted;
c            2= each nonzero column scaled to have length one;
c            3= identity scaling; scaling effectively suppressed.
c         .
c        call dbocls( )
c
c     use of this option adds 2 to the required length of iopt(*).
c
c   4
c   -
c     this options allows the user to provide arbitrary (positive)
c     column scaling for the matrix. if lp is the processing pointer
c     for iopt(*),
c
c        iopt(lp)=4
c        iopt(lp+1)=ioff
c        x(ncols+ioff),...,x(ncols+ioff+ncols-1)
c        = positive scale factors for cols. of e.
c         .
c        call dbocls( )
c
c     use of this option adds 2 to the required length of iopt(*)
c     and ncols to the required length of x(*).
c
c   5
c   -
c     this option allows the user to provide an option array to the
c     low-level subprogram sbols( ). if lp is the processing pointer
c     for iopt(*),
c
c        iopt(lp)=5
c        iopt(lp+1)= position in iopt(*) where option array
c                    data for sbols( ) begins.
c         .
c        call dbocls( )
c
c     use of this option adds 2 to the required length of iopt(*).
c
c   6
c   -
c     this option allows the user to provide an option array to the
c     low-level subprogram sbolsm( ). if lp is the processing pointer
c     for iopt(*),
c
c        iopt(lp)=6
c        iopt(lp+1)= position in iopt(*) where option array
c                    data for sbolsm( ) begins.
c         .
c        call dbocls( )
c
c     use of this option adds 2 to the required length of iopt(*).
c
c   7
c   -
c     move the processing pointer (either forward or backward) to the
c     location iopt(lp+1). the processing pointer moves to locations
c     lp+2 if option number 7 is used with the value -7.  for
c     example to skip over locations 3,...,ncols+2,
c
c       iopt(1)=7
c       iopt(2)=ncols+3
c       (iopt(i), i=3,...,ncols+2 are not defined here.)
c       iopt(ncols+3)=99
c       call dbocls( )
c
c     caution: misuse of this option can yield some very hard-to-find
c     bugs. use it with care. it is intended to be used for passing
c     option arrays to other subprograms.
c
c   8
c   -
c     this option allows the user to suppress the algorithmic feature
c     of dbocls( ) that processes the constraint equations c*x = y and
c     resolves infeasibilities. the steps normally done are to solve
c     c*x - y = 0 in a least squares sense using the stated bounds on
c     both x and y. then the "reachable" vector y = c*x is computed
c     using the solution x obtained. finally the stated bounds for y are
c     enlarged to include c*x. to suppress the feature:
c
c
c       iopt(lp)=8
c         .
c       call dbocls( )
c
c     use of this option adds 1 to the required length of iopt(*).
c
c   9
c   -
c     this option allows the user to suppress the pretriangularizing
c     step of the least squares matrix that is done within dbocls( ).
c     this is primarily a means of enhancing the subprogram efficiency
c     and has little effect on accuracy. to suppress the step, set:
c
c       iopt(lp)=9
c         .
c       call dbocls( )
c
c     use of this option adds 1 to the required length of iopt(*).
c
c   99
c   --
c     there are no more options to change.
c
c     only option numbers -99, -9,-8,...,-1, 1,2,...,9, and 99 are
c     permitted. other values are errors. options -99,-1,...,-9 mean
c     that the respective options 99,1,...,9 are left at their default
c     values. an example is the option to suppress the preprocessing of
c     constraints:
c
c       iopt(1)=-8 option is recognized but not changed
c       iopt(2)=99
c       call dbocls( )
c
c    error messages for dbocls()
c    ----- -------- --- --------
c
c warning in...
c dbocls(). the row dimension of w(,)=(i1) must be .ge. the number
c of effective rows=(i2).
c           in above message, i1=         1
c           in above message, i2=         2
c error number =        41
c
c warning in...
c dbocls(). the column dimension of w(,)=(i1) must be .ge. ncols+
c mcon+1=(i2).
c           in above message, i1=         2
c           in above message, i2=         3
c error number =        42
c
c warning in...
c dbocls(). the dimensions of the arrays bl(),bu(), and ind()=(i1)
c must be .ge. ncols+mcon=(i2).
c           in above message, i1=         1
c           in above message, i2=         2
c error number =        43
c
c warning in...
c dbocls(). the dimension of x()=(i1) must be
c .ge. the reqd.length=(i2).
c           in above message, i1=         1
c           in above message, i2=         2
c error number =        44
c
c warning in...
c dbocls(). the .
c dbocls() the dimension of iw()=(i1) must be .ge. 2*ncols+2*mcon=(i2).
c           in above message, i1=         1
c           in above message, i2=         4
c error number =        46
c
c warning in...
c dbocls(). the dimension of iopt()=(i1) must be .ge. the reqd.
c len.=(i2).
c           in above message, i1=        16
c           in above message, i2=        18
c error number =        47
c
c warning in...
c dbocls(). iscale option=(i1) must be 1-3.
c           in above message, i1=         0
c error number =        48
c
c warning in...
c dbocls(). offset past x(ncols) (i1) for user-provided column scaling
c must be positive.
c           in above message, i1=         0
c error number =        49
c
c warning in...
c dbocls(). each provided col. scale factor must be positive.
c  component (i1) now = (r1).
c           in above message, i1=         1
c           in above message, r1=    0.
c error number =        50
c
c warning in...
c dbocls(). the option number=(i1) is not defined.
c           in above message, i1=      1001
c error number =        51
c
c warning in...
c dbocls(). no. of rows=(i1) must be .ge. 0 .and. .le. mdw-mcon=(i2).
c           in above message, i1=         2
c           in above message, i2=         1
c error number =        52
c
c warning in...
c dbocls(). mdw=(i1) must be positive.
c           in above message, i1=         0
c error number =        53
c
c warning in...
c dbocls(). mcon=(i1) must be nonnegative.
c           in above message, i1=        -1
c error number =        54
c
c warning in...
c dbocls(). ncols=(i1) the no. of variables must be positive.
c           in above message, i1=         0
c error number =        55
c
c warning in...
c dbocls(). for j=(i1), ind(j)=(i2) must be 1-4.
c           in above message, i1=         1
c           in above message, i2=         0
c error number =        56
c
c warning in...
c dbocls(). for j=(i1), bound bl(j)=(r1) is .gt. bu(j)=(r2).
c           in above message, i1=         1
c           in above message, r1=     .1000000000e+01
c           in above message, r2=    0.
c error number =        57
c           linear constraints, snla rept. sand82-1517, aug., (1982).
c
c***references  r. j. hanson, linear least squares with bounds and
c                 linear constraints, report sand82-1517, sandia
c                 laboratories, august 1982.
c***routines called  d1mach, dasum, dbols, dcopy, ddot, dnrm2, dscal,
c                    xermsg
c***revision history  (yymmdd)
c   821220  date written
c   891006  cosmetic changes to prologue.  (wrb)
c   891006  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900510  convert xerrwv calls to xermsg calls.  (rwc)
c   910819  added variable m for mout+mcon in reference to dbols.  (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  dbocls
c     revised 850604-0900
c     revised yymmdd-hhmm
c
c    purpose
c    -------
c     this is the main subprogram that solves the least squares
c     problem consisting of linear constraints
c
c              c*x = y
c
c     and least squares equations
c
c              e*x = f
c
c     in this formulation the vectors x and y are both unknowns.
c     further, x and y may both have user-specified bounds on each
c     component.  the user must have dimension statements of the
c     form
c
c     dimension w(mdw,ncols+mcon+1), bl(ncols+mcon),bu(ncols+mcon),
c               x(2*(ncols+mcon)+2+nx), rw(6*ncols+5*mcon)
c
c     integer ind(ncols+mcon), iopt(16+ni), iw(2*(ncols+mcon))
c
c     to change this subprogram from single to double precision begin
c     editing at the card 'c++'.
c     change this subprogram to dbocls and the strings
c     /sdot/ to /ddot/, /snrm2/ to /dnrm2/, /srelpr/ to /drelpr/,
c     /r1mach/ to /d1mach/, /e0/ to /d0/, /scopy/ to /dcopy/,
c     /sscal/ to /dscal/, /sasum/ to /dasum/, /sbols/ to /dbols/,
c     /real            / to /double precision/.
c ++
      double precision w(mdw,*),bl(*),bu(*),x(*),rw(*)
      double precision anorm, cnorm, one, rnorm, rnormc, drelpr
      double precision t, t1, t2, ddot, dnrm2, wt, zero
      double precision dasum, d1mach
c     this variable remains typed real.
      integer ind(*),iopt(*),iw(*),jopt(05)
      logical checkl,filter,accum,pretri
      character*8 xern1, xern2
      character*16 xern3, xern4
      save igo,accum,checkl
      data igo/0/
c***first executable statement  dbocls
      nerr = 0
      mode = 0
      if (igo.eq.0) then
c     do(check validity of input data)
c     procedure(check validity of input data)
c
c     see that mdw is .gt.0. gross check only.
          if (mdw.le.0) then
              write (xern1, '(i8)') mdw
              call xermsg ('slatec', 'dbocls', 'mdw = ' // xern1 //
     *           ' must be positive.', 53, 1)
c     do(return to user program unit)
              go to 260
          endif
c
c     see that number of constraints is nonnegative.
          if (mcon.lt.0) then
              write (xern1, '(i8)') mcon
              call xermsg ('slatec', 'dbocls', 'mcon = ' // xern1 //
     *           ' must be non-negative', 54, 1)
c     do(return to user program unit)
              go to 260
          endif
c
c     see that number of unknowns is positive.
          if (ncols.le.0) then
              write (xern1, '(i8)') ncols
              call xermsg ('slatec', 'dbocls', 'ncols = ' // xern1 //
     *           ' the no. of variables, must be positive.', 55, 1)
c     do(return to user program unit)
              go to 260
          endif
c
c     see that constraint indicators are all well-defined.
          do 10 j = 1,ncols + mcon
              if (ind(j).lt.1 .or. ind(j).gt.4) then
                  write (xern1, '(i8)') j
                  write (xern2, '(i8)') ind(j)
                  call xermsg ('slatec', 'dbocls', 'ind(' // xern1 //
     *               ') = ' // xern2 // ' must be 1-4.', 56, 1)
c     do(return to user program unit)
                  go to 260
              endif
   10     continue
c
c     see that bounds are consistent.
          do 20 j = 1,ncols + mcon
              if (ind(j).eq.3) then
                  if (bl(j).gt.bu(j)) then
                     write (xern1, '(i8)') j
                     write (xern3, '(1pe15.6)') bl(j)
                     write (xern4, '(1pe15.6)') bu(j)
                     call xermsg ('slatec', 'dbocls', 'bound bl(' //
     *                  xern1 // ') = ' // xern3 // ' is .gt. bu(' //
     *                  xern1 // ') = ' // xern4, 57, 1)
c     do(return to user program unit)
                      go to 260
                  endif
              endif
   20     continue
c     end procedure
c     do(process option array)
c     procedure(process option array)
          zero = 0.d0
          one = 1.d0
          drelpr = d1mach(4)
          checkl = .false.
          filter = .true.
          lenx = 2* (ncols+mcon) + 2
          iscale = 1
          igo = 1
          accum = .false.
          pretri = .true.
          lopt = 0
          mopt = 0
          lp = 0
          lds = 0
c     do forever
   30     continue
          lp = lp + lds
          ip = iopt(lp+1)
          jp = abs(ip)
c
c     test for no more options to change.
          if (ip.eq.99) then
              if (lopt.eq.0) lopt = - (lp+2)
              if (mopt.eq.0) mopt = - (abs(lopt)+7)
              if (lopt.lt.0) then
                  lbou = abs(lopt)
              else
                  lbou = lopt - 15
              endif
c
c     send col. scaling to dbols().
              iopt(lbou) = 4
              iopt(lbou+1) = 1
c
c     pass an option array for dbolsm().
              iopt(lbou+2) = 5
c
c     loc. of option array for dbolsm( ).
              iopt(lbou+3) = 8
c
c     skip to start of user-given option array for dbols().
              iopt(lbou+4) = 6
              iopt(lbou+6) = 99
              if (lopt.gt.0) then
                  iopt(lbou+5) = lopt - lbou + 1
              else
                  iopt(lbou+4) = -iopt(lbou+4)
              endif
              if (mopt.lt.0) then
                  lboum = abs(mopt)
              else
                  lboum = mopt - 8
              endif
c
c     change pretriangularization factor in dbolsm().
              iopt(lboum) = 5
              iopt(lboum+1) = ncols + mcon + 1
c
c     pass weight to dbolsm() for rank test.
              iopt(lboum+2) = 6
              iopt(lboum+3) = ncols + mcon + 2
              iopt(lboum+4) = mcon
c
c     skip to user-given option array for dbolsm( ).
              iopt(lboum+5) = 1
              iopt(lboum+7) = 99
              if (mopt.gt.0) then
                  iopt(lboum+6) = mopt - lboum + 1
              else
                  iopt(lboum+5) = -iopt(lboum+5)
              endif
c     exit forever
              go to 50
          else if (jp.eq.99) then
              lds = 1
c     cycle forever
              go to 50
          else if (jp.eq.1) then
              if (ip.gt.0) then
c
c     set up direction flag location, row stacking pointer
c     location, and location for number of new rows.
                  locacc = lp + 2
c
c                  iopt(locacc-1)=option number for seq. accumulation.
c     contents..   iopt(locacc  )=user direction flag, 1 or 2.
c                  iopt(locacc+1)=row stacking pointer.
c                  iopt(locacc+2)=number of new rows to process.
c     user action with this option..
c      (set up option data for seq. accumulation in iopt(*).)
c      (move block of equations into w(*,*)  starting at first
c       row of w(*,*) below the rows for the constraint matrix c.
c       set iopt(locacc+2)=no. of least squares equations in block.
c              loop
c              call dbocls()
c
c                  if(iopt(locacc) .eq. 1) then
c                      stack equas. into w(*,*), starting at
c                      row iopt(locacc+1).
c                       into w(*,*).
c                       set iopt(locacc+2)=no. of equas.
c                      if last block of equas., set iopt(locacc)=2.
c                  else if iopt(locacc) .eq. 2) then
c                      (process is over. exit loop.)
c                  else
c                      (error condition. should not happen.)
c                  end if
c              end loop
                  iopt(locacc+1) = mcon + 1
                  accum = .true.
                  iopt(locacc) = igo
              endif
              lds = 4
c     cycle forever
              go to 30
          else if (jp.eq.2) then
              if (ip.gt.0) then
c
c     get actual lengths of arrays for checking against needs.
                  locdim = lp + 2
c
c     lmdw.ge.mcon+max(mout,ncols), if mcon.gt.0 .and filter
c     lmdw.ge.mcon+mout, otherwise
c
c     lndw.ge.ncols+mcon+1
c     llb .ge.ncols+mcon
c     llx .ge.2*(ncols+mcon)+2+extra reqd. in options.
c     llrw.ge.6*ncols+5*mcon
c     lliw.ge.2*(ncols+mcon)
c     liop.ge. amount reqd. for option array.
                  lmdw = iopt(locdim)
                  lndw = iopt(locdim+1)
                  llb = iopt(locdim+2)
                  llx = iopt(locdim+3)
                  llrw = iopt(locdim+4)
                  lliw = iopt(locdim+5)
                  liopt = iopt(locdim+6)
                  checkl = .true.
              endif
              lds = 8
c     cycle forever
              go to 30
c
c     option to modify the column scaling.
          else if (jp.eq.3) then
              if (ip.gt.0) then
                  iscale = iopt(lp+2)
c
c     see that iscale is 1 thru 3.
                  if (iscale.lt.1 .or. iscale.gt.3) then
                      write (xern1, '(i8)') iscale
                      call xermsg ('slatec', 'dbocls',
     *                   'iscale option = ' // xern1 // ' must be 1-3',
     *                   48, 1)
c     do(return to user program unit)
                      go to 260
                  endif
              endif
              lds = 2
c     cycle forever
              go to 30
c
c     in this option the user has provided scaling.  the
c     scale factors for the columns begin in x(ncols+iopt(lp+2)).
          else if (jp.eq.4) then
              if (ip.gt.0) then
                  iscale = 4
                  if (iopt(lp+2).le.0) then
                      write (xern1, '(i8)') iopt(lp+2)
                      call xermsg ('slatec', 'dbocls',
     *                   'offset past x(ncols) (' // xern1 //
     *           ') for user-provided column scaling must be positive.',
     *                   49, 1)
c     do(return to user program unit)
                      go to 260
                  endif
                  call dcopy(ncols,x(ncols+iopt(lp+2)),1,rw,1)
                  lenx = lenx + ncols
                  do 40 j = 1,ncols
                      if (rw(j).le.zero) then
                          write (xern1, '(i8)') j
                          write (xern3, '(1pe15.6)') rw(j)
                          call xermsg ('slatec', 'dbocls',
     *            'each provided column scale factor must be ' //
     *            'positive.$$component ' // xern1 // ' now = ' //
     *            xern3, 50, 1)
c     do(return to user program unit)
                          go to 260
                      endif
   40             continue
              endif
              lds = 2
c     cycle forever
              go to 30
c
c     in this option an option array is provided to dbols().
          else if (jp.eq.5) then
              if (ip.gt.0) then
                  lopt = iopt(lp+2)
              endif
              lds = 2
c     cycle forever
              go to 30
c
c     in this option an option array is provided to dbolsm().
          else if (jp.eq.6) then
              if (ip.gt.0) then
                  mopt = iopt(lp+2)
              endif
              lds = 2
c     cycle forever
              go to 30
c
c     this option uses the next loc of iopt(*) as a
c     pointer value to skip to next.
          else if (jp.eq.7) then
              if (ip.gt.0) then
                  lp = iopt(lp+2) - 1
                  lds = 0
              else
                  lds = 2
              endif
c     cycle forever
              go to 30
c
c     this option avoids the constraint resolving phase for
c     the linear constraints c*x=y.
          else if (jp.eq.8) then
              filter = .not. (ip.gt.0)
              lds = 1
c     cycle forever
              go to 30
c
c     this option suppresses pre-triangularization of the least
c     squares equations.
          else if (jp.eq.9) then
              pretri = .not. (ip.gt.0)
              lds = 1
c     cycle forever
              go to 30
c
c     no valid option number was noted. this is an error condition.
          else
              write (xern1, '(i8)') jp
              call xermsg ('slatec', 'dbocls', 'option number = ' //
     *           xern1 // ' is not defined.', 51, 1)
c     do(return to user program unit)
              go to 260
          endif
c     end forever
c     end procedure
   50     continue
          if (checkl) then
c     do(check lengths of arrays)
c     procedure(check lengths of arrays)
c
c     this feature allows the user to make sure that the
c     arrays are long enough for the intended problem size and use.
           if(filter .and. .not.accum) then
                mdwl=mcon+max(mrows,ncols)
           else
                mdwl=mcon+ncols+1
           endif
              if (lmdw.lt.mdwl) then
                  write (xern1, '(i8)') lmdw
                  write (xern2, '(i8)') mdwl
                  call xermsg ('slatec', 'dbocls',
     *               'the row dimension of w(,) = ' // xern1 //
     *               ' must be .ge. the number of effective rows = ' //
     *               xern2, 41, 1)
c     do(return to user program unit)
                  go to 260
              endif
              if (lndw.lt.ncols+mcon+1) then
                  write (xern1, '(i8)') lndw
                  write (xern2, '(i8)') ncols+mcon+1
                  call xermsg ('slatec', 'dbocls',
     *               'the column dimension of w(,) = ' // xern1 //
     *               ' must be .ge. ncols+mcon+1 = ' // xern2, 42, 1)
c     do(return to user program unit)
                  go to 260
              endif
              if (llb.lt.ncols+mcon) then
                  write (xern1, '(i8)') llb
                  write (xern2, '(i8)') ncols+mcon
                  call xermsg ('slatec', 'dbocls',
     *           'the dimensions of the arrays bs(), bu(), and ind() = '
     *               // xern1 // ' must be .ge. ncols+mcon = ' // xern2,
     *               43, 1)
c     do(return to user program unit)
                  go to 260
              endif
              if (llx.lt.lenx) then
                  write (xern1, '(i8)') llx
                  write (xern2, '(i8)') lenx
                  call xermsg ('slatec', 'dbocls',
     *              'the dimension of x() = ' // xern1 //
     *              ' must be .ge. the required length = ' // xern2,
     *              44, 1)
c     do(return to user program unit)
                  go to 260
              endif
              if (llrw.lt.6*ncols+5*mcon) then
                  write (xern1, '(i8)') llrw
                  write (xern2, '(i8)') 6*ncols+5*mcon
                  call xermsg ('slatec', 'dbocls',
     *               'the dimension of rw() = ' // xern1 //
     *               ' must be .ge. 6*ncols+5*mcon = ' // xern2, 45, 1)
c     do(return to user program unit)
                  go to 260
              endif
              if (lliw.lt.2*ncols+2*mcon) then
                  write (xern1, '(i8)') lliw
                  write (xern2, '(i8)') 2*ncols+2*mcon
                  call xermsg ('slatec', 'dbocls',
     *               'the dimension of iw() = ' // xern1 //
     *               ' must be .ge. 2*ncols+2*mcon = ' // xern2, 46, 1)
c     do(return to user program unit)
                  go to 260
              endif
              if (liopt.lt.lp+17) then
                  write (xern1, '(i8)') liopt
                  write (xern2, '(i8)') lp+17
                  call xermsg ('slatec', 'dbocls',
     *               'the dimension of iopt() = ' // xern1 //
     *               ' must be .ge. the required len = ' // xern2, 47,1)
c     do(return to user program unit)
                  go to 260
              endif
c     end procedure
          endif
      endif
c
c     optionally go back to the user for accumulation of least squares
c     equations and directions for processing these equations.
c     do(accumulate least squares equations)
c     procedure(accumulate least squares equations)
      if (accum) then
          mrows = iopt(locacc+1) - 1 - mcon
          inrows = iopt(locacc+2)
          mnew = mrows + inrows
          if (mnew.lt.0 .or. mnew+mcon.gt.mdw) then
              write (xern1, '(i8)') mnew
              write (xern2, '(i8)') mdw-mcon
              call xermsg ('slatec', 'dbocls', 'no. of rows = ' //
     *           xern1 //  ' must be .ge. 0 .and. .le. mdw-mcon = ' //
     *           xern2, 52, 1)
c    (return to user program unit)
              go to 260
          endif
      endif
c
c     use the software of dbols( ) for the triangularization of the
c     least squares matrix.  this may involve a systaltic interchange
c     of processing pointers between the calling and called (dbols())
c     program units.
      jopt(01) = 1
      jopt(02) = 2
      jopt(04) = mrows
      jopt(05) = 99
      irw = ncols + 1
      iiw = 1
      if (accum .or. pretri) then
          call dbols(w(mcon+1,1),mdw,mout,ncols,bl,bu,ind,jopt,x,rnorm,
     *               mode,rw(irw),iw(iiw))
      else
          mout = mrows
      endif
      if (accum) then
          accum = iopt(locacc) .eq. 1
          iopt(locacc+1) = jopt(03) + mcon
          mrows = min(ncols+1,mnew)
      endif
c     end procedure
      if (accum) return
c     do(solve constrained and bounded least squares problem)
c     procedure(solve constrained and bounded least squares problem)
c
c     move right hand side of least squares equations.
      call dcopy(mout,w(mcon+1,ncols+1),1,w(mcon+1,ncols+mcon+1),1)
      if (mcon.gt.0 .and. filter) then
c
c     project the linear constraints into a reachable set.
          do 60 i = 1,mcon
              call dcopy(ncols,w(i,1),mdw,w(mcon+1,ncols+i),1)
   60     continue
c
c      place (-)identity matrix after constraint data.
          do 70 j = ncols + 1,ncols + mcon + 1
              w(1,j) = zero
              call dcopy(mcon,w(1,j),0,w(1,j),1)
   70     continue
          w(1,ncols+1) = -one
          call dcopy(mcon,w(1,ncols+1),0,w(1,ncols+1),mdw+1)
c
c     obtain a 'feasible point' for the linear constraints.
          jopt(01) = 99
          irw = ncols + 1
          iiw = 1
          call dbols(w,mdw,mcon,ncols+mcon,bl,bu,ind,jopt,x,rnormc,
     *               modec,rw(irw),iw(iiw))
c
c     enlarge the bounds set, if required, to include points that
c     can be reached.
          do 130 j = ncols + 1,ncols + mcon
              icase = ind(j)
              if (icase.lt.4) then
                  t = ddot(ncols,w(mcon+1,j),1,x,1)
              endif
              go to (80,90,100,110),icase
              go to 120
c     case 1
   80         bl(j) = min(t,bl(j))
              go to 120
c     case 2
   90         bu(j) = max(t,bu(j))
              go to 120
c     case 3
  100         bl(j) = min(t,bl(j))
              bu(j) = max(t,bu(j))
              go to 120
c     case 4
  110         continue
  120         continue
  130     continue
c
c     move constraint data back to the original area.
          do 140 j = ncols + 1,ncols + mcon
              call dcopy(ncols,w(mcon+1,j),1,w(j-ncols,1),mdw)
  140     continue
      endif
      if (mcon.gt.0) then
          do 150 j = ncols + 1,ncols + mcon
              w(mcon+1,j) = zero
              call dcopy(mout,w(mcon+1,j),0,w(mcon+1,j),1)
  150     continue
c
c     put in (-)identity matrix (possibly) once again.
          do 160 j = ncols + 1,ncols + mcon + 1
              w(1,j) = zero
              call dcopy(mcon,w(1,j),0,w(1,j),1)
  160     continue
          w(1,ncols+1) = -one
          call dcopy(mcon,w(1,ncols+1),0,w(1,ncols+1),mdw+1)
      endif
c
c     compute nominal column scaling for the unweighted matrix.
      cnorm = zero
      anorm = zero
      do 170 j = 1,ncols
          t1 = dasum(mcon,w(1,j),1)
          t2 = dasum(mout,w(mcon+1,1),1)
          t = t1 + t2
          if (t.eq.zero) t = one
          cnorm = max(cnorm,t1)
          anorm = max(anorm,t2)
          x(ncols+mcon+j) = one/t
  170 continue
      go to (180,190,210,220),iscale
      go to 230
c     case 1
  180 continue
      go to 230
c     case 2
c
c     scale cols. (before weighting) to have length one.
  190 do 200 j = 1,ncols
          t = dnrm2(mcon+mout,w(1,j),1)
          if (t.eq.zero) t = one
          x(ncols+mcon+j) = one/t
  200 continue
      go to 230
c     case 3
c
c     suppress scaling (use unit matrix).
  210 x(ncols+mcon+1) = one
      call dcopy(ncols,x(ncols+mcon+1),0,x(ncols+mcon+1),1)
      go to 230
c     case 4
c
c     the user has provided scaling.
  220 call dcopy(ncols,rw,1,x(ncols+mcon+1),1)
  230 continue
      do 240 j = ncols + 1,ncols + mcon
          x(ncols+mcon+j) = one
  240 continue
c
c     weight the least squares equations.
      wt = drelpr
      if (anorm.gt.zero) wt = wt/anorm
      if (cnorm.gt.zero) wt = wt*cnorm
      do 250 i = 1,mout
          call dscal(ncols,wt,w(i+mcon,1),mdw)
  250 continue
      call dscal(mout,wt,w(mcon+1,mcon+ncols+1),1)
      lrw = 1
      liw = 1
c
c     set the new triangularization factor.
      x(2* (ncols+mcon)+1) = zero
c
c     set the weight to use in components .gt. mcon,
c     when making linear independence test.
      x(2* (ncols+mcon)+2) = one/wt
      m = mout+mcon
      call dbols(w,mdw,m,ncols+mcon,bl,bu,ind,iopt(lbou),x,
     *           rnorm,mode,rw(lrw),iw(liw))
      rnorm = rnorm/wt
c     end procedure
c     procedure(return to user program unit)
  260 if (mode.ge.0) mode = -nerr
      igo = 0
      return
c     end program
      end
