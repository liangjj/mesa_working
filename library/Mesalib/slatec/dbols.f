*deck dbols
      subroutine dbols (w, mdw, mrows, ncols, bl, bu, ind, iopt, x,
     +   rnorm, mode, rw, iw)
c***begin prologue  dbols
c***purpose  solve the problem
c                 e*x = f (in the least  squares  sense)
c            with bounds on selected x values.
c***library   slatec
c***category  k1a2a, g2e, g2h1, g2h2
c***type      double precision (sbols-s, dbols-d)
c***keywords  bounds, constraints, inequality, least squares, linear
c***author  hanson, r. j., (snla)
c***description
c
c   **** all input and output real variables are double precision ****
c
c     the user must have dimension statements of the form:
c
c       dimension w(mdw,ncols+1), bl(ncols), bu(ncols),
c      * x(ncols+nx), rw(5*ncols)
c       integer ind(ncols), iopt(1+ni), iw(2*ncols)
c
c     (here nx=number of extra locations required for option 4; nx=0
c     for no options; nx=ncols if this option is in use. here ni=number
c     of extra locations required for options 1-6; ni=0 for no
c     options.)
c
c   input
c   -----
c
c    --------------------
c    w(mdw,*),mrows,ncols
c    --------------------
c     the array w(*,*) contains the matrix [e:f] on entry. the matrix
c     [e:f] has mrows rows and ncols+1 columns. this data is placed in
c     the array w(*,*) with e occupying the first ncols columns and the
c     right side vector f in column ncols+1. the row dimension, mdw, of
c     the array w(*,*) must satisfy the inequality mdw .ge. mrows.
c     other values of mdw are errors. the values of mrows and ncols
c     must be positive. other values are errors. there is an exception
c     to this when using option 1 for accumulation of blocks of
c     equations. in that case mrows is an output variable only, and the
c     matrix data for [e:f] is placed in w(*,*), one block of rows at a
c     time.  mrows contains the number of rows in the matrix after
c     triangularizing several blocks of equations. this is an output
c     parameter only when option 1 is used. see iopt(*) contents
c     for details about option 1.
c
c    ------------------
c    bl(*),bu(*),ind(*)
c    ------------------
c     these arrays contain the information about the bounds that the
c     solution values are to satisfy. the value of ind(j) tells the
c     type of bound and bl(j) and bu(j) give the explicit values for
c     the respective upper and lower bounds.
c
c    1.    for ind(j)=1, require x(j) .ge. bl(j).
c          (the value of bu(j) is not used.)
c    2.    for ind(j)=2, require x(j) .le. bu(j).
c          (the value of bl(j) is not used.)
c    3.    for ind(j)=3, require x(j) .ge. bl(j) and
c                                x(j) .le. bu(j).
c    4.    for ind(j)=4, no bounds on x(j) are required.
c          (the values of bl(j) and bu(j) are not used.)
c
c     values other than 1,2,3 or 4 for ind(j) are errors. in the case
c     ind(j)=3 (upper and lower bounds) the condition bl(j) .gt. bu(j)
c     is an error.
c
c    -------
c    iopt(*)
c    -------
c     this is the array where the user can specify nonstandard options
c     for dbolsm( ). most of the time this feature can be ignored by
c     setting the input value iopt(1)=99. occasionally users may have
c     needs that require use of the following subprogram options. for
c     details about how to use the options see below: iopt(*) contents.
c
c     option number   brief statement of purpose
c     ------ ------   ----- --------- -- -------
c           1         return to user for accumulation of blocks
c                     of least squares equations.
c           2         check lengths of all arrays used in the
c                     subprogram.
c           3         standard scaling of the data matrix, e.
c           4         user provides column scaling for matrix e.
c           5         provide option array to the low-level
c                     subprogram dbolsm( ).
c           6         move the iopt(*) processing pointer.
c          99         no more options to change.
c
c    ----
c    x(*)
c    ----
c     this array is used to pass data associated with option 4. ignore
c     this parameter if this option is not used. otherwise see below:
c     iopt(*) contents.
c
c    output
c    ------
c
c    ----------
c    x(*),rnorm
c    ----------
c     the array x(*) contains a solution (if mode .ge.0 or .eq.-22) for
c     the constrained least squares problem. the value rnorm is the
c     minimum residual vector length.
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
c     of the cases -37,-36,...,-22, or -17,...,-2. values .lt. -1
c     correspond to an abnormal completion of the subprogram. to
c     understand the abnormal completion codes see below: error
c     messages for dbols( ). an approximate solution will be returned
c     to the user only when max. iterations is reached, mode=-22.
c     values for mode=-37,...,-22 come from the low-level subprogram
c     dbolsm(). see the section error messages for dbolsm() in the
c     documentation for dbolsm().
c
c    -----------
c    rw(*),iw(*)
c    -----------
c     these are working arrays with 5*ncols and 2*ncols entries.
c     (normally the user can ignore the contents of these arrays,
c     but they must be dimensioned properly.)
c
c    iopt(*) contents
c    ------- --------
c     the option array allows a user to modify internal variables in
c     the subprogram without recompiling the source code. a central
c     goal of the initial software design was to do a good job for most
c     people. thus the use of options will be restricted to a select
c     group of users. the processing of the option array proceeds as
c     follows: a pointer, here called lp, is initially set to the value
c     1. this value is updated as each option is processed. at the
c     pointer position the option number is extracted and used for
c     locating other information that allows for options to be changed.
c     the portion of the array iopt(*) that is used for each option is
c     fixed; the user and the subprogram both know how many locations
c     are needed for each option. a great deal of error checking is
c     done by the subprogram on the contents of the option array.
c     nevertheless it is still possible to give the subprogram optional
c     input that is meaningless. for example option 4 uses the
c     locations x(ncols+ioff),...,x(ncols+ioff+ncols-1) for passing
c     scaling data. the user must manage the allocation of these
c     locations.
c
c   1
c   -
c     this option allows the user to solve problems with a large number
c     of rows compared to the number of variables. the idea is that the
c     subprogram returns to the user (perhaps many times) and receives
c     new least squares equations from the calling program unit.
c     eventually the user signals "that's all" and then computes the
c     solution with one final call to subprogram dbols( ). the value of
c     mrows is an output variable when this option is used. its value
c     is always in the range 0 .le. mrows .le. ncols+1. it is equal to
c     the number of rows after the triangularization of the entire set
c     of equations. if lp is the processing pointer for iopt(*), the
c     usage for the sequential processing of blocks of equations is
c
c        iopt(lp)=1
c        move block of equations to w(*,*) starting at
c        the first row of w(*,*).
c        iopt(lp+3)=# of rows in the block; user defined
c
c     the user now calls dbols( ) in a loop. the value of iopt(lp+1)
c     directs the user's action. the value of iopt(lp+2) points to
c     where the subsequent rows are to be placed in w(*,*).
c
c      .<loop
c      . call dbols()
c      . if(iopt(lp+1) .eq. 1) then
c      .    iopt(lp+3)=# of rows in the new block; user defined
c      .    place new block of iopt(lp+3) rows in
c      .    w(*,*) starting at row iopt(lp+2).
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
c
c   2
c   -
c     this option is useful for checking the lengths of all arrays used
c     by dbols() against their actual requirements for this problem.
c     the idea is simple: the user's program unit passes the declared
c     dimension information of the arrays. these values are compared
c     against the problem-dependent needs within the subprogram. if any
c     of the dimensions are too small an error message is printed and a
c     negative value of mode is returned, -11 to -17. the printed error
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
c        call dbols()
c
c     use of this option adds 8 to the required length of iopt(*).
c
c   3
c   -
c     this option changes the type of scaling for the data matrix e.
c     nominally each nonzero column of e is scaled so that the
c     magnitude of its largest entry is equal to the value one. if lp
c     is the processing pointer for iopt(*),
c
c        iopt(lp)=3
c        iopt(lp+1)=1,2 or 3
c            1= nominal scaling as noted;
c            2= each nonzero column scaled to have length one;
c            3= identity scaling; scaling effectively suppressed.
c         .
c        call dbols()
c
c     use of this option adds 2 to the required length of iopt(*).
c
c   4
c   -
c     this option allows the user to provide arbitrary (positive)
c     column scaling for the matrix e. if lp is the processing pointer
c     for iopt(*),
c
c        iopt(lp)=4
c        iopt(lp+1)=ioff
c        x(ncols+ioff),...,x(ncols+ioff+ncols-1)
c        = positive scale factors for cols. of e.
c         .
c        call dbols()
c
c     use of this option adds 2 to the required length of iopt(*) and
c     ncols to the required length of x(*).
c
c   5
c   -
c     this option allows the user to provide an option array to the
c     low-level subprogram dbolsm(). if lp is the processing pointer
c     for iopt(*),
c
c        iopt(lp)=5
c        iopt(lp+1)= position in iopt(*) where option array
c                    data for dbolsm() begins.
c         .
c        call dbols()
c
c     use of this option adds 2 to the required length of iopt(*).
c
c   6
c   -
c     move the processing pointer (either forward or backward) to the
c     location iopt(lp+1). the processing point is moved to entry
c     lp+2 of iopt(*) if the option is left with -6 in iopt(lp).  for
c     example to skip over locations 3,...,ncols+2 of iopt(*),
c
c       iopt(1)=6
c       iopt(2)=ncols+3
c       (iopt(i), i=3,...,ncols+2 are not defined here.)
c       iopt(ncols+3)=99
c       call dbols()
c
c     caution: misuse of this option can yield some very hard
c     -to-find bugs.  use it with care.
c
c   99
c   --
c     there are no more options to change.
c
c     only option numbers -99, -6,-5,...,-1, 1,2,...,6, and 99 are
c     permitted. other values are errors. options -99,-1,...,-6 mean
c     that the respective options 99,1,...,6 are left at their default
c     values. an example is the option to modify the (rank) tolerance:
c
c       iopt(1)=-3 option is recognized but not changed
c       iopt(2)=2  scale nonzero cols. to have length one
c       iopt(3)=99
c
c    error messages for dbols()
c    ----- -------- --- -------
c
c warning in...
c dbols(). mdw=(i1) must be positive.
c           in above message, i1=         0
c error number =         2
c (normally a return to the user takes place following this message.)
c
c warning in...
c dbols(). ncols=(i1) the no. of variables must be positive.
c           in above message, i1=         0
c error number =         3
c (normally a return to the user takes place following this message.)
c
c warning in...
c dbols(). for j=(i1), ind(j)=(i2) must be 1-4.
c           in above message, i1=         1
c           in above message, i2=         0
c error number =         4
c (normally a return to the user takes place following this message.)
c
c warning in...
c dbols(). for j=(i1), bound bl(j)=(r1) is .gt. bu(j)=(r2).
c           in above message, i1=         1
c           in above message, r1=    0.
c           in above message, r2=    above message, i1=         0
c error number =         6
c (normally a return to the user takes place following this message.)
c
c warning in...
c dbols(). iscale option=(i1) must be 1-3.
c           in above message, i1=         0
c error number =         7
c (normally a return to the user takes place following this message.)
c
c warning in...
c dbols(). offset past x(ncols) (i1) for user-provided  column scaling
c must be positive.
c           in above message, i1=         0
c error number =         8
c (normally a return to the user takes place following this message.)
c
c warning in...
c dbols(). each provided col. scale factor must be positive.
c component (i1) now = (r1).
c           in above message, i1=        nd. .le. mdw=(i2).
c           in above message, i1=         1
c           in above message, i2=         0
c error number =        10
c (normally a return to the user takes place following this message.)
c
c warning in...
c dbols().the row dimension of w(,)=(i1) must be .ge.the number of rows=
c (i2).
c           in above message, i1=         0
c           in above message, i2=         1
c error number =        11
c (normally a return to the user takes place following this message.)
c
c warning in...
c dbols(). the column dimension of w(,)=(i1) must be .ge. ncols+1=(i2).
c           in above message, i1=         0
c           in above message, i2=         2
c error number =        12
c (normally a return to the user takes place following this message.)
c
c warning in...
c dbols().the dimensions of the arrays bl(),bu(), and ind()=(i1) must be
c .ge. ncols=(i2).
c           in above message, i1=         0
c           in above message, i2=         1
c error number =        13
c (normally a return to the user takes place following this message.)
c
c warning in...
c dbols(). the dimension of x()=(i1) must be .ge. the reqd. length=(i2).
c           in above message, i1=         0
c           in above message, i2=         2
c error number =        14
c (normally a return to the user takes place following this message.)
c
c warning in...
c dbols(). the dimension of rw()=(i1) must be .ge. 5*ncols=(i2).
c           in above message, i1=         0
c           in above message, i2=         3
c error number =        15
c (normally a return to the user takes place following this message.)
c
c warning in...
c dbols() the dimension of iw()=(i1) must be .ge. 2*ncols=(i2).
c           in above message, i1=         0
c           in above message, i2=         2
c error number =        16
c (normally a return to the user takes place following this message.)
c
c warning in...
c dbols() the dimension of iopt()=(i1) must be .ge. the reqd. len.=(i2).
c           in above message, i1=         0
c           in above message, i2=         1
c error number =        17
c (normally a return to the user takes place following this message.)
c
c***references  r. j. hanson, linear least squares with bounds and
c                 linear constraints, report sand82-1517, sandia
c                 laboratories, august 1982.
c***routines called  dbolsm, dcopy, dnrm2, drot, drotg, idamax, xermsg
c***revision history  (yymmdd)
c   821220  date written
c   891006  cosmetic changes to prologue.  (wrb)
c   891006  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900510  convert xerrwv calls to xermsg calls.  (rwc)
c   920501  reformatted the references section.  (wrb)
c***end prologue  dbols
c
c     solve linear least squares system with bounds on
c     selected variables.
c     revised 850329-1400
c     revised yymmdd-hhmm
c     to change this subprogram from single to double precision begin
c     editing at the card 'c++'.
c     change this subprogram name to dbols and the strings
c     /scopy/ to /dcopy/, /sbol/ to /dbol/,
c     /snrm2/ to /dnrm2/, /isamax/ to /idamax/,
c     /srotg/ to /drotg/, /srot/ to /drot/, /e0/ to /d0/,
c     /real            / to /double precision/.
c ++
      double precision w(mdw,*),bl(*),bu(*),x(*),rw(*)
      double precision sc, ss, one, dnrm2, rnorm, zero
c
c     this variable should remain type real.
      integer ind(*),iopt(*),iw(*)
      logical checkl
      character*8 xern1, xern2
      character*16 xern3, xern4
      save igo,locacc,lopt,iscale
      data igo/0/
c***first executable statement  dbols
      nerr = 0
      mode = 0
      if (igo.eq.0) then
c     do(check validity of input data)
c     procedure(check validity of input data)
c
c     see that mdw is .gt.0. gross check only.
          if (mdw.le.0) then
              write (xern1, '(i8)') mdw
              call xermsg ('slatec', 'dbols', 'mdw = ' // xern1 //
     *           ' must be positive.', 2, 1)
c     do(return to user program unit)
              go to 190
          endif
c
c     see that number of unknowns is positive.
          if (ncols.le.0) then
              write (xern1, '(i8)') ncols
              call xermsg ('slatec', 'dbols', 'ncols = ' // xern1 //
     *           ' the no. of variables must be positive.', 3, 1)
c     do(return to user program unit)
              go to 190
          endif
c
c     see that constraint indicators are all well-defined.
          do 10 j = 1,ncols
              if (ind(j).lt.1 .or. ind(j).gt.4) then
                  write (xern1, '(i8)') j
                  write (xern2, '(i8)') ind(j)
                  call xermsg ('slatec', 'dbols', 'ind(' // xern1 //
     *               ') = ' // xern2 // ' must be 1-4.', 4, 1)
c     do(return to user program unit)
                  go to 190
              endif
   10     continue
c
c     see that bounds are consistent.
          do 20 j = 1,ncols
              if (ind(j).eq.3) then
                  if (bl(j).gt.bu(j)) then
                      write (xern1, '(i8)') j
                      write (xern3, '(1pe15.6)') bl(j)
                      write (xern4, '(1pe15.6)') bu(j)
                      call xermsg ('slatec', 'dbols', 'bound bl(' //
     *                   xern1 // ') = ' // xern3 // ' is .gt. bu(' //
     *                   xern1 // ') = ' // xern4, 5, 1)
c     do(return to user program unit)
                      go to 190
                  endif
              endif
   20     continue
c     end procedure
c     do(process option array)
c     procedure(process option array)
          zero = 0.d0
          one = 1.d0
          checkl = .false.
          lenx = ncols
          iscale = 1
          igo = 2
          lopt = 0
          lp = 0
          lds = 0
   30     continue
          lp = lp + lds
          ip = iopt(lp+1)
          jp = abs(ip)
c
c     test for no more options.
          if (ip.eq.99) then
              if (lopt.eq.0) lopt = lp + 1
              go to 50
          else if (jp.eq.99) then
              lds = 1
              go to 30
          else if (jp.eq.1) then
              if (ip.gt.0) then
c
c     set up direction flag, row stacking pointer
c     location, and location for number of new rows.
                  locacc = lp + 2
c
c                  iopt(locacc-1)=option number for seq. accumulation.
c     contents..   iopt(locacc  )=user direction flag, 1 or 2.
c                  iopt(locacc+1)=row stacking pointer.
c                  iopt(locacc+2)=number of new rows to process.
c     user action with this option..
c      (set up option data for seq. accumulation in iopt(*).
c      must also start process with iopt(locacc)=1.)
c      (move block of equations into w(*,*)  starting at first
c       row of w(*,*).  set iopt(locacc+2)=no. of rows in block.)
c              loop
c              call dbols()
c
c                  if(iopt(locacc) .eq. 1) then
c                      stack equas., starting at row iopt(locacc+1),
c                       into w(*,*).
c                       set iopt(locacc+2)=no. of equas.
c                      if last block of equas., set iopt(locacc)=2.
c                  else if iopt(locacc) .eq. 2) then
c                      (process is over. exit loop.)
c                  else
c                      (error condition. should not happen.)
c                  end if
c              end loop
c              set iopt(locacc-1)=-option number for seq. accumulation.
c              call dbols( )
                  iopt(locacc+1) = 1
                  igo = 1
              endif
              lds = 4
              go to 30
          else if (jp.eq.2) then
              if (ip.gt.0) then
c
c     get actual lengths of arrays for checking against needs.
                  locdim = lp + 2
c
c     lmdw.ge.mrows
c     lndw.ge.ncols+1
c     llb .ge.ncols
c     llx .ge.ncols+extra reqd. in options.
c     llrw.ge.5*ncols
c     lliw.ge.2*ncols
c     liop.ge. amount reqd. for ioption array.
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
                      call xermsg ('slatec', 'dbols', 'iscale option = '
     *                   // xern1 // ' must be 1-3', 7, 1)
c     do(return to user program unit)
                      go to 190
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
                      call xermsg ('slatec', 'dbols',
     *                   'offset past x(ncols) (' // xern1 //
     *                   ') for user-provided column scaling must ' //
     *                   'be positive.',  8, 1)
c     do(return to user program unit)
                      go to 190
                  endif
                  call dcopy(ncols,x(ncols+iopt(lp+2)),1,rw,1)
                  lenx = lenx + ncols
                  do 40 j = 1,ncols
                      if (rw(j).le.zero) then
                          write (xern1, '(i8)') j
                          write (xern3, '(1pe15.6)') rw(j)
                          call xermsg ('slatec', 'dbols',
     *                       'each provided column scale factor ' //
     *                       'must be positive.$$component ' // xern1 //
     *                       ' now = ' // xern3, 9, 1)
                          go to 190
                      endif
   40             continue
              endif
              lds = 2
c     cycle forever
              go to 30
c
c     in this option an option array is provided to dbolsm().
          else if (jp.eq.5) then
              if (ip.gt.0) then
                  lopt = iopt(lp+2)
              endif
              lds = 2
c     cycle forever
              go to 30
c
c     this option uses the next loc of iopt(*) as an
c     increment to skip.
          else if (jp.eq.6) then
              if (ip.gt.0) then
                  lp = iopt(lp+2) - 1
                  lds = 0
              else
                  lds = 2
              endif
c     cycle forever
              go to 30
c
c     no valid option number was noted. this is an error condition.
          else
              write (xern1, '(i8)') jp
              call xermsg ('slatec', 'dbols', 'the option number = ' //
     *           xern1 // ' is not defined.', 6, 1)
c     do(return to user program unit)
              go to 190
          endif
   50     continue
c     end procedure
          if (checkl) then
c     do(check lengths of arrays)
c     procedure(check lengths of arrays)
c
c     this feature allows the user to make sure that the
c     arrays are long enough for the intended problem size and use.
              if (lmdw.lt.mrows) then
                  write (xern1, '(i8)') lmdw
                  write (xern2, '(i8)') mrows
                  call xermsg ('slatec', 'dbols',
     *               'the row dimension of w(,) = ' // xern1 //
     *               ' must be .ge. the number of rows = ' // xern2,
     *               11, 1)
c     do(return to user program unit)
                  go to 190
              endif
              if (lndw.lt.ncols+1) then
                  write (xern1, '(i8)') lndw
                  write (xern2, '(i8)') ncols+1
                  call xermsg ('slatec', 'dbols',
     *               'the column dimension of w(,) = ' // xern1 //
     *               ' must be .ge. ncols+1 = ' // xern2, 12, 1)
                  go to 190
              endif
              if (llb.lt.ncols) then
                  write (xern1, '(i8)') llb
                  write (xern2, '(i8)') ncols
                  call xermsg ('slatec', 'dbols',
     *           'the dimensions of the arrays bl(), bu(), and ind() = '
     *               // xern1 // ' must be .ge. ncols = ' // xern2,
     *               13, 1)
c     do(return to user program unit)
                  go to 190
              endif
              if (llx.lt.lenx) then
                  write (xern1, '(i8)') llx
                  write (xern2, '(i8)') lenx
                  call xermsg ('slatec', 'dbols',
     *               'the dimension of x() = ' // xern1 //
     *               ' must be .ge. the required length = ' // xern2,
     *               14, 1)
c     do(return to user program unit)
                  go to 190
              endif
              if (llrw.lt.5*ncols) then
                  write (xern1, '(i8)') llrw
                  write (xern2, '(i8)') 5*ncols
                  call xermsg ('slatec', 'dbols',
     *               'the dimension of rw() = ' // xern1 //
     *               ' must be .ge. 5*ncols = ' // xern2, 15, 1)
c     do(return to user program unit)
                  go to 190
              endif
              if (lliw.lt.2*ncols) then
                  write (xern1, '(i8)') lliw
                  write (xern2, '(i8)') 2*ncols
                  call xermsg ('slatec', 'dbols',
     *               'the dimension of iw() = ' // xern1 //
     *               ' must be .ge. 2*ncols = ' // xern2, 16, 1)
c     do(return to user program unit)
                  go to 190
              endif
              if (liopt.lt.lp+1) then
                  write (xern1, '(i8)') liopt
                  write (xern2, '(i8)') lp+1
                  call xermsg ('slatec', 'dbols',
     *               'the dimension of iopt() = ' // xern1 //
     *               ' must be .ge. the required len = ' // xern2, 17,1)
c     do(return to user program unit)
                  go to 190
              endif
c     end procedure
          endif
      endif
      go to (60,90),igo
      go to 180
c
c     go back to the user for accumulation of least squares
c     equations and directions to quit processing.
c     case 1
   60 continue
c     do(accumulate least squares equations)
c     procedure(accumulate least squares equations)
      mrows = iopt(locacc+1) - 1
      inrows = iopt(locacc+2)
      mnew = mrows + inrows
      if (mnew.lt.0 .or. mnew.gt.mdw) then
          write (xern1, '(i8)') mnew
          write (xern2, '(i8)') mdw
          call xermsg ('slatec', 'dbols', 'no. of rows = ' // xern1 //
     *       ' must be .ge. 0 .and. .le. mdw = ' // xern2, 10, 1)
c     do(return to user program unit)
          go to 190
      endif
      do 80 j = 1,min(ncols+1,mnew)
          do 70 i = mnew,max(mrows,j) + 1,-1
              ibig = idamax(i-j,w(j,j),1) + j - 1
c
c     pivot for increased stability.
              call drotg(w(ibig,j),w(i,j),sc,ss)
              call drot(ncols+1-j,w(ibig,j+1),mdw,w(i,j+1),mdw,sc,ss)
              w(i,j) = zero
   70     continue
   80 continue
      mrows = min(ncols+1,mnew)
      iopt(locacc+1) = mrows + 1
      igo = iopt(locacc)
c     end procedure
      if (igo.eq.2) then
          igo = 0
      endif
      go to 180
c     case 2
   90 continue
c     do(initialize variables and data values)
c     procedure(initialize variables and data values)
      do 150 j = 1,ncols
          go to (100,110,120,130),iscale
          go to 140
  100     continue
c     case 1
c
c     this is the nominal scaling. each nonzero
c     col. has max. norm equal to one.
          ibig = idamax(mrows,w(1,j),1)
          rw(j) = abs(w(ibig,j))
          if (rw(j).eq.zero) then
              rw(j) = one
          else
              rw(j) = one/rw(j)
          endif
          go to 140
  110     continue
c     case 2
c
c     this choice of scaling makes each nonzero column
c     have euclidean length equal to one.
          rw(j) = dnrm2(mrows,w(1,j),1)
          if (rw(j).eq.zero) then
              rw(j) = one
          else
              rw(j) = one/rw(j)
          endif
          go to 140
  120     continue
c     case 3
c
c     this case effectively suppresses scaling by setting
c     the scaling matrix to the identity matrix.
          rw(1) = one
          call dcopy(ncols,rw,0,rw,1)
          go to 160
  130     continue
c     case 4
          go to 160
  140     continue
  150 continue
  160 continue
c     end procedure
c     do(solve bounded least squares problem)
c     procedure(solve bounded least squares problem)
c
c     initialize ibasis(*), j=1,ncols, and ibb(*), j=1,ncols,
c     to =j,and =1, for use in dbolsm( ).
      do 170 j = 1,ncols
          iw(j) = j
          iw(j+ncols) = 1
          rw(3*ncols+j) = bl(j)
          rw(4*ncols+j) = bu(j)
  170 continue
      call dbolsm(w,mdw,mrows,ncols,rw(3*ncols+1),rw(4*ncols+1),ind,
     .            iopt(lopt),x,rnorm,mode,rw(ncols+1),rw(2*ncols+1),rw,
     .            iw,iw(ncols+1))
c     end procedure
      igo = 0
  180 continue
      return
c     procedure(return to user program unit)
  190 if(mode.ge.0)mode = -nerr
      igo = 0
      return
c     end procedure
      end
