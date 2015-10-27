*deck sbolsm
      subroutine sbolsm (w, mdw, minput, ncols, bl, bu, ind, iopt, x,
     +   rnorm, mode, rw, ww, scl, ibasis, ibb)
c***begin prologue  sbolsm
c***subsidiary
c***purpose  subsidiary to sbocls and sbols
c***library   slatec
c***type      single precision (sbolsm-s, dbolsm-d)
c***author  (unknown)
c***description
c
c          solve e*x = f (least squares sense) with bounds on
c            selected x values.
c     the user must have dimension statements of the form:
c
c       dimension w(mdw,ncols+1), bl(ncols), bu(ncols),
c      * x(ncols+nx), rw(ncols), ww(ncols), scl(ncols)
c       integer ind(ncols), iopt(1+ni), ibasis(ncols), ibb(ncols)
c
c     (here nx=number of extra locations required for options 1,...,7;
c     nx=0 for no options; here ni=number of extra locations possibly
c     required for options 1-7; ni=0 for no options; ni=14 if all the
c     options are simultaneously in use.)
c
c    input
c    -----
c
c    --------------------
c    w(mdw,*),minput,ncols
c    --------------------
c     the array w(*,*) contains the matrix [e:f] on entry. the matrix
c     [e:f] has minput rows and ncols+1 columns. this data is placed in
c     the array w(*,*) with e occupying the first ncols columns and the
c     right side vector f in column ncols+1. the row dimension, mdw, of
c     the array w(*,*) must satisfy the inequality mdw .ge. minput.
c     other values of mdw are errors. the values of minput and ncols
c     must be positive. other values are errors.
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
c    2.    for ind(j)=2, require x(j) .le. bu(j).
c    3.    for ind(j)=3, require x(j) .ge. bl(j) and
c                                x(j) .le. bu(j).
c    4.    for ind(j)=4, no bounds on x(j) are required.
c     the values of bl(*),bl(*) are modified by the subprogram. values
c     other than 1,2,3 or 4 for ind(j) are errors. in the case ind(j)=3
c     (upper and lower bounds) the condition bl(j) .gt. bu(j) is an
c     error.
c
c    -------
c    iopt(*)
c    -------
c     this is the array where the user can specify nonstandard options
c     for sbolsm. most of the time this feature can be ignored by
c     setting the input value iopt(1)=99. occasionally users may have
c     needs that require use of the following subprogram options. for
c     details about how to use the options see below: iopt(*) contents.
c
c     option number   brief statement of purpose
c     ----- ------   ----- --------- -- -------
c           1         move the iopt(*) processing pointer.
c           2         change rank determination tolerance.
c           3         change blow-up factor that determines the
c                     size of variables being dropped from active
c                     status.
c           4         reset the maximum number of iterations to use
c                     in solving the problem.
c           5         the data matrix is triangularized before the
c                     problem is solved whenever (ncols/minput) .lt.
c                     fac. change the value of fac.
c           6         redefine the weighting matrix used for
c                     linear independence checking.
c           7         debug output is desired.
c          99         no more options to change.
c
c    ----
c    x(*)
c    ----
c     this array is used to pass data associated with options 1,2,3 and
c     5. ignore this input parameter if none of these options are used.
c     otherwise see below: iopt(*) contents.
c
c    ----------------
c    ibasis(*),ibb(*)
c    ----------------
c     these arrays must be initialized by the user. the values
c         ibasis(j)=j, j=1,...,ncols
c         ibb(j)   =1, j=1,...,ncols
c     are appropriate except when using nonstandard features.
c
c    ------
c    scl(*)
c    ------
c     this is the array of scaling factors to use on the columns of the
c     matrix e. these values must be defined by the user. to suppress
c     any column scaling set scl(j)=1.0, j=1,...,ncols.
c
c    output
c    ------
c
c    ----------
c    x(*),rnorm
c    ----------
c     the array x(*) contains a solution (if mode .ge. 0 or .eq. -22)
c     for the constrained least squares problem. the value rnorm is the
c     minimum residual vector length.
c
c    ----
c    mode
c    ----
c     the sign of mode determines whether the subprogram has completed
c     normally, or encountered an error condition or abnormal status.
c     a value of mode .ge. 0 signifies that the subprogram has completed
c     normally. the value of mode (.ge. 0) is the number of variables
c     in an active status: not at a bound nor at the value zero, for
c     the case of free variables. a negative value of mode will be one
c     of the 18 cases -38,-37,...,-22, or -1. values .lt. -1 correspond
c     to an abnormal completion of the subprogram. to understand the
c     abnormal completion codes see below: error messages for sbolsm
c     an approximate solution will be returned to the user only when
c     maximum iterations is reached, mode=-22.
c
c    -----------
c    rw(*),ww(*)
c    -----------
c     these are working arrays each with ncols entries. the array rw(*)
c     contains the working (scaled, nonactive) solution values. the
c     array ww(*) contains the working (scaled, active) gradient vector
c     values.
c
c    ----------------
c    ibasis(*),ibb(*)
c    ----------------
c     these arrays contain information about the status of the solution
c     when mode .ge. 0. the indices ibasis(k), k=1,...,mode, show the
c     nonactive variables; indices ibasis(k), k=mode+1,..., ncols are
c     the active variables. the value (ibb(j)-1) is the number of times
c     variable j was reflected from its upper bound. (normally the user
c     can ignore these parameters.)
c
c    iopt(*) contents
c    ------- --------
c     the option array allows a user to modify internal variables in
c     the subprogram without recompiling the source code. a central
c     goal of the initial software design was to do a good job for most
c     people. thus the use of options will be restricted to a select
c     group of users. the processing of the option array proceeds as
c     follows: a pointer, here called lp, is initially set to the value
c     1. the value is updated as the options are processed.  at the
c     pointer position the option number is extracted and used for
c     locating other information that allows for options to be changed.
c     the portion of the array iopt(*) that is used for each option is
c     fixed; the user and the subprogram both know how many locations
c     are needed for each option. a great deal of error checking is
c     done by the subprogram on the contents of the option array.
c     nevertheless it is still possible to give the subprogram optional
c     input that is meaningless. for example, some of the options use
c     the location x(ncols+ioff) for passing data. the user must manage
c     the allocation of these locations when more than one piece of
c     option data is being passed to the subprogram.
c
c   1
c   -
c     move the processing pointer (either forward or backward) to the
c     location iopt(lp+1). the processing pointer is moved to location
c     lp+2 of iopt(*) in case iopt(lp)=-1.  for example to skip over
c     locations 3,...,ncols+2 of iopt(*),
c
c       iopt(1)=1
c       iopt(2)=ncols+3
c       (iopt(i), i=3,...,ncols+2 are not defined here.)
c       iopt(ncols+3)=99
c       call sbolsm
c
c     caution: misuse of this option can yield some very hard-to-find
c     bugs.  use it with care.
c
c   2
c   -
c     the algorithm that solves the bounded least squares problem
c     iteratively drops columns from the active set. this has the
c     effect of joining a new column vector to the qr factorization of
c     the rectangular matrix consisting of the partially triangularized
c     nonactive columns. after triangularizing this matrix a test is
c     made on the size of the pivot element. the column vector is
c     rejected as dependent if the magnitude of the pivot element is
c     .le. tol* magnitude of the column in components strictly above
c     the pivot element. nominally the value of this (rank) tolerance
c     is tol = sqrt(r1mach(4)). to change only the value of tol, for
c     example,
c
c       x(ncols+1)=tol
c       iopt(1)=2
c       iopt(2)=1
c       iopt(3)=99
c       call sbolsm
c
c     generally, if lp is the processing pointer for iopt(*),
c
c       x(ncols+ioff)=tol
c       iopt(lp)=2
c       iopt(lp+1)=ioff
c        .
c       call sbolsm
c
c     the required length of iopt(*) is increased by 2 if option 2 is
c     used; the required length of x(*) is increased by 1. a value of
c     ioff .le. 0 is an error. a value of tol .le. r1mach(4) gives a
c     warning message; it is not considered an error.
c
c   3
c   -
c     a solution component is left active (not used) if, roughly
c     speaking, it seems too large. mathematically the new component is
c     left active if the magnitude is .ge.((vector norm of f)/(matrix
c     norm of e))/blowup. nominally the factor blowup = sqrt(r1mach(4)).
c     to change only the value of blowup, for example,
c
c       x(ncols+2)=blowup
c       iopt(1)=3
c       iopt(2)=2
c       iopt(3)=99
c       call sbolsm
c
c     generally, if lp is the processing pointer for iopt(*),
c
c       x(ncols+ioff)=blowup
c       iopt(lp)=3
c       iopt(lp+1)=ioff
c        .
c       call sbolsm
c
c     the required length of iopt(*) is increased by 2 if option 3 is
c     used; the required length of x(*) is increased by 1. a value of
c     ioff .le. 0 is an error. a value of blowup .le. 0.0 is an error.
c
c   4
c   -
c     normally the algorithm for solving the bounded least squares
c     problem requires between ncols/3 and ncols drop-add steps to
c     converge. (this remark is based on examining a small number of
c     test cases.) the amount of arithmetic for such problems is
c     typically about twice that required for linear least squares if
c     there are no bounds and if plane rotations are used in the
c     solution method. convergence of the algorithm, while
c     mathematically certain, can be much slower than indicated. to
c     avoid this potential but unlikely event itmax drop-add steps are
c     permitted. nominally itmax=5*(max(minput,ncols)). to change the
c     value of itmax, for example,
c
c       iopt(1)=4
c       iopt(2)=itmax
c       iopt(3)=99
c       call sbolsm
c
c     generally, if lp is the processing pointer for iopt(*),
c
c       iopt(lp)=4
c       iopt(lp+1)=itmax
c        .
c       call sbolsm
c
c     the value of itmax must be .gt. 0. other values are errors. use
c     of this option increases the required length of iopt(*) by 2.
c
c   5
c   -
c     for purposes of increased efficiency the minput by ncols+1 data
c     matrix [e:f] is triangularized as a first step whenever minput
c     satisfies fac*minput .gt. ncols. nominally fac=0.75. to change the
c     value of fac,
c
c       x(ncols+3)=fac
c       iopt(1)=5
c       iopt(2)=3
c       iopt(3)=99
c       call sbolsm
c
c     generally, if lp is the processing pointer for iopt(*),
c
c       x(ncols+ioff)=fac
c       iopt(lp)=5
c       iopt(lp+1)=ioff
c        .
c       call sbolsm
c
c     the value of fac must be nonnegative. other values are errors.
c     resetting fac=0.0 suppresses the initial triangularization step.
c     use of this option increases the required length of iopt(*) by 2;
c     the required length of of x(*) is increased by 1.
c
c   6
c   -
c     the norm used in testing the magnitudes of the pivot element
c     compared to the mass of the column above the pivot line can be
c     changed. the type of change that this option allows is to weight
c     the components with an index larger than mval by the parameter
c     wt. normally mval=0 and wt=1. to change both the values mval and
c     wt, where lp is the processing pointer for iopt(*),
c
c       x(ncols+ioff)=wt
c       iopt(lp)=6
c       iopt(lp+1)=ioff
c       iopt(lp+2)=mval
c
c     use of this option increases the required length of iopt(*) by 3.
c     the length of x(*) is increased by 1. values of mval must be
c     nonnegative and not greater than minput. other values are errors.
c     the value of wt must be positive. any other value is an error. if
c     either error condition is present a message will be printed.
c
c   7
c   -
c     debug output, showing the detailed add-drop steps for the
c     constrained least squares problem, is desired. this option is
c     intended to be used to locate suspected bugs.
c
c   99
c   --
c     there are no more options to change.
c
c     the values for options are 1,...,7,99, and are the only ones
c     permitted. other values are errors. options -99,-1,...,-7 mean
c     that the repective options 99,1,...,7 are left at their default
c     values. an example is the option to modify the (rank) tolerance:
c
c       x(ncols+1)=tol
c       iopt(1)=-2
c       iopt(2)=1
c       iopt(3)=99
c
c    error messages for sbolsm
c    ----- -------- --- ---------
c    -22    more than itmax = ... iterations solving bounded least
c           squares problem.
c
c    -23    the option number = ... is not defined.
c
c    -24    the offset = ... beyond postion ncols = ... must be positive
c           for option number 2.
c
c    -25    the tolerance for rank determination = ... is less than
c           machine precision = ....
c
c    -26    the offset = ... beyond position ncols = ... must be postive
c           for option number 3.
c
c    -27    the reciprocal of the blow-up factor for rejecting variables
c           must be positive. now = ....
c
c    -28    the maximum number of iterations = ... must be positive.
c
c    -29    the offset = ... beyond position ncols = ... must be postive
c           for option number 5.
c
c    -30    the factor (ncols/minput) where pretriangularizing is
c           performed must be nonnegative. now = ....
c
c    -31    the number of rows = ... must be positive.
c
c    -32    the number of columns = ... must be postive.
c
c    -33    the row dimension of w(,) = ... must be .ge. the number of
c           rows = ....
c
c    -34    for j = ... the constraint indicator must be 1-4.
c
c    -35    for j = ... the lower bound = ... is .gt. the upper bound =
c           ....
c
c    -36    the input order of columns = ... is not between 1 and ncols
c           = ....
c
c    -37    the bound polarity flag in component j = ... must be
c           positive. now = ....
c
c    -38    the row separator to apply weighting (...) must lie between
c           0 and minput = .... weight = ... must be positive.
c
c***see also  sbocls, sbols
c***routines called  ivout, r1mach, saxpy, scopy, sdot, smout, snrm2,
c                    srot, srotg, sswap, svout, xermsg
c***revision history  (yymmdd)
c   821220  date written
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   900510  convert xerrwv calls to xermsg calls.  (rwc)
c   920422  fixed usage of minput.  (wrb)
c   901009  editorial changes, code now reads from top to bottom.  (rwc)
c***end prologue  sbolsm
c
c     purpose
c     -------
c     this is the main subprogram that solves the bounded
c     least squares problem.  the problem solved here is:
c
c     solve e*x =  f  (least squares sense)
c     with bounds on selected x values.
c
c     to change this subprogram from single to double precision begin
c     editing at the card 'c++'.
c     change the subprogram name to dbolsm and the strings
c     /saxpy/ to /daxpy/, /scopy/ to /dcopy/,
c     /sdot/ to /ddot/, /snrm2/ to /dnrm2/,
c     /srot/ to /drot/, /srotg/ to /drotg/, /r1mach/ to /d1mach/,
c     /svout/ to /dvout/, /smout/ to /dmout/,
c     /sswap/ to /dswap/, /e0/ to /d0/,
c     /real            / to /double precision/.
c++
c
      real w(mdw,*),bl(*),bu(*)
      real x(*),rw(*),ww(*),scl(*)
      real alpha,beta,bou,colabv,colblo
      real cl1,cl2,cl3,one,big
      real fac,rnorm,sc,ss,t,tolind,wt
      real two,t1,t2,wbig,wlarge,wmag,xnew
      real zero,sdot,snrm2
      real r1mach,tolsze
      integer ibasis(*),ibb(*),ind(*),iopt(*)
      logical found,constr
      character*8 xern1, xern2
      character*16 xern3, xern4
c
      parameter (zero=0.0e0, one=1.0e0, two=2.0e0)
c
      inext(idum) = min(idum+1,mrows)
c***first executable statement  sbolsm
c
c     verify that the problem dimensions are defined properly.
c
      if (minput.le.0) then
          write (xern1, '(i8)') minput
          call xermsg ('slatec', 'sbolsm', 'the number of rows = ' //
     *       xern1 // ' must be positive.', 31, 1)
          mode = -31
          return
      endif
c
      if (ncols.le.0) then
          write (xern1, '(i8)') ncols
          call xermsg ('slatec', 'sbolsm', 'the number of columns = ' //
     *       xern1 // ' must be positive.', 32, 1)
          mode = -32
          return
      endif
c
      if (mdw.lt.minput) then
          write (xern1, '(i8)') mdw
          write (xern2, '(i8)') minput
          call xermsg ('slatec', 'sbolsm',
     *       'the row dimension of w(,) = ' // xern1 //
     *       ' must be .ge. the number of rows = ' // xern2, 33, 1)
          mode = -33
          return
      endif
c
c     verify that bound information is correct.
c
      do 10 j = 1,ncols
          if (ind(j).lt.1 .or. ind(j).gt.4) then
              write (xern1, '(i8)') j
              write (xern2, '(i8)') ind(j)
              call xermsg ('slatec', 'sbolsm', 'for j = ' // xern1 //
     *           ' the constraint indicator must be 1-4', 34, 1)
              mode = -34
              return
          endif
   10 continue
c
      do 20 j = 1,ncols
          if (ind(j).eq.3) then
              if (bu(j).lt.bl(j)) then
                  write (xern1, '(i8)') j
                  write (xern3, '(1pe15.6)') bl(j)
                  write (xern4, '(1pe15.6)') bu(j)
                  call xermsg ('slatec', 'sbolsm', 'for j = ' // xern1
     *               // ' the lower bound = ' // xern3 //
     *               ' is .gt. the upper bound = ' // xern4, 35, 1)
                  mode = -35
                  return
              endif
          endif
   20 continue
c
c     check that permutation and polarity arrays have been set.
c
      do 30 j = 1,ncols
          if (ibasis(j).lt.1 .or. ibasis(j).gt.ncols) then
              write (xern1, '(i8)') ibasis(j)
              write (xern2, '(i8)') ncols
              call xermsg ('slatec', 'sbolsm',
     *           'the input order of columns = ' // xern1 //
     *           ' is not between 1 and ncols = ' // xern2, 36, 1)
              mode = -36
              return
          endif
c
          if (ibb(j).le.0) then
              write (xern1, '(i8)') j
              write (xern2, '(i8)') ibb(j)
              call xermsg ('slatec', 'sbolsm',
     *           'the bound polarity flag in component j = ' // xern1 //
     *           ' must be positive.$$now = ' // xern2, 37, 1)
              mode = -37
              return
          endif
   30 continue
c
c     process the option array.
c
      fac = 0.75e0
      tolind = sqrt(r1mach(4))
      tolsze = sqrt(r1mach(4))
      itmax = 5*max(minput,ncols)
      wt = one
      mval = 0
      iprint = 0
c
c     changes to some parameters can occur through the option array,
c     iopt(*).  process this array looking carefully for input data
c     errors.
c
      lp = 0
      lds = 0
c
c     test for no more options.
c
  590 lp = lp + lds
      ip = iopt(lp+1)
      jp = abs(ip)
      if (ip.eq.99) then
          go to 470
      else if (jp.eq.99) then
          lds = 1
      else if (jp.eq.1) then
c
c         move the iopt(*) processing pointer.
c
          if (ip.gt.0) then
              lp = iopt(lp+2) - 1
              lds = 0
          else
              lds = 2
          endif
      else if (jp.eq.2) then
c
c         change tolerance for rank determination.
c
          if (ip.gt.0) then
              ioff = iopt(lp+2)
              if (ioff.le.0) then
                  write (xern1, '(i8)') ioff
                  write (xern2, '(i8)') ncols
                  call xermsg ('slatec', 'sbolsm', 'the offset = ' //
     *               xern1 // ' beyond position ncols = ' // xern2 //
     *               ' must be positive for option number 2.', 24, 1)
                  mode = -24
                  return
              endif
c
              tolind = x(ncols+ioff)
              if (tolind.lt.r1mach(4)) then
                  write (xern3, '(1pe15.6)') tolind
                  write (xern4, '(1pe15.6)') r1mach(4)
                  call xermsg ('slatec', 'sbolsm',
     *               'the tolerance for rank determination = ' // xern3
     *               // ' is less than machine precision = ' // xern4,
     *               25, 0)
                  mode = -25
              endif
          endif
          lds = 2
      else if (jp.eq.3) then
c
c         change blowup factor for allowing variables to become
c         inactive.
c
          if (ip.gt.0) then
              ioff = iopt(lp+2)
              if (ioff.le.0) then
                  write (xern1, '(i8)') ioff
                  write (xern2, '(i8)') ncols
                  call xermsg ('slatec', 'sbolsm', 'the offset = ' //
     *               xern1 // ' beyond position ncols = ' // xern2 //
     *               ' must be positive for option number 3.', 26, 1)
                  mode = -26
                  return
              endif
c
              tolsze = x(ncols+ioff)
              if (tolsze.le.zero) then
                  write (xern3, '(1pe15.6)') tolsze
                  call xermsg ('slatec', 'sbolsm', 'the reciprocal ' //
     *               'of the blow-up factor for rejecting variables ' //
     *               'must be positive.$$now = ' // xern3, 27, 1)
                  mode = -27
                  return
              endif
          endif
          lds = 2
      else if (jp.eq.4) then
c
c         change the maximum number of iterations allowed.
c
          if (ip.gt.0) then
              itmax = iopt(lp+2)
              if (itmax.le.0) then
                  write (xern1, '(i8)') itmax
                  call xermsg ('slatec', 'sbolsm',
     *               'the maximum number of iterations = ' // xern1 //
     *               ' must be positive.', 28, 1)
                  mode = -28
                  return
              endif
          endif
          lds = 2
      else if (jp.eq.5) then
c
c         change the factor for pretriangularizing the data matrix.
c
          if (ip.gt.0) then
              ioff = iopt(lp+2)
              if (ioff.le.0) then
                  write (xern1, '(i8)') ioff
                  write (xern2, '(i8)') ncols
                  call xermsg ('slatec', 'sbolsm', 'the offset = ' //
     *               xern1 // ' beyond position ncols = ' // xern2 //
     *               ' must be positive for option number 5.', 29, 1)
                  mode = -29
                  return
              endif
c
              fac = x(ncols+ioff)
              if (fac.lt.zero) then
                  write (xern3, '(1pe15.6)') fac
                  call xermsg ('slatec', 'sbolsm',
     *               'the factor (ncols/minput) where pre-' //
     *               'triangularizing is performed must be non-' //
     *               'negative.$$now = ' // xern3, 30, 0)
                  mode = -30
                  return
              endif
          endif
          lds = 2
      else if (jp.eq.6) then
c
c         change the weighting factor (from 1.0) to apply to components
c         numbered .gt. mval (initially set to 1.)  this trick is needed
c         for applications of this subprogram to the heavily weighted
c         least squares problem that come from equality constraints.
c
          if (ip.gt.0) then
              ioff = iopt(lp+2)
              mval = iopt(lp+3)
              wt = x(ncols+ioff)
          endif
c
          if (mval.lt.0 .or. mval.gt.minput .or. wt.le.zero) then
              write (xern1, '(i8)') mval
              write (xern2, '(i8)') minput
              write (xern3, '(1pe15.6)') wt
              call xermsg ('slatec', 'sbolsm',
     *           'the row separator to apply weighting (' // xern1 //
     *           ') must lie between 0 and minput = ' // xern2 //
     *           '.$$weight = ' // xern3 // ' must be positive.', 38, 0)
              mode = -38
              return
          endif
          lds = 3
      else if (jp.eq.7) then
c
c         turn on debug output.
c
          if (ip.gt.0) iprint = 1
          lds = 2
      else
          write (xern1, '(i8)') ip
          call xermsg ('slatec', 'sbolsm', 'the option number = ' //
     *       xern1 // ' is not defined.', 23, 1)
          mode = -23
          return
      endif
      go to 590
c
c     pretriangularize rectangular arrays of certain sizes for
c     increased efficiency.
c
  470 if (fac*minput.gt.ncols) then
          do 490 j = 1,ncols+1
              do 480 i = minput,j+mval+1,-1
                  call srotg(w(i-1,j),w(i,j),sc,ss)
                  w(i,j) = zero
                  call srot(ncols-j+1,w(i-1,j+1),mdw,w(i,j+1),mdw,sc,ss)
  480         continue
  490     continue
          mrows = ncols + mval + 1
      else
          mrows = minput
      endif
c
c     set the x(*) array to zero so all components are defined.
c
      call scopy(ncols,zero,0,x,1)
c
c     the arrays ibasis(*) and ibb(*) are initialized by the calling
c     program and the column scaling is defined in the calling program.
c     'big' is plus infinity on this machine.
c
      big = r1mach(2)
      do 550 j = 1,ncols
          if (ind(j).eq.1) then
              bu(j) = big
          else if (ind(j).eq.2) then
              bl(j) = -big
          else if (ind(j).eq.4) then
              bl(j) = -big
              bu(j) = big
          endif
  550 continue
c
      do 570 j = 1,ncols
          if ((bl(j).le.zero.and.zero.le.bu(j).and.abs(bu(j)).lt.
     *        abs(bl(j))) .or. bu(j).lt.zero) then
              t = bu(j)
              bu(j) = -bl(j)
              bl(j) = -t
              scl(j) = -scl(j)
              do 560 i = 1,mrows
                  w(i,j) = -w(i,j)
  560         continue
          endif
c
c         indices in set t(=tight) are denoted by negative values
c         of ibasis(*).
c
          if (bl(j).ge.zero) then
              ibasis(j) = -ibasis(j)
              t = -bl(j)
              bu(j) = bu(j) + t
              call saxpy(mrows,t,w(1,j),1,w(1,ncols+1),1)
          endif
  570 continue
c
      nsetb = 0
      iter = 0
c
      if (iprint.gt.0) then
          call smout(mrows,ncols+1,mdw,w,'('' pretri. input matrix'')',
     *               -4)
          call svout(ncols,bl,'('' lower bounds'')',-4)
          call svout(ncols,bu,'('' upper bounds'')',-4)
      endif
c
  580 iter = iter + 1
      if (iter.gt.itmax) then
         write (xern1, '(i8)') itmax
         call xermsg ('slatec', 'sbolsm', 'more than itmax = ' // xern1
     *      // ' iterations solving bounded least squares problem.',
     *      22, 1)
         mode = -22
c
c        rescale and translate variables.
c
         igopr = 1
         go to 130
      endif
c
c     find a variable to become non-active.
c                                                 t
c     compute (negative) of gradient vector, w = e *(f-e*x).
c
      call scopy(ncols,zero,0,ww,1)
      do 200 j = nsetb+1,ncols
          jcol = abs(ibasis(j))
          ww(j) = sdot(mrows-nsetb,w(inext(nsetb),j),1,
     *            w(inext(nsetb),ncols+1),1)*abs(scl(jcol))
  200 continue
c
      if (iprint.gt.0) then
          call svout(ncols,ww,'('' gradient values'')',-4)
          call ivout(ncols,ibasis,'('' internal variable order'')',-4)
          call ivout(ncols,ibb,'('' bound polarity'')',-4)
      endif
c
c     if active set = number of total rows, quit.
c
  210 if (nsetb.eq.mrows) then
          found = .false.
          go to 120
      endif
c
c     choose an extremal component of gradient vector for a candidate
c     to become non-active.
c
      wlarge = -big
      wmag = -big
      do 220 j = nsetb+1,ncols
          t = ww(j)
          if (t.eq.big) go to 220
          itemp = ibasis(j)
          jcol = abs(itemp)
          t1 = snrm2(mval-nsetb,w(inext(nsetb),j),1)
          if (itemp.lt.0) then
              if (mod(ibb(jcol),2).eq.0) t = -t
              if (t.lt.zero) go to 220
              if (mval.gt.nsetb) t = t1
              if (t.gt.wlarge) then
                  wlarge = t
                  jlarge = j
              endif
          else
              if (mval.gt.nsetb) t = t1
              if (abs(t).gt.wmag) then
                  wmag = abs(t)
                  jmag = j
              endif
          endif
  220 continue
c
c     choose magnitude of largest component of gradient for candidate.
c
      jbig = 0
      wbig = zero
      if (wlarge.gt.zero) then
          jbig = jlarge
          wbig = wlarge
      endif
c
      if (wmag.ge.wbig) then
          jbig = jmag
          wbig = wmag
      endif
c
      if (jbig.eq.0) then
          found = .false.
          if (iprint.gt.0) then
              call ivout(0,i,'('' found no variable to enter'')',-4)
          endif
          go to 120
      endif
c
c     see if the incoming column is sufficiently independent.  this
c     test is made before an elimination is performed.
c
      if (iprint.gt.0)
     *    call ivout(1,jbig,'('' try to bring in this col.'')',-4)
c
      if (mval.le.nsetb) then
          cl1 = snrm2(mval,w(1,jbig),1)
          cl2 = abs(wt)*snrm2(nsetb-mval,w(inext(mval),jbig),1)
          cl3 = abs(wt)*snrm2(mrows-nsetb,w(inext(nsetb),jbig),1)
          call srotg(cl1,cl2,sc,ss)
          colabv = abs(cl1)
          colblo = cl3
      else
          cl1 = snrm2(nsetb,w(1,jbig),1)
          cl2 = snrm2(mval-nsetb,w(inext(nsetb),jbig),1)
          cl3 = abs(wt)*snrm2(mrows-mval,w(inext(mval),jbig),1)
          colabv = cl1
          call srotg(cl2,cl3,sc,ss)
          colblo = abs(cl2)
      endif
c
      if (colblo.le.tolind*colabv) then
          ww(jbig) = big
          if (iprint.gt.0)
     *        call ivout(0,i,'('' variable is dependent, not used.'')',
     *           -4)
          go to 210
      endif
c
c     swap matrix columns nsetb+1 and jbig, plus pointer information,
c     and gradient values.
c
      nsetb = nsetb + 1
      if (nsetb.ne.jbig) then
          call sswap(mrows,w(1,nsetb),1,w(1,jbig),1)
          call sswap(1,ww(nsetb),1,ww(jbig),1)
          itemp = ibasis(nsetb)
          ibasis(nsetb) = ibasis(jbig)
          ibasis(jbig) = itemp
      endif
c
c     eliminate entries below the pivot line in column nsetb.
c
      if (mrows.gt.nsetb) then
          do 230 i = mrows,nsetb+1,-1
              if (i.eq.mval+1) go to 230
              call srotg(w(i-1,nsetb),w(i,nsetb),sc,ss)
              w(i,nsetb) = zero
              call srot(ncols-nsetb+1,w(i-1,nsetb+1),mdw,w(i,nsetb+1),
     *                  mdw,sc,ss)
  230     continue
c
          if (mval.ge.nsetb .and. mval.lt.mrows) then
              call srotg(w(nsetb,nsetb),w(mval+1,nsetb),sc,ss)
              w(mval+1,nsetb) = zero
              call srot(ncols-nsetb+1,w(nsetb,nsetb+1),mdw,
     *                  w(mval+1,nsetb+1),mdw,sc,ss)
          endif
      endif
c
      if (w(nsetb,nsetb).eq.zero) then
          ww(nsetb) = big
          nsetb = nsetb - 1
          if (iprint.gt.0) then
              call ivout(0,i,'('' pivot is zero, not used.'')',-4)
          endif
          go to 210
      endif
c
c     check that new variable is moving in the right direction.
c
      itemp = ibasis(nsetb)
      jcol = abs(itemp)
      xnew = (w(nsetb,ncols+1)/w(nsetb,nsetb))/abs(scl(jcol))
      if (itemp.lt.0) then
c
c         if(ww(nsetb).ge.zero.and.xnew.le.zero) exit(quit)
c         if(ww(nsetb).le.zero.and.xnew.ge.zero) exit(quit)
c
          if ((ww(nsetb).ge.zero.and.xnew.le.zero) .or.
     *        (ww(nsetb).le.zero.and.xnew.ge.zero)) go to 240
      endif
      found = .true.
      go to 120
c
  240 ww(nsetb) = big
      nsetb = nsetb - 1
      if (iprint.gt.0)
     *    call ivout(0,i,'('' variable has bad direction, not used.'')',
     *       -4)
      go to 210
c
c     solve the triangular system.
c
  270 call scopy(nsetb,w(1,ncols+1),1,rw,1)
      do 280 j = nsetb,1,-1
          rw(j) = rw(j)/w(j,j)
          jcol = abs(ibasis(j))
          t = rw(j)
          if (mod(ibb(jcol),2).eq.0) rw(j) = -rw(j)
          call saxpy(j-1,-t,w(1,j),1,rw,1)
          rw(j) = rw(j)/abs(scl(jcol))
  280 continue
c
      if (iprint.gt.0) then
          call svout(nsetb,rw,'('' soln. values'')',-4)
          call ivout(nsetb,ibasis,'('' cols. used'')',-4)
      endif
c
      if (lgopr.eq.2) then
          call scopy(nsetb,rw,1,x,1)
          do 450 j = 1,nsetb
              itemp = ibasis(j)
              jcol = abs(itemp)
              if (itemp.lt.0) then
                  bou = zero
              else
                  bou = bl(jcol)
              endif
c
              if ((-bou).ne.big) bou = bou/abs(scl(jcol))
              if (x(j).le.bou) then
                  jdrop1 = j
                  go to 340
              endif
c
              bou = bu(jcol)
              if (bou.ne.big) bou = bou/abs(scl(jcol))
              if (x(j).ge.bou) then
                  jdrop2 = j
                  go to 340
              endif
  450     continue
          go to 340
      endif
c
c     see if the unconstrained solution (obtained by solving the
c     triangular system) satisfies the problem bounds.
c
      alpha = two
      beta = two
      x(nsetb) = zero
      do 310 j = 1,nsetb
          itemp = ibasis(j)
          jcol = abs(itemp)
          t1 = two
          t2 = two
          if (itemp.lt.0) then
              bou = zero
          else
              bou = bl(jcol)
          endif
          if ((-bou).ne.big) bou = bou/abs(scl(jcol))
          if (rw(j).le.bou) t1 = (x(j)-bou)/ (x(j)-rw(j))
          bou = bu(jcol)
          if (bou.ne.big) bou = bou/abs(scl(jcol))
          if (rw(j).ge.bou) t2 = (bou-x(j))/ (rw(j)-x(j))
c
c     if not, then compute a step length so that the variables remain
c     feasible.
c
          if (t1.lt.alpha) then
              alpha = t1
              jdrop1 = j
          endif
c
          if (t2.lt.beta) then
              beta = t2
              jdrop2 = j
          endif
  310 continue
c
      constr = alpha .lt. two .or. beta .lt. two
      if (.not.constr) then
c
c         accept the candidate because it satisfies the stated bounds
c         on the variables.
c
          call scopy(nsetb,rw,1,x,1)
          go to 580
      endif
c
c     take a step that is as large as possible with all variables
c     remaining feasible.
c
      do 330 j = 1,nsetb
          x(j) = x(j) + min(alpha,beta)* (rw(j)-x(j))
  330 continue
c
      if (alpha.le.beta) then
          jdrop2 = 0
      else
          jdrop1 = 0
      endif
c
  340 if (jdrop1+jdrop2.le.0 .or. nsetb.le.0) go to 580
  350 jdrop = jdrop1 + jdrop2
      itemp = ibasis(jdrop)
      jcol = abs(itemp)
      if (jdrop2.gt.0) then
c
c         variable is at an upper bound.  subtract multiple of this
c         column from right hand side.
c
          t = bu(jcol)
          if (itemp.gt.0) then
              bu(jcol) = t - bl(jcol)
              bl(jcol) = -t
              itemp = -itemp
              scl(jcol) = -scl(jcol)
              do 360 i = 1,jdrop
                  w(i,jdrop) = -w(i,jdrop)
  360         continue
          else
              ibb(jcol) = ibb(jcol) + 1
              if (mod(ibb(jcol),2).eq.0) t = -t
          endif
c
c     variable is at a lower bound.
c
      else
          if (itemp.lt.zero) then
              t = zero
          else
              t = -bl(jcol)
              bu(jcol) = bu(jcol) + t
              itemp = -itemp
          endif
      endif
c
      call saxpy(jdrop,t,w(1,jdrop),1,w(1,ncols+1),1)
c
c     move certain columns left to achieve upper hessenberg form.
c
      call scopy(jdrop,w(1,jdrop),1,rw,1)
      do 370 j = jdrop+1,nsetb
          ibasis(j-1) = ibasis(j)
          x(j-1) = x(j)
          call scopy(j,w(1,j),1,w(1,j-1),1)
  370 continue
c
      ibasis(nsetb) = itemp
      w(1,nsetb) = zero
      call scopy(mrows-jdrop,w(1,nsetb),0,w(jdrop+1,nsetb),1)
      call scopy(jdrop,rw,1,w(1,nsetb),1)
c
c     transform the matrix from upper hessenberg form to upper
c     triangular form.
c
      nsetb = nsetb - 1
      do 390 i = jdrop,nsetb
c
c         look for small pivots and avoid mixing weighted and
c         nonweighted rows.
c
          if (i.eq.mval) then
              t = zero
              do 380 j = i,nsetb
                  jcol = abs(ibasis(j))
                  t1 = abs(w(i,j)*scl(jcol))
                  if (t1.gt.t) then
                      jbig = j
                      t = t1
                  endif
  380         continue
              go to 400
          endif
          call srotg(w(i,i),w(i+1,i),sc,ss)
          w(i+1,i) = zero
          call srot(ncols-i+1,w(i,i+1),mdw,w(i+1,i+1),mdw,sc,ss)
  390 continue
      go to 430
c
c     the triangularization is completed by giving up the hessenberg
c     form and triangularizing a rectangular matrix.
c
  400 call sswap(mrows,w(1,i),1,w(1,jbig),1)
      call sswap(1,ww(i),1,ww(jbig),1)
      call sswap(1,x(i),1,x(jbig),1)
      itemp = ibasis(i)
      ibasis(i) = ibasis(jbig)
      ibasis(jbig) = itemp
      jbig = i
      do 420 j = jbig,nsetb
          do 410 i = j+1,mrows
              call srotg(w(j,j),w(i,j),sc,ss)
              w(i,j) = zero
              call srot(ncols-j+1,w(j,j+1),mdw,w(i,j+1),mdw,sc,ss)
  410     continue
  420 continue
c
c     see if the remaining coefficients are feasible.  they should be
c     because of the way min(alpha,beta) was chosen.  any that are not
c     feasible will be set to their bounds and appropriately translated.
c
  430 jdrop1 = 0
      jdrop2 = 0
      lgopr = 2
      go to 270
c
c     find a variable to become non-active.
c
  120 if (found) then
          lgopr = 1
          go to 270
      endif
c
c     rescale and translate variables.
c
      igopr = 2
  130 call scopy(nsetb,x,1,rw,1)
      call scopy(ncols,zero,0,x,1)
      do 140 j = 1,nsetb
          jcol = abs(ibasis(j))
          x(jcol) = rw(j)*abs(scl(jcol))
  140 continue
c
      do 150 j = 1,ncols
          if (mod(ibb(j),2).eq.0) x(j) = bu(j) - x(j)
  150 continue
c
      do 160 j = 1,ncols
          jcol = ibasis(j)
          if (jcol.lt.0) x(-jcol) = bl(-jcol) + x(-jcol)
  160 continue
c
      do 170 j = 1,ncols
          if (scl(j).lt.zero) x(j) = -x(j)
  170 continue
c
      i = max(nsetb,mval)
      rnorm = snrm2(mrows-i,w(inext(i),ncols+1),1)
c
      if (igopr.eq.2) mode = nsetb
      return
      end
