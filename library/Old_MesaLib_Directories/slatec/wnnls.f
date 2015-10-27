*deck wnnls
      subroutine wnnls (w, mdw, me, ma, n, l, prgopt, x, rnorm, mode,
     +   iwork, work)
c***begin prologue  wnnls
c***purpose  solve a linearly constrained least squares problem with
c            equality constraints and nonnegativity constraints on
c            selected variables.
c***library   slatec
c***category  k1a2a
c***type      single precision (wnnls-s, dwnnls-d)
c***keywords  constrained least squares, curve fitting, data fitting,
c             equality constraints, inequality constraints,
c             nonnegativity constraints, quadratic programming
c***author  hanson, r. j., (snla)
c           haskell, k. h., (snla)
c***description
c
c     abstract
c
c     this subprogram solves a linearly constrained least squares
c     problem.  suppose there are given matrices e and a of
c     respective dimensions me by n and ma by n, and vectors f
c     and b of respective lengths me and ma.  this subroutine
c     solves the problem
c
c               ex = f, (equations to be exactly satisfied)
c
c               ax = b, (equations to be approximately satisfied,
c                        in the least squares sense)
c
c               subject to components l+1,...,n nonnegative
c
c     any values me.ge.0, ma.ge.0 and 0.le. l .le.n are permitted.
c
c     the problem is reposed as problem wnnls
c
c               (wt*e)x = (wt*f)
c               (   a)    (   b), (least squares)
c               subject to components l+1,...,n nonnegative.
c
c     the subprogram chooses the heavy weight (or penalty parameter) wt.
c
c     the parameters for wnnls are
c
c     input..
c
c     w(*,*),mdw,  the array w(*,*) is double subscripted with first
c     me,ma,n,l    dimensioning parameter equal to mdw.  for this
c                  discussion let us call m = me + ma.  then mdw
c                  must satisfy mdw.ge.m.  the condition mdw.lt.m
c                  is an error.
c
c                  the array w(*,*) contains the matrices and vectors
c
c                       (e  f)
c                       (a  b)
c
c                  in rows and columns 1,...,m and 1,...,n+1
c                  respectively.  columns 1,...,l correspond to
c                  unconstrained variables x(1),...,x(l).  the
c                  remaining variables are constrained to be
c                  nonnegative. the condition l.lt.0 or l.gt.n is
c                  an error.
c
c     prgopt(*)    this real-valued array is the option vector.
c                  if the user is satisfied with the nominal
c                  subprogram features set
c
c                  prgopt(1)=1 (or prgopt(1)=1.0)
c
c                  otherwise prgopt(*) is a linked list consisting of
c                  groups of data of the following form
c
c                  link
c                  key
c                  data set
c
c                  the parameters link and key are each one word.
c                  the data set can be comprised of several words.
c                  the number of items depends on the value of key.
c                  the value of link points to the first
c                  entry of the next group of data within
c                  prgopt(*).  the exception is when there are
c                  no more options to change.  in that
c                  case link=1 and the values key and data set
c                  are not referenced. the general layout of
c                  prgopt(*) is as follows.
c
c               ...prgopt(1)=link1 (link to first entry of next group)
c               .  prgopt(2)=key1 (key to the option change)
c               .  prgopt(3)=data value (data value for this change)
c               .       .
c               .       .
c               .       .
c               ...prgopt(link1)=link2 (link to the first entry of
c               .                       next group)
c               .  prgopt(link1+1)=key2 (key to the option change)
c               .  prgopt(link1+2)=data value
c               ...     .
c               .       .
c               .       .
c               ...prgopt(link)=1 (no more options to change)
c
c                  values of link that are nonpositive are errors.
c                  a value of link.gt.nlink=100000 is also an error.
c                  this helps prevent using invalid but positive
c                  values of link that will probably extend
c                  beyond the program limits of prgopt(*).
c                  unrecognized values of key are ignored.  the
c                  order of the options is arbitrary and any number
c                  of options can be changed with the following
c                  restriction.  to prevent cycling in the
c                  processing of the option array a count of the
c                  number of options changed is maintained.
c                  whenever this count exceeds nopt=1000 an error
c                  message is printed and the subprogram returns.
c
c                  options..
c
c                  key=6
c                         scale the nonzero columns of the
c                  entire data matrix
c                  (e)
c                  (a)
c                  to have length one. the data set for
c                  this option is a single value.  it must
c                  be nonzero if unit length column scaling is
c                  desired.
c
c                  key=7
c                         scale columns of the entire data matrix
c                  (e)
c                  (a)
c                  with a user-provided diagonal matrix.
c                  the data set for this option consists
c                  of the n diagonal scaling factors, one for
c                  each matrix column.
c
c                  key=8
c                         change the rank determination tolerance from
c                  the nominal value of sqrt(srelpr).  this quantity
c                  can be no smaller than srelpr, the arithmetic-
c                  storage precision.  the quantity used
c                  here is internally restricted to be at
c                  least srelpr.  the data set for this option
c                  is the new tolerance.
c
c                  key=9
c                         change the blow-up parameter from the
c                  nominal value of sqrt(srelpr).  the reciprocal of
c                  this parameter is used in rejecting solution
c                  components as too large when a variable is
c                  first brought into the active set.  too large
c                  means that the proposed component times the
c                  reciprocal of the parameter is not less than
c                  the ratio of the norms of the right-side
c                  vector and the data matrix.
c                  this parameter can be no smaller than srelpr,
c                  the arithmetic-storage precision.
c
c                  for example, suppose we want to provide
c                  a diagonal matrix to scale the problem
c                  matrix and change the tolerance used for
c                  determining linear dependence of dropped col
c                  vectors.  for these options the dimensions of
c                  prgopt(*) must be at least n+6.  the fortran
c                  statements defining these options would
c                  be as follows.
c
c                  prgopt(1)=n+3 (link to entry n+3 in prgopt(*))
c                  prgopt(2)=7 (user-provided scaling key)
c
c                  call scopy(n,d,1,prgopt(3),1) (copy the n
c                  scaling factors from a user array called d(*)
c                  into prgopt(3)-prgopt(n+2))
c
c                  prgopt(n+3)=n+6 (link to entry n+6 of prgopt(*))
c                  prgopt(n+4)=8 (linear dependence tolerance key)
c                  prgopt(n+5)=... (new value of the tolerance)
c
c                  prgopt(n+6)=1 (no more options to change)
c
c
c     iwork(1),    the amounts of working storage actually allocated
c     iwork(2)     for the working arrays work(*) and iwork(*),
c                  respectively.  these quantities are compared with
c                  the actual amounts of storage needed for wnnls( ).
c                  insufficient storage allocated for either work(*)
c                  or iwork(*) is considered an error.  this feature
c                  was included in wnnls( ) because miscalculating
c                  the storage formulas for work(*) and iwork(*)
c                  might very well lead to subtle and hard-to-find
c                  execution errors.
c
c                  the length of work(*) must be at least
c
c                  lw = me+ma+5*n
c                  this test will not be made if iwork(1).le.0.
c
c                  the length of iwork(*) must be at least
c
c                  liw = me+ma+n
c                  this test will not be made if iwork(2).le.0.
c
c     output..
c
c     x(*)         an array dimensioned at least n, which will
c                  contain the n components of the solution vector
c                  on output.
c
c     rnorm        the residual norm of the solution.  the value of
c                  rnorm contains the residual vector length of the
c                  equality constraints and least squares equations.
c
c     mode         the value of mode indicates the success or failure
c                  of the subprogram.
c
c                  mode = 0  subprogram completed successfully.
c
c                       = 1  max. number of iterations (equal to
c                            3*(n-l)) exceeded. nearly all problems
c                            should complete in fewer than this
c                            number of iterations. an approximate
c                            solution and its corresponding residual
c                            vector length are in x(*) and rnorm.
c
c                       = 2  usage error occurred.  the offending
c                            condition is noted with the error
c                            processing subprogram, xermsg( ).
c
c     user-designated
c     working arrays..
c
c     work(*)      a real-valued working array of length at least
c                  m + 5*n.
c
c     iwork(*)     an integer-valued working array of length at least
c                  m+n.
c
c***references  k. h. haskell and r. j. hanson, an algorithm for
c                 linear least squares problems with equality and
c                 nonnegativity constraints, report sand77-0552, sandia
c                 laboratories, june 1978.
c               k. h. haskell and r. j. hanson, selected algorithms for
c                 the linearly constrained least squares problem - a
c                 users guide, report sand78-1290, sandia laboratories,
c                 august 1979.
c               k. h. haskell and r. j. hanson, an algorithm for
c                 linear least squares problems with equality and
c                 nonnegativity constraints, mathematical programming
c                 21 (1981), pp. 98-118.
c               r. j. hanson and k. h. haskell, two algorithms for the
c                 linearly constrained least squares problem, acm
c                 transactions on mathematical software, september 1982.
c               c. l. lawson and r. j. hanson, solving least squares
c                 problems, prentice-hall, inc., 1974.
c***routines called  wnlsm, xermsg
c***revision history  (yymmdd)
c   790701  date written
c   890206  revision date from version 3.2
c   890618  completely restructured and revised.  (wrb & rwc)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900510  convert xerrwv calls to xermsg calls.  (rwc)
c   920501  reformatted the references section.  (wrb)
c***end prologue  wnnls
      real              prgopt(*), rnorm, w(mdw,*), work(*), x(*)
      integer iwork(*)
      character*8 xern1
c
c
c***first executable statement  wnnls
      mode = 0
      if (ma+me.le.0 .or. n.le.0) return
      if (iwork(1).gt.0) then
         lw = me + ma + 5*n
         if (iwork(1).lt.lw) then
            write (xern1, '(i8)') lw
            call xermsg ('slatec', 'wnnls', 'insufficient storage ' //
     *         'allocated for work(*), need lw = ' // xern1, 2, 1)
            mode = 2
            return
         endif
      endif
c
      if (iwork(2).gt.0) then
         liw = me + ma + n
         if (iwork(2).lt.liw) then
            write (xern1, '(i8)') liw
            call xermsg ('slatec', 'wnnls', 'insufficient storage ' //
     *         'allocated for iwork(*), need liw = ' // xern1, 2, 1)
            mode = 2
            return
         endif
      endif
c
      if (mdw.lt.me+ma) then
         call xermsg ('slatec', 'wnnls',
     *      'the value mdw.lt.me+ma is an error', 1, 1)
         mode = 2
         return
      endif
c
      if (l.lt.0 .or. l.gt.n) then
         call xermsg ('slatec', 'wnnls',
     *      'l.ge.0 .and. l.le.n is required', 2, 1)
         mode = 2
         return
      endif
c
c     the purpose of this subroutine is to break up the arrays
c     work(*) and iwork(*) into separate work arrays
c     required by the main subroutine wnlsm( ).
c
      l1 = n + 1
      l2 = l1 + n
      l3 = l2 + me + ma
      l4 = l3 + n
      l5 = l4 + n
c
      call wnlsm(w, mdw, me, ma, n, l, prgopt, x, rnorm, mode, iwork,
     *           iwork(l1), work(1), work(l1), work(l2), work(l3),
     *           work(l4), work(l5))
      return
      end
