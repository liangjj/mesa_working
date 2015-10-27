*deck lsei
      subroutine lsei (w, mdw, me, ma, mg, n, prgopt, x, rnorme, rnorml,
     +   mode, ws, ip)
c***begin prologue  lsei
c***purpose  solve a linearly constrained least squares problem with
c            equality and inequality constraints, and optionally compute
c            a covariance matrix.
c***library   slatec
c***category  k1a2a, d9
c***type      single precision (lsei-s, dlsei-d)
c***keywords  constrained least squares, curve fitting, data fitting,
c             equality constraints, inequality constraints,
c             quadratic programming
c***author  hanson, r. j., (snla)
c           haskell, k. h., (snla)
c***description
c
c     abstract
c
c     this subprogram solves a linearly constrained least squares
c     problem with both equality and inequality constraints, and, if the
c     user requests, obtains a covariance matrix of the solution
c     parameters.
c
c     suppose there are given matrices e, a and g of respective
c     dimensions me by n, ma by n and mg by n, and vectors f, b and h of
c     respective lengths me, ma and mg.  this subroutine solves the
c     linearly constrained least squares problem
c
c                   ex = f, (e me by n) (equations to be exactly
c                                       satisfied)
c                   ax = b, (a ma by n) (equations to be
c                                       approximately satisfied,
c                                       least squares sense)
c                   gx .ge. h,(g mg by n) (inequality constraints)
c
c     the inequalities gx .ge. h mean that every component of the
c     product gx must be .ge. the corresponding component of h.
c
c     in case the equality constraints cannot be satisfied, a
c     generalized inverse solution residual vector length is obtained
c     for f-ex.  this is the minimal length possible for f-ex.
c
c     any values me .ge. 0, ma .ge. 0, or mg .ge. 0 are permitted.  the
c     rank of the matrix e is estimated during the computation.  we call
c     this value kranke.  it is an output parameter in ip(1) defined
c     below.  using a generalized inverse solution of ex=f, a reduced
c     least squares problem with inequality constraints is obtained.
c     the tolerances used in these tests for determining the rank
c     of e and the rank of the reduced least squares problem are
c     given in sandia tech. rept. sand-78-1290.  they can be
c     modified by the user if new values are provided in
c     the option list of the array prgopt(*).
c
c     the user must dimension all arrays appearing in the call list..
c     w(mdw,n+1),prgopt(*),x(n),ws(2*(me+n)+k+(mg+2)*(n+7)),ip(mg+2*n+2)
c     where k=max(ma+mg,n).  this allows for a solution of a range of
c     problems in the given working space.  the dimension of ws(*)
c     given is a necessary overestimate.  once a particular problem
c     has been run, the output parameter ip(3) gives the actual
c     dimension required for that problem.
c
c     the parameters for lsei( ) are
c
c     input..
c
c     w(*,*),mdw,   the array w(*,*) is doubly subscripted with
c     me,ma,mg,n    first dimensioning parameter equal to mdw.
c                   for this discussion let us call m = me+ma+mg.  then
c                   mdw must satisfy mdw .ge. m.  the condition
c                   mdw .lt. m is an error.
c
c                   the array w(*,*) contains the matrices and vectors
c
c                                  (e  f)
c                                  (a  b)
c                                  (g  h)
c
c                   in rows and columns 1,...,m and 1,...,n+1
c                   respectively.
c
c                   the integers me, ma, and mg are the
c                   respective matrix row dimensions
c                   of e, a and g.  each matrix has n columns.
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
c                  case, link=1 and the values key and data set
c                  are not referenced.  the general layout of
c                  prgopt(*) is as follows.
c
c               ...prgopt(1) = link1 (link to first entry of next group)
c               .  prgopt(2) = key1 (key to the option change)
c               .  prgopt(3) = data value (data value for this change)
c               .       .
c               .       .
c               .       .
c               ...prgopt(link1)   = link2 (link to the first entry of
c               .                       next group)
c               .  prgopt(link1+1) = key2 (key to the option change)
c               .  prgopt(link1+2) = data value
c               ...     .
c               .       .
c               .       .
c               ...prgopt(link) = 1 (no more options to change)
c
c                  values of link that are nonpositive are errors.
c                  a value of link .gt. nlink=100000 is also an error.
c                  this helps prevent using invalid but positive
c                  values of link that will probably extend
c                  beyond the program limits of prgopt(*).
c                  unrecognized values of key are ignored.  the
c                  order of the options is arbitrary and any number
c                  of options can be changed with the following
c                  restriction.  to prevent cycling in the
c                  processing of the option array, a count of the
c                  number of options changed is maintained.
c                  whenever this count exceeds nopt=1000, an error
c                  message is printed and the subprogram returns.
c
c                  options..
c
c                  key=1
c                         compute in w(*,*) the n by n
c                  covariance matrix of the solution variables
c                  as an output parameter.  nominally the
c                  covariance matrix will not be computed.
c                  (this requires no user input.)
c                  the data set for this option is a single value.
c                  it must be nonzero when the covariance matrix
c                  is desired.  if it is zero, the covariance
c                  matrix is not computed.  when the covariance matrix
c                  is computed, the first dimensioning parameter
c                  of the array w(*,*) must satisfy mdw .ge. max(m,n).
c
c                  key=10
c                         suppress scaling of the inverse of the
c                  normal matrix by the scale factor rnorm**2/
c                  max(1, no. of degrees of freedom).  this option
c                  only applies when the option for computing the
c                  covariance matrix (key=1) is used.  with key=1 and
c                  key=10 used as options the unscaled inverse of the
c                  normal matrix is returned in w(*,*).
c                  the data set for this option is a single value.
c                  when it is nonzero no scaling is done.  when it is
c                  zero scaling is done.  the nominal case is to do
c                  scaling so if option (key=1) is used alone, the
c                  matrix will be scaled on output.
c
c                  key=2
c                         scale the nonzero columns of the
c                         entire data matrix.
c                  (e)
c                  (a)
c                  (g)
c
c                  to have length one.  the data set for this
c                  option is a single value.  it must be
c                  nonzero if unit length column scaling
c                  is desired.
c
c                  key=3
c                         scale columns of the entire data matrix
c                  (e)
c                  (a)
c                  (g)
c
c                  with a user-provided diagonal matrix.
c                  the data set for this option consists
c                  of the n diagonal scaling factors, one for
c                  each matrix column.
c
c                  key=4
c                         change the rank determination tolerance for
c                  the equality constraint equations from
c                  the nominal value of sqrt(srelpr).  this quantity can
c                  be no smaller than srelpr, the arithmetic-
c                  storage precision.  the quantity srelpr is the
c                  largest positive number such that t=1.+srelpr
c                  satisfies t .eq. 1.  the quantity used
c                  here is internally restricted to be at
c                  least srelpr.  the data set for this option
c                  is the new tolerance.
c
c                  key=5
c                         change the rank determination tolerance for
c                  the reduced least squares equations from
c                  the nominal value of sqrt(srelpr).  this quantity can
c                  be no smaller than srelpr, the arithmetic-
c                  storage precision.  the quantity used
c                  here is internally restricted to be at
c                  least srelpr.  the data set for this option
c                  is the new tolerance.
c
c                  for example, suppose we want to change
c                  the tolerance for the reduced least squares
c                  problem, compute the covariance matrix of
c                  the solution parameters, and provide
c                  column scaling for the data matrix.  for
c                  these options the dimension of prgopt(*)
c                  must be at least n+9.  the fortran statements
c                  defining these options would be as follows:
c
c                  prgopt(1)=4 (link to entry 4 in prgopt(*))
c                  prgopt(2)=1 (covariance matrix key)
c                  prgopt(3)=1 (covariance matrix wanted)
c
c                  prgopt(4)=7 (link to entry 7 in prgopt(*))
c                  prgopt(5)=5 (least squares equas.  tolerance key)
c                  prgopt(6)=... (new value of the tolerance)
c
c                  prgopt(7)=n+9 (link to entry n+9 in prgopt(*))
c                  prgopt(8)=3 (user-provided column scaling key)
c
c                  call scopy (n, d, 1, prgopt(9), 1)  (copy the n
c                    scaling factors from the user array d(*)
c                    to prgopt(9)-prgopt(n+8))
c
c                  prgopt(n+9)=1 (no more options to change)
c
c                  the contents of prgopt(*) are not modified
c                  by the subprogram.
c                  the options for wnnls( ) can also be included
c                  in this array.  the values of key recognized
c                  by wnnls( ) are 6, 7 and 8.  their functions
c                  are documented in the usage instructions for
c                  subroutine wnnls( ).  normally these options
c                  do not need to be modified when using lsei( ).
c
c     ip(1),       the amounts of working storage actually
c     ip(2)        allocated for the working arrays ws(*) and
c                  ip(*), respectively.  these quantities are
c                  compared with the actual amounts of storage
c                  needed by lsei( ).  insufficient storage
c                  allocated for either ws(*) or ip(*) is an
c                  error.  this feature was included in lsei( )
c                  because miscalculating the storage formulas
c                  for ws(*) and ip(*) might very well lead to
c                  subtle and hard-to-find execution errors.
c
c                  the length of ws(*) must be at least
c
c                  lw = 2*(me+n)+k+(mg+2)*(n+7)
c
c                  where k = max(ma+mg,n)
c                  this test will not be made if ip(1).le.0.
c
c                  the length of ip(*) must be at least
c
c                  lip = mg+2*n+2
c                  this test will not be made if ip(2).le.0.
c
c     output..
c
c     x(*),rnorme,  the array x(*) contains the solution parameters
c     rnorml        if the integer output flag mode = 0 or 1.
c                   the definition of mode is given directly below.
c                   when mode = 0 or 1, rnorme and rnorml
c                   respectively contain the residual vector
c                   euclidean lengths of f - ex and b - ax.  when
c                   mode=1 the equality constraint equations ex=f
c                   are contradictory, so rnorme .ne. 0.  the residual
c                   vector f-ex has minimal euclidean length.  for
c                   mode .ge. 2, none of these parameters is defined.
c
c     mode          integer flag that indicates the subprogram
c                   status after completion.  if mode .ge. 2, no
c                   solution has been computed.
c
c                   mode =
c
c                   0  both equality and inequality constraints
c                      are compatible and have been satisfied.
c
c                   1  equality constraints are contradictory.
c                      a generalized inverse solution of ex=f was used
c                      to minimize the residual vector length f-ex.
c                      in this sense, the solution is still meaningful.
c
c                   2  inequality constraints are contradictory.
c
c                   3  both equality and inequality constraints
c                      are contradictory.
c
c                   the following interpretation of
c                   mode=1,2 or 3 must be made.  the
c                   sets consisting of all solutions
c                   of the equality constraints ex=f
c                   and all vectors satisfying gx .ge. h
c                   have no points in common.  (in
c                   particular this does not say that
c                   each individual set has no points
c                   at all, although this could be the
c                   case.)
c
c                   4  usage error occurred.  the value
c                      of mdw is .lt. me+ma+mg, mdw is
c                      .lt. n and a covariance matrix is
c                      requested, or the option vector
c                      prgopt(*) is not properly defined,
c                      or the lengths of the working arrays
c                      ws(*) and ip(*), when specified in
c                      ip(1) and ip(2) respectively, are not
c                      long enough.
c
c     w(*,*)        the array w(*,*) contains the n by n symmetric
c                   covariance matrix of the solution parameters,
c                   provided this was requested on input with
c                   the option vector prgopt(*) and the output
c                   flag is returned with mode = 0 or 1.
c
c     ip(*)         the integer working array has three entries
c                   that provide rank and working array length
c                   information after completion.
c
c                      ip(1) = rank of equality constraint
c                              matrix.  define this quantity
c                              as kranke.
c
c                      ip(2) = rank of reduced least squares
c                              problem.
c
c                      ip(3) = the amount of storage in the
c                              working array ws(*) that was
c                              actually used by the subprogram.
c                              the formula given above for the length
c                              of ws(*) is a necessary overestimate.
c                              if exactly the same problem matrices
c                              are used in subsequent executions,
c                              the declared dimension of ws(*) can
c                              be reduced to this output value.
c     user designated
c     working arrays..
c
c     ws(*),ip(*)              these are respectively type real
c                              and type integer working arrays.
c                              their required minimal lengths are
c                              given above.
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
c***routines called  h12, lsi, r1mach, sasum, saxpy, scopy, sdot, snrm2,
c                    sscal, sswap, xermsg
c***revision history  (yymmdd)
c   790701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890618  completely restructured and extensively revised (wrb & rwc)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900510  convert xerrwv calls to xermsg calls.  (rwc)
c   920501  reformatted the references section.  (wrb)
c***end prologue  lsei
      integer ip(3), ma, mdw, me, mg, mode, n
      real             prgopt(*), rnorme, rnorml, w(mdw,*), ws(*), x(*)
c
      external h12, lsi, r1mach, sasum, saxpy, scopy, sdot, snrm2,
     *   sscal, sswap, xermsg
      real             r1mach, sasum, sdot, snrm2
c
      real             enorm, fnorm, gam, rb, rn, rnmax, size, sn,
     *   snmax, srelpr, t, tau, uj, up, vj, xnorm, xnrme
      integer i, imax, j, jp1, k, key, kranke, last, lchk, link, m,
     *   mapke1, mdeqc, mend, mep1, n1, n2, next, nlink, nopt, np1,
     *   ntimes
      logical cov, first
      character*8 xern1, xern2, xern3, xern4
      save first, srelpr
c
      data first /.true./
c***first executable statement  lsei
c
c     set the nominal tolerance used in the code for the equality
c     constraint equations.
c
      if (first) srelpr = r1mach(4)
      first = .false.
      tau = sqrt(srelpr)
c
c     check that enough storage was allocated in ws(*) and ip(*).
c
      mode = 4
      if (min(n,me,ma,mg) .lt. 0) then
         write (xern1, '(i8)') n
         write (xern2, '(i8)') me
         write (xern3, '(i8)') ma
         write (xern4, '(i8)') mg
         call xermsg ('slatec', 'lsei', 'all of the variables n, me,' //
     *      ' ma, mg must be .ge. 0$$entered routine with' //
     *      '$$n  = ' // xern1 //
     *      '$$me = ' // xern2 //
     *      '$$ma = ' // xern3 //
     *      '$$mg = ' // xern4, 2, 1)
         return
      endif
c
      if (ip(1).gt.0) then
         lchk = 2*(me+n) + max(ma+mg,n) + (mg+2)*(n+7)
         if (ip(1).lt.lchk) then
            write (xern1, '(i8)') lchk
            call xermsg ('slatec', 'lsei', 'insufficient storage ' //
     *         'allocated for ws(*), need lw = ' // xern1, 2, 1)
            return
         endif
      endif
c
      if (ip(2).gt.0) then
         lchk = mg + 2*n + 2
         if (ip(2).lt.lchk) then
            write (xern1, '(i8)') lchk
            call xermsg ('slatec', 'lsei', 'insufficient storage ' //
     *         'allocated for ip(*), need lip = ' // xern1, 2, 1)
            return
         endif
      endif
c
c     compute number of possible right multiplying householder
c     transformations.
c
      m = me + ma + mg
      if (n.le.0 .or. m.le.0) then
         mode = 0
         rnorme = 0
         rnorml = 0
         return
      endif
c
      if (mdw.lt.m) then
         call xermsg ('slatec', 'lsei', 'mdw.lt.me+ma+mg is an error',
     +      2, 1)
         return
      endif
c
      np1 = n + 1
      kranke = min(me,n)
      n1 = 2*kranke + 1
      n2 = n1 + n
c
c     set nominal values.
c
c     the nominal column scaling used in the code is
c     the identity scaling.
c
      call scopy (n, 1.e0, 0, ws(n1), 1)
c
c     no covariance matrix is nominally computed.
c
      cov = .false.
c
c     process option vector.
c     define bound for number of options to change.
c
      nopt = 1000
      ntimes = 0
c
c     define bound for positive values of link.
c
      nlink = 100000
      last = 1
      link = prgopt(1)
      if (link.eq.0 .or. link.gt.nlink) then
         call xermsg ('slatec', 'lsei',
     +      'the option vector is undefined', 2, 1)
         return
      endif
c
  100 if (link.gt.1) then
         ntimes = ntimes + 1
         if (ntimes.gt.nopt) then
            call xermsg ('slatec', 'lsei',
     +         'the links in the option vector are cycling.', 2, 1)
            return
         endif
c
         key = prgopt(last+1)
         if (key.eq.1) then
            cov = prgopt(last+2) .ne. 0.e0
         elseif (key.eq.2 .and. prgopt(last+2).ne.0.e0) then
            do 110 j = 1,n
               t = snrm2(m,w(1,j),1)
               if (t.ne.0.e0) t = 1.e0/t
               ws(j+n1-1) = t
  110       continue
         elseif (key.eq.3) then
            call scopy (n, prgopt(last+2), 1, ws(n1), 1)
         elseif (key.eq.4) then
            tau = max(srelpr,prgopt(last+2))
         endif
c
         next = prgopt(link)
         if (next.le.0 .or. next.gt.nlink) then
         call xermsg ('slatec', 'lsei',
     +      'the option vector is undefined', 2, 1)
            return
         endif
c
         last = link
         link = next
         go to 100
      endif
c
      do 120 j = 1,n
         call sscal (m, ws(n1+j-1), w(1,j), 1)
  120 continue
c
      if (cov .and. mdw.lt.n) then
         call xermsg ('slatec', 'lsei',
     +      'mdw .lt. n when cov matrix needed, is an error', 2, 1)
         return
      endif
c
c     problem definition and option vector ok.
c
      mode = 0
c
c     compute norm of equality constraint matrix and right side.
c
      enorm = 0.e0
      do 130 j = 1,n
         enorm = max(enorm,sasum(me,w(1,j),1))
  130 continue
c
      fnorm = sasum(me,w(1,np1),1)
      snmax = 0.e0
      rnmax = 0.e0
      do 150 i = 1,kranke
c
c        compute maximum ratio of vector lengths. partition is at
c        column i.
c
         do 140 k = i,me
            sn = sdot(n-i+1,w(k,i),mdw,w(k,i),mdw)
            rn = sdot(i-1,w(k,1),mdw,w(k,1),mdw)
            if (rn.eq.0.e0 .and. sn.gt.snmax) then
               snmax = sn
               imax = k
            elseif (k.eq.i .or. sn*rnmax.gt.rn*snmax) then
               snmax = sn
               rnmax = rn
               imax = k
            endif
  140    continue
c
c        interchange rows if necessary.
c
         if (i.ne.imax) call sswap (np1, w(i,1), mdw, w(imax,1), mdw)
         if (snmax.gt.rnmax*tau**2) then
c
c        eliminate elements i+1,...,n in row i.
c
            call h12 (1, i, i+1, n, w(i,1), mdw, ws(i), w(i+1,1), mdw,
     +                1, m-i)
         else
            kranke = i - 1
            go to 160
         endif
  150 continue
c
c     save diagonal terms of lower trapezoidal matrix.
c
  160 call scopy (kranke, w, mdw+1, ws(kranke+1), 1)
c
c     use householder transformation from left to achieve
c     kranke by kranke upper triangular form.
c
      if (kranke.lt.me) then
         do 170 k = kranke,1,-1
c
c           apply transformation to matrix cols. 1,...,k-1.
c
            call h12 (1, k, kranke+1, me, w(1,k), 1, up, w, 1, mdw, k-1)
c
c           apply to rt side vector.
c
            call h12 (2, k, kranke+1, me, w(1,k), 1, up, w(1,np1), 1, 1,
     +                1)
  170    continue
      endif
c
c     solve for variables 1,...,kranke in new coordinates.
c
      call scopy (kranke, w(1, np1), 1, x, 1)
      do 180 i = 1,kranke
         x(i) = (x(i)-sdot(i-1,w(i,1),mdw,x,1))/w(i,i)
  180 continue
c
c     compute residuals for reduced problem.
c
      mep1 = me + 1
      rnorml = 0.e0
      do 190 i = mep1,m
         w(i,np1) = w(i,np1) - sdot(kranke,w(i,1),mdw,x,1)
         sn = sdot(kranke,w(i,1),mdw,w(i,1),mdw)
         rn = sdot(n-kranke,w(i,kranke+1),mdw,w(i,kranke+1),mdw)
         if (rn.le.sn*tau**2 .and. kranke.lt.n)
     *      call scopy (n-kranke, 0.e0, 0, w(i,kranke+1), mdw)
  190 continue
c
c     compute equality constraint equations residual length.
c
      rnorme = snrm2(me-kranke,w(kranke+1,np1),1)
c
c     move reduced problem data upward if kranke.lt.me.
c
      if (kranke.lt.me) then
         do 200 j = 1,np1
            call scopy (m-me, w(me+1,j), 1, w(kranke+1,j), 1)
  200    continue
      endif
c
c     compute solution of reduced problem.
c
      call lsi(w(kranke+1, kranke+1), mdw, ma, mg, n-kranke, prgopt,
     +         x(kranke+1), rnorml, mode, ws(n2), ip(2))
c
c     test for consistency of equality constraints.
c
      if (me.gt.0) then
         mdeqc = 0
         xnrme = sasum(kranke,w(1,np1),1)
         if (rnorme.gt.tau*(enorm*xnrme+fnorm)) mdeqc = 1
         mode = mode + mdeqc
c
c        check if solution to equality constraints satisfies inequality
c        constraints when there are no degrees of freedom left.
c
         if (kranke.eq.n .and. mg.gt.0) then
            xnorm = sasum(n,x,1)
            mapke1 = ma + kranke + 1
            mend = ma + kranke + mg
            do 210 i = mapke1,mend
               size = sasum(n,w(i,1),mdw)*xnorm + abs(w(i,np1))
               if (w(i,np1).gt.tau*size) then
                  mode = mode + 2
                  go to 290
               endif
  210       continue
         endif
      endif
c
c     replace diagonal terms of lower trapezoidal matrix.
c
      if (kranke.gt.0) then
         call scopy (kranke, ws(kranke+1), 1, w, mdw+1)
c
c        reapply transformation to put solution in original coordinates.
c
         do 220 i = kranke,1,-1
            call h12 (2, i, i+1, n, w(i,1), mdw, ws(i), x, 1, 1, 1)
  220    continue
c
c        compute covariance matrix of equality constrained problem.
c
         if (cov) then
            do 270 j = min(kranke,n-1),1,-1
               rb = ws(j)*w(j,j)
               if (rb.ne.0.e0) rb = 1.e0/rb
               jp1 = j + 1
               do 230 i = jp1,n
                  w(i,j) = rb*sdot(n-j,w(i,jp1),mdw,w(j,jp1),mdw)
  230          continue
c
               gam = 0.5e0*rb*sdot(n-j,w(jp1,j),1,w(j,jp1),mdw)
               call saxpy (n-j, gam, w(j,jp1), mdw, w(jp1,j), 1)
               do 250 i = jp1,n
                  do 240 k = i,n
                     w(i,k) = w(i,k) + w(j,i)*w(k,j) + w(i,j)*w(j,k)
                     w(k,i) = w(i,k)
  240             continue
  250          continue
               uj = ws(j)
               vj = gam*uj
               w(j,j) = uj*vj + uj*vj
               do 260 i = jp1,n
                  w(j,i) = uj*w(i,j) + vj*w(j,i)
  260          continue
               call scopy (n-j, w(j, jp1), mdw, w(jp1,j), 1)
  270       continue
         endif
      endif
c
c     apply the scaling to the covariance matrix.
c
      if (cov) then
         do 280 i = 1,n
            call sscal (n, ws(i+n1-1), w(i,1), mdw)
            call sscal (n, ws(i+n1-1), w(1,i), 1)
  280    continue
      endif
c
c     rescale solution vector.
c
  290 if (mode.le.1) then
         do 300 j = 1,n
            x(j) = x(j)*ws(n1+j-1)
  300    continue
      endif
c
      ip(1) = kranke
      ip(3) = ip(3) + 2*kranke + n
      return
      end
