*deck pchce
      subroutine pchce (ic, vc, n, x, h, slope, d, incfd, ierr)
c***begin prologue  pchce
c***subsidiary
c***purpose  set boundary conditions for pchic
c***library   slatec (pchip)
c***type      single precision (pchce-s, dpchce-d)
c***author  fritsch, f. n., (llnl)
c***description
c
c          pchce:  pchic end derivative setter.
c
c    called by pchic to set end derivatives as requested by the user.
c    it must be called after interior derivative values have been set.
c                      -----
c
c    to facilitate two-dimensional applications, includes an increment
c    between successive values of the d-array.
c
c ----------------------------------------------------------------------
c
c  calling sequence:
c
c        parameter  (incfd = ...)
c        integer  ic(2), n, ierr
c        real  vc(2), x(n), h(n), slope(n), d(incfd,n)
c
c        call  pchce (ic, vc, n, x, h, slope, d, incfd, ierr)
c
c   parameters:
c
c     ic -- (input) integer array of length 2 specifying desired
c           boundary conditions:
c           ic(1) = ibeg, desired condition at beginning of data.
c           ic(2) = iend, desired condition at end of data.
c           ( see prologue to pchic for details. )
c
c     vc -- (input) real array of length 2 specifying desired boundary
c           values.  vc(1) need be set only if ic(1) = 2 or 3 .
c                    vc(2) need be set only if ic(2) = 2 or 3 .
c
c     n -- (input) number of data points.  (assumes n.ge.2)
c
c     x -- (input) real array of independent variable values.  (the
c           elements of x are assumed to be strictly increasing.)
c
c     h -- (input) real array of interval lengths.
c     slope -- (input) real array of data slopes.
c           if the data are (x(i),y(i)), i=1(1)n, then these inputs are:
c                  h(i) =  x(i+1)-x(i),
c              slope(i) = (y(i+1)-y(i))/h(i),  i=1(1)n-1.
c
c     d -- (input) real array of derivative values at the data points.
c           the value corresponding to x(i) must be stored in
c                d(1+(i-1)*incfd),  i=1(1)n.
c          (output) the value of d at x(1) and/or x(n) is changed, if
c           necessary, to produce the requested boundary conditions.
c           no other entries in d are changed.
c
c     incfd -- (input) increment between successive values in d.
c           this argument is provided primarily for 2-d applications.
c
c     ierr -- (output) error flag.
c           normal return:
c              ierr = 0  (no errors).
c           warning errors:
c              ierr = 1  if ibeg.lt.0 and d(1) had to be adjusted for
c                        monotonicity.
c              ierr = 2  if iend.lt.0 and d(1+(n-1)*incfd) had to be
c                        adjusted for monotonicity.
c              ierr = 3  if both of the above are true.
c
c    -------
c    warning:  this routine does no validity-checking of arguments.
c    -------
c
c  fortran intrinsics used:  abs.
c
c***see also  pchic
c***routines called  pchdf, pchst, xermsg
c***revision history  (yymmdd)
c   820218  date written
c   820805  converted to slatec library version.
c   870707  minor corrections made to prologue..
c   890411  added save statements (vers. 3.2).
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900328  added type section.  (wrb)
c   910408  updated author section in prologue.  (wrb)
c   930503  improved purpose.  (fnf)
c***end prologue  pchce
c
c  programming notes:
c     1. the function  pchst(arg1,arg2)  is assumed to return zero if
c        either argument is zero, +1 if they are of the same sign, and
c        -1 if they are of opposite sign.
c     2. one could reduce the number of arguments and amount of local
c        storage, at the expense of reduced code clarity, by passing in
c        the array wk (rather than splitting it into h and slope) and
c        increasing its length enough to incorporate stemp and xtemp.
c     3. the two monotonicity checks only use the sufficient conditions.
c        thus, it is possible (but unlikely) for a boundary condition to
c        be changed, even though the original interpolant was monotonic.
c        (at least the result is a continuous function of the data.)
c**end
c
c  declare arguments.
c
      integer  ic(2), n, incfd, ierr
      real  vc(2), x(*), h(*), slope(*), d(incfd,*)
c
c  declare local variables.
c
      integer  ibeg, iend, ierf, index, j, k
      real  half, stemp(3), three, two, xtemp(4), zero
      save zero, half, two, three
      real  pchdf, pchst
c
c  initialize.
c
      data  zero /0./,  half /0.5/,  two /2./,  three /3./
c
c***first executable statement  pchce
      ibeg = ic(1)
      iend = ic(2)
      ierr = 0
c
c  set to default boundary conditions if n is too small.
c
      if ( abs(ibeg).gt.n )  ibeg = 0
      if ( abs(iend).gt.n )  iend = 0
c
c  treat beginning boundary condition.
c
      if (ibeg .eq. 0)  go to 2000
      k = abs(ibeg)
      if (k .eq. 1)  then
c        boundary value provided.
         d(1,1) = vc(1)
      else if (k .eq. 2)  then
c        boundary second derivative provided.
         d(1,1) = half*( (three*slope(1) - d(1,2)) - half*vc(1)*h(1) )
      else if (k .lt. 5)  then
c        use k-point derivative formula.
c        pick up first k points, in reverse order.
         do 10  j = 1, k
            index = k-j+1
c           index runs from k down to 1.
            xtemp(j) = x(index)
            if (j .lt. k)  stemp(j) = slope(index-1)
   10    continue
c                 -----------------------------
         d(1,1) = pchdf (k, xtemp, stemp, ierf)
c                 -----------------------------
         if (ierf .ne. 0)  go to 5001
      else
c        use 'not a knot' condition.
         d(1,1) = ( three*(h(1)*slope(2) + h(2)*slope(1))
     *             - two*(h(1)+h(2))*d(1,2) - h(1)*d(1,3) ) / h(2)
      endif
c
      if (ibeg .gt. 0)  go to 2000
c
c  check d(1,1) for compatibility with monotonicity.
c
      if (slope(1) .eq. zero)  then
         if (d(1,1) .ne. zero)  then
            d(1,1) = zero
            ierr = ierr + 1
         endif
      else if ( pchst(d(1,1),slope(1)) .lt. zero)  then
         d(1,1) = zero
         ierr = ierr + 1
      else if ( abs(d(1,1)) .gt. three*abs(slope(1)) )  then
         d(1,1) = three*slope(1)
         ierr = ierr + 1
      endif
c
c  treat end boundary condition.
c
 2000 continue
      if (iend .eq. 0)  go to 5000
      k = abs(iend)
      if (k .eq. 1)  then
c        boundary value provided.
         d(1,n) = vc(2)
      else if (k .eq. 2)  then
c        boundary second derivative provided.
         d(1,n) = half*( (three*slope(n-1) - d(1,n-1)) +
     *                                           half*vc(2)*h(n-1) )
      else if (k .lt. 5)  then
c        use k-point derivative formula.
c        pick up last k points.
         do 2010  j = 1, k
            index = n-k+j
c           index runs from n+1-k up to n.
            xtemp(j) = x(index)
            if (j .lt. k)  stemp(j) = slope(index)
 2010    continue
c                 -----------------------------
         d(1,n) = pchdf (k, xtemp, stemp, ierf)
c                 -----------------------------
         if (ierf .ne. 0)  go to 5001
      else
c        use 'not a knot' condition.
         d(1,n) = ( three*(h(n-1)*slope(n-2) + h(n-2)*slope(n-1))
     *             - two*(h(n-1)+h(n-2))*d(1,n-1) - h(n-1)*d(1,n-2) )
     *                                                         / h(n-2)
      endif
c
      if (iend .gt. 0)  go to 5000
c
c  check d(1,n) for compatibility with monotonicity.
c
      if (slope(n-1) .eq. zero)  then
         if (d(1,n) .ne. zero)  then
            d(1,n) = zero
            ierr = ierr + 2
         endif
      else if ( pchst(d(1,n),slope(n-1)) .lt. zero)  then
         d(1,n) = zero
         ierr = ierr + 2
      else if ( abs(d(1,n)) .gt. three*abs(slope(n-1)) )  then
         d(1,n) = three*slope(n-1)
         ierr = ierr + 2
      endif
c
c  normal return.
c
 5000 continue
      return
c
c  error return.
c
 5001 continue
c     error return from pchdf.
c   *** this case should never occur ***
      ierr = -1
      call xermsg ('slatec', 'pchce', 'error return from pchdf', ierr,
     +   1)
      return
c------------- last line of pchce follows ------------------------------
      end
