*deck dpchsp
      subroutine dpchsp (ic, vc, n, x, f, d, incfd, wk, nwk, ierr)
c***begin prologue  dpchsp
c***purpose  set derivatives needed to determine the hermite represen-
c            tation of the cubic spline interpolant to given data, with
c            specified boundary conditions.
c***library   slatec (pchip)
c***category  e1a
c***type      double precision (pchsp-s, dpchsp-d)
c***keywords  cubic hermite interpolation, pchip,
c             piecewise cubic interpolation, spline interpolation
c***author  fritsch, f. n., (llnl)
c             lawrence livermore national laboratory
c             p.o. box 808  (l-316)
c             livermore, ca  94550
c             fts 532-4275, (510) 422-4275
c***description
c
c          dpchsp:   piecewise cubic hermite spline
c
c     computes the hermite representation of the cubic spline inter-
c     polant to the data given in x and f satisfying the boundary
c     conditions specified by ic and vc.
c
c     to facilitate two-dimensional applications, includes an increment
c     between successive values of the f- and d-arrays.
c
c     the resulting piecewise cubic hermite function may be evaluated
c     by dpchfe or dpchfd.
c
c     note:  this is a modified version of c. de boor's cubic spline
c            routine cubspl.
c
c ----------------------------------------------------------------------
c
c  calling sequence:
c
c        parameter  (incfd = ...)
c        integer  ic(2), n, nwk, ierr
c        double precision  vc(2), x(n), f(incfd,n), d(incfd,n), wk(nwk)
c
c        call  dpchsp (ic, vc, n, x, f, d, incfd, wk, nwk, ierr)
c
c   parameters:
c
c     ic -- (input) integer array of length 2 specifying desired
c           boundary conditions:
c           ic(1) = ibeg, desired condition at beginning of data.
c           ic(2) = iend, desired condition at end of data.
c
c           ibeg = 0  to set d(1) so that the third derivative is con-
c              tinuous at x(2).  this is the "not a knot" condition
c              provided by de boor's cubic spline routine cubspl.
c              < this is the default boundary condition. >
c           ibeg = 1  if first derivative at x(1) is given in vc(1).
c           ibeg = 2  if second derivative at x(1) is given in vc(1).
c           ibeg = 3  to use the 3-point difference formula for d(1).
c                     (reverts to the default b.c. if n.lt.3 .)
c           ibeg = 4  to use the 4-point difference formula for d(1).
c                     (reverts to the default b.c. if n.lt.4 .)
c          notes:
c           1. an error return is taken if ibeg is out of range.
c           2. for the "natural" boundary condition, use ibeg=2 and
c              vc(1)=0.
c
c           iend may take on the same values as ibeg, but applied to
c           derivative at x(n).  in case iend = 1 or 2, the value is
c           given in vc(2).
c
c          notes:
c           1. an error return is taken if iend is out of range.
c           2. for the "natural" boundary condition, use iend=2 and
c              vc(2)=0.
c
c     vc -- (input) real*8 array of length 2 specifying desired boundary
c           values, as indicated above.
c           vc(1) need be set only if ic(1) = 1 or 2 .
c           vc(2) need be set only if ic(2) = 1 or 2 .
c
c     n -- (input) number of data points.  (error return if n.lt.2 .)
c
c     x -- (input) real*8 array of independent variable values.  the
c           elements of x must be strictly increasing:
c                x(i-1) .lt. x(i),  i = 2(1)n.
c           (error return if not.)
c
c     f -- (input) real*8 array of dependent variable values to be
c           interpolated.  f(1+(i-1)*incfd) is value corresponding to
c           x(i).
c
c     d -- (output) real*8 array of derivative values at the data
c           points.  these values will determine the cubic spline
c           interpolant with the requested boundary conditions.
c           the value corresponding to x(i) is stored in
c                d(1+(i-1)*incfd),  i=1(1)n.
c           no other entries in d are changed.
c
c     incfd -- (input) increment between successive values in f and d.
c           this argument is provided primarily for 2-d applications.
c           (error return if  incfd.lt.1 .)
c
c     wk -- (scratch) real*8 array of working storage.
c
c     nwk -- (input) length of work array.
c           (error return if nwk.lt.2*n .)
c
c     ierr -- (output) error flag.
c           normal return:
c              ierr = 0  (no errors).
c           "recoverable" errors:
c              ierr = -1  if n.lt.2 .
c              ierr = -2  if incfd.lt.1 .
c              ierr = -3  if the x-array is not strictly increasing.
c              ierr = -4  if ibeg.lt.0 or ibeg.gt.4 .
c              ierr = -5  if iend.lt.0 of iend.gt.4 .
c              ierr = -6  if both of the above are true.
c              ierr = -7  if nwk is too small.
c               note:  the above errors are checked in the order listed,
c                   and following arguments have **not** been validated.
c             (the d-array has not been changed in any of these cases.)
c              ierr = -8  in case of trouble solving the linear system
c                         for the interior derivative values.
c             (the d-array may have been changed in this case.)
c             (             do **not** use it!                )
c
c***references  carl de boor, a practical guide to splines, springer-
c                 verlag, new york, 1978, pp. 53-59.
c***routines called  dpchdf, xermsg
c***revision history  (yymmdd)
c   820503  date written
c   820804  converted to slatec library version.
c   870707  corrected xerror calls for d.p. name(s).
c   890206  corrected xerror calls.
c   890411  added save statements (vers. 3.2).
c   890703  corrected category record.  (wrb)
c   890831  modified array declarations.  (wrb)
c   891006  cosmetic changes to prologue.  (wrb)
c   891006  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   920429  revised format and order of references.  (wrb,fnf)
c***end prologue  dpchsp
c  programming notes:
c
c     to produce a single precision version, simply:
c        a. change dpchsp to pchsp wherever it occurs,
c        b. change the double precision declarations to real, and
c        c. change the constants zero, half, ... to single precision.
c
c  declare arguments.
c
      integer  ic(2), n, incfd, nwk, ierr
      double precision  vc(2), x(*), f(incfd,*), d(incfd,*), wk(2,*)
c
c  declare local variables.
c
      integer  ibeg, iend, index, j, nm1
      double precision  g, half, one, stemp(3), three, two, xtemp(4),
     *  zero
      save zero, half, one, two, three
      double precision  dpchdf
c
      data  zero /0.d0/, half/.5d0/, one/1.d0/, two/2.d0/, three/3.d0/
c
c  validity-check arguments.
c
c***first executable statement  dpchsp
      if ( n.lt.2 )  go to 5001
      if ( incfd.lt.1 )  go to 5002
      do 1  j = 2, n
         if ( x(j).le.x(j-1) )  go to 5003
    1 continue
c
      ibeg = ic(1)
      iend = ic(2)
      ierr = 0
      if ( (ibeg.lt.0).or.(ibeg.gt.4) )  ierr = ierr - 1
      if ( (iend.lt.0).or.(iend.gt.4) )  ierr = ierr - 2
      if ( ierr.lt.0 )  go to 5004
c
c  function definition is ok -- go on.
c
      if ( nwk .lt. 2*n )  go to 5007
c
c  compute first differences of x sequence and store in wk(1,.). also,
c  compute first divided difference of data and store in wk(2,.).
      do 5  j=2,n
         wk(1,j) = x(j) - x(j-1)
         wk(2,j) = (f(1,j) - f(1,j-1))/wk(1,j)
    5 continue
c
c  set to default boundary conditions if n is too small.
c
      if ( ibeg.gt.n )  ibeg = 0
      if ( iend.gt.n )  iend = 0
c
c  set up for boundary conditions.
c
      if ( (ibeg.eq.1).or.(ibeg.eq.2) )  then
         d(1,1) = vc(1)
      else if (ibeg .gt. 2)  then
c        pick up first ibeg points, in reverse order.
         do 10  j = 1, ibeg
            index = ibeg-j+1
c           index runs from ibeg down to 1.
            xtemp(j) = x(index)
            if (j .lt. ibeg)  stemp(j) = wk(2,index)
   10    continue
c                 --------------------------------
         d(1,1) = dpchdf (ibeg, xtemp, stemp, ierr)
c                 --------------------------------
         if (ierr .ne. 0)  go to 5009
         ibeg = 1
      endif
c
      if ( (iend.eq.1).or.(iend.eq.2) )  then
         d(1,n) = vc(2)
      else if (iend .gt. 2)  then
c        pick up last iend points.
         do 15  j = 1, iend
            index = n-iend+j
c           index runs from n+1-iend up to n.
            xtemp(j) = x(index)
            if (j .lt. iend)  stemp(j) = wk(2,index+1)
   15    continue
c                 --------------------------------
         d(1,n) = dpchdf (iend, xtemp, stemp, ierr)
c                 --------------------------------
         if (ierr .ne. 0)  go to 5009
         iend = 1
      endif
c
c --------------------( begin coding from cubspl )--------------------
c
c  **** a tridiagonal linear system for the unknown slopes s(j) of
c  f  at x(j), j=1,...,n, is generated and then solved by gauss elim-
c  ination, with s(j) ending up in d(1,j), all j.
c     wk(1,.) and wk(2,.) are used for temporary storage.
c
c  construct first equation from first boundary condition, of the form
c             wk(2,1)*s(1) + wk(1,1)*s(2) = d(1,1)
c
      if (ibeg .eq. 0)  then
         if (n .eq. 2)  then
c           no condition at left end and n = 2.
            wk(2,1) = one
            wk(1,1) = one
            d(1,1) = two*wk(2,2)
         else
c           not-a-knot condition at left end and n .gt. 2.
            wk(2,1) = wk(1,3)
            wk(1,1) = wk(1,2) + wk(1,3)
            d(1,1) =((wk(1,2) + two*wk(1,1))*wk(2,2)*wk(1,3)
     *                        + wk(1,2)**2*wk(2,3)) / wk(1,1)
         endif
      else if (ibeg .eq. 1)  then
c        slope prescribed at left end.
         wk(2,1) = one
         wk(1,1) = zero
      else
c        second derivative prescribed at left end.
         wk(2,1) = two
         wk(1,1) = one
         d(1,1) = three*wk(2,2) - half*wk(1,2)*d(1,1)
      endif
c
c  if there are interior knots, generate the corresponding equations and
c  carry out the forward pass of gauss elimination, after which the j-th
c  equation reads    wk(2,j)*s(j) + wk(1,j)*s(j+1) = d(1,j).
c
      nm1 = n-1
      if (nm1 .gt. 1)  then
         do 20 j=2,nm1
            if (wk(2,j-1) .eq. zero)  go to 5008
            g = -wk(1,j+1)/wk(2,j-1)
            d(1,j) = g*d(1,j-1)
     *                  + three*(wk(1,j)*wk(2,j+1) + wk(1,j+1)*wk(2,j))
            wk(2,j) = g*wk(1,j-1) + two*(wk(1,j) + wk(1,j+1))
   20    continue
      endif
c
c  construct last equation from second boundary condition, of the form
c           (-g*wk(2,n-1))*s(n-1) + wk(2,n)*s(n) = d(1,n)
c
c     if slope is prescribed at right end, one can go directly to back-
c     substitution, since arrays happen to be set up just right for it
c     at this point.
      if (iend .eq. 1)  go to 30
c
      if (iend .eq. 0)  then
         if (n.eq.2 .and. ibeg.eq.0)  then
c           not-a-knot at right endpoint and at left endpoint and n = 2.
            d(1,2) = wk(2,2)
            go to 30
         else if ((n.eq.2) .or. (n.eq.3 .and. ibeg.eq.0))  then
c           either (n=3 and not-a-knot also at left) or (n=2 and *not*
c           not-a-knot at left end point).
            d(1,n) = two*wk(2,n)
            wk(2,n) = one
            if (wk(2,n-1) .eq. zero)  go to 5008
            g = -one/wk(2,n-1)
         else
c           not-a-knot and n .ge. 3, and either n.gt.3 or  also not-a-
c           knot at left end point.
            g = wk(1,n-1) + wk(1,n)
c           do not need to check following denominators (x-differences).
            d(1,n) = ((wk(1,n)+two*g)*wk(2,n)*wk(1,n-1)
     *                  + wk(1,n)**2*(f(1,n-1)-f(1,n-2))/wk(1,n-1))/g
            if (wk(2,n-1) .eq. zero)  go to 5008
            g = -g/wk(2,n-1)
            wk(2,n) = wk(1,n-1)
         endif
      else
c        second derivative prescribed at right endpoint.
         d(1,n) = three*wk(2,n) + half*wk(1,n)*d(1,n)
         wk(2,n) = two
         if (wk(2,n-1) .eq. zero)  go to 5008
         g = -one/wk(2,n-1)
      endif
c
c  complete forward pass of gauss elimination.
c
      wk(2,n) = g*wk(1,n-1) + wk(2,n)
      if (wk(2,n) .eq. zero)   go to 5008
      d(1,n) = (g*d(1,n-1) + d(1,n))/wk(2,n)
c
c  carry out back substitution
c
   30 continue
      do 40 j=nm1,1,-1
         if (wk(2,j) .eq. zero)  go to 5008
         d(1,j) = (d(1,j) - wk(1,j)*d(1,j+1))/wk(2,j)
   40 continue
c --------------------(  end  coding from cubspl )--------------------
c
c  normal return.
c
      return
c
c  error returns.
c
 5001 continue
c     n.lt.2 return.
      ierr = -1
      call xermsg ('slatec', 'dpchsp',
     +   'number of data points less than two', ierr, 1)
      return
c
 5002 continue
c     incfd.lt.1 return.
      ierr = -2
      call xermsg ('slatec', 'dpchsp', 'increment less than one', ierr,
     +   1)
      return
c
 5003 continue
c     x-array not strictly increasing.
      ierr = -3
      call xermsg ('slatec', 'dpchsp',
     +   'x-array not strictly increasing', ierr, 1)
      return
c
 5004 continue
c     ic out of range return.
      ierr = ierr - 3
      call xermsg ('slatec', 'dpchsp', 'ic out of range', ierr, 1)
      return
c
 5007 continue
c     nwk too small return.
      ierr = -7
      call xermsg ('slatec', 'dpchsp', 'work array too small', ierr, 1)
      return
c
 5008 continue
c     singular system.
c   *** theoretically, this can only occur if successive x-values   ***
c   *** are equal, which should already have been caught (ierr=-3). ***
      ierr = -8
      call xermsg ('slatec', 'dpchsp', 'singular linear system', ierr,
     +   1)
      return
c
 5009 continue
c     error return from dpchdf.
c   *** this case should never occur ***
      ierr = -9
      call xermsg ('slatec', 'dpchsp', 'error return from dpchdf',
     +   ierr, 1)
      return
c------------- last line of dpchsp follows -----------------------------
      end
