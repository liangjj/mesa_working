*deck pchci
      subroutine pchci (n, h, slope, d, incfd)
c***begin prologue  pchci
c***subsidiary
c***purpose  set interior derivatives for pchic
c***library   slatec (pchip)
c***type      single precision (pchci-s, dpchci-d)
c***author  fritsch, f. n., (llnl)
c***description
c
c          pchci:  pchic initial derivative setter.
c
c    called by pchic to set derivatives needed to determine a monotone
c    piecewise cubic hermite interpolant to the data.
c
c    default boundary conditions are provided which are compatible
c    with monotonicity.  if the data are only piecewise monotonic, the
c    interpolant will have an extremum at each point where monotonicity
c    switches direction.
c
c    to facilitate two-dimensional applications, includes an increment
c    between successive values of the d-array.
c
c    the resulting piecewise cubic hermite function should be identical
c    (within roundoff error) to that produced by pchim.
c
c ----------------------------------------------------------------------
c
c  calling sequence:
c
c        parameter  (incfd = ...)
c        integer  n
c        real  h(n), slope(n), d(incfd,n)
c
c        call  pchci (n, h, slope, d, incfd)
c
c   parameters:
c
c     n -- (input) number of data points.
c           if n=2, simply does linear interpolation.
c
c     h -- (input) real array of interval lengths.
c     slope -- (input) real array of data slopes.
c           if the data are (x(i),y(i)), i=1(1)n, then these inputs are:
c                  h(i) =  x(i+1)-x(i),
c              slope(i) = (y(i+1)-y(i))/h(i),  i=1(1)n-1.
c
c     d -- (output) real array of derivative values at the data points.
c           if the data are monotonic, these values will determine a
c           a monotone cubic hermite function.
c           the value corresponding to x(i) is stored in
c                d(1+(i-1)*incfd),  i=1(1)n.
c           no other entries in d are changed.
c
c     incfd -- (input) increment between successive values in d.
c           this argument is provided primarily for 2-d applications.
c
c    -------
c    warning:  this routine does no validity-checking of arguments.
c    -------
c
c  fortran intrinsics used:  abs, max, min.
c
c***see also  pchic
c***routines called  pchst
c***revision history  (yymmdd)
c   820218  date written
c   820601  modified end conditions to be continuous functions of
c           data when monotonicity switches in next interval.
c   820602  1. modified formulas so end conditions are less prone
c             to over/underflow problems.
c           2. minor modification to hsum calculation.
c   820805  converted to slatec library version.
c   890411  added save statements (vers. 3.2).
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910408  updated author section in prologue.  (wrb)
c   930503  improved purpose.  (fnf)
c***end prologue  pchci
c
c  programming notes:
c     1. the function  pchst(arg1,arg2)  is assumed to return zero if
c        either argument is zero, +1 if they are of the same sign, and
c        -1 if they are of opposite sign.
c**end
c
c  declare arguments.
c
      integer  n, incfd
      real  h(*), slope(*), d(incfd,*)
c
c  declare local variables.
c
      integer  i, nless1
      real  del1, del2, dmax, dmin, drat1, drat2, hsum, hsumt3, three,
     *      w1, w2, zero
      save zero, three
      real  pchst
c
c  initialize.
c
      data  zero /0./,  three /3./
c***first executable statement  pchci
      nless1 = n - 1
      del1 = slope(1)
c
c  special case n=2 -- use linear interpolation.
c
      if (nless1 .gt. 1)  go to 10
      d(1,1) = del1
      d(1,n) = del1
      go to 5000
c
c  normal case  (n .ge. 3).
c
   10 continue
      del2 = slope(2)
c
c  set d(1) via non-centered three-point formula, adjusted to be
c     shape-preserving.
c
      hsum = h(1) + h(2)
      w1 = (h(1) + hsum)/hsum
      w2 = -h(1)/hsum
      d(1,1) = w1*del1 + w2*del2
      if ( pchst(d(1,1),del1) .le. zero)  then
         d(1,1) = zero
      else if ( pchst(del1,del2) .lt. zero)  then
c        need do this check only if monotonicity switches.
         dmax = three*del1
         if (abs(d(1,1)) .gt. abs(dmax))  d(1,1) = dmax
      endif
c
c  loop through interior points.
c
      do 50  i = 2, nless1
         if (i .eq. 2)  go to 40
c
         hsum = h(i-1) + h(i)
         del1 = del2
         del2 = slope(i)
   40    continue
c
c        set d(i)=0 unless data are strictly monotonic.
c
         d(1,i) = zero
         if ( pchst(del1,del2) .le. zero)  go to 50
c
c        use brodlie modification of butland formula.
c
         hsumt3 = hsum+hsum+hsum
         w1 = (hsum + h(i-1))/hsumt3
         w2 = (hsum + h(i)  )/hsumt3
         dmax = max( abs(del1), abs(del2) )
         dmin = min( abs(del1), abs(del2) )
         drat1 = del1/dmax
         drat2 = del2/dmax
         d(1,i) = dmin/(w1*drat1 + w2*drat2)
c
   50 continue
c
c  set d(n) via non-centered three-point formula, adjusted to be
c     shape-preserving.
c
      w1 = -h(n-1)/hsum
      w2 = (h(n-1) + hsum)/hsum
      d(1,n) = w1*del1 + w2*del2
      if ( pchst(d(1,n),del2) .le. zero)  then
         d(1,n) = zero
      else if ( pchst(del1,del2) .lt. zero)  then
c        need do this check only if monotonicity switches.
         dmax = three*del2
         if (abs(d(1,n)) .gt. abs(dmax))  d(1,n) = dmax
      endif
c
c  normal return.
c
 5000 continue
      return
c------------- last line of pchci follows ------------------------------
      end
