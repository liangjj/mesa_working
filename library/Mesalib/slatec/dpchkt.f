*deck dpchkt
      subroutine dpchkt (n, x, knotyp, t)
c***begin prologue  dpchkt
c***subsidiary
c***purpose  compute b-spline knot sequence for dpchbs.
c***library   slatec (pchip)
c***category  e3
c***type      double precision (pchkt-s, dpchkt-d)
c***author  fritsch, f. n., (llnl)
c***description
c
c     set a knot sequence for the b-spline representation of a pch
c     function with breakpoints x.  all knots will be at least double.
c     endknots are set as:
c        (1) quadruple knots at endpoints if knotyp=0;
c        (2) extrapolate the length of end interval if knotyp=1;
c        (3) periodic if knotyp=2.
c
c  input arguments:  n, x, knotyp.
c  output arguments:  t.
c
c  restrictions/assumptions:
c     1. n.ge.2 .  (not checked)
c     2. x(i).lt.x(i+1), i=1,...,n .  (not checked)
c     3. 0.le.knotyp.le.2 .  (acts like knotyp=0 for any other value.)
c
c***see also  dpchbs
c***routines called  (none)
c***revision history  (yymmdd)
c   870701  date written
c   900405  converted fortran to upper case.
c   900410  converted prologue to slatec 4.0 format.
c   900410  minor cosmetic changes.
c   900430  produced double precision version.
c   930514  changed nknots from an output to an input variable.  (fnf)
c   930604  removed unused variable nknots from argument list.  (fnf)
c***end prologue  dpchkt
c
c*internal notes:
c
c  since this is subsidiary to dpchbs, which validates its input before
c  calling, it is unnecessary for such validation to be done here.
c
c**end
c
c  declare arguments.
c
      integer  n, knotyp
      double precision  x(*), t(*)
c
c  declare local variables.
c
      integer  j, k, ndim
      double precision  hbeg, hend
c***first executable statement  dpchkt
c
c  initialize.
c
      ndim = 2*n
c
c  set interior knots.
c
      j = 1
      do 20  k = 1, n
         j = j + 2
         t(j) = x(k)
         t(j+1) = t(j)
   20 continue
c     assertion:  at this point t(3),...,t(ndim+2) have been set and
c                 j=ndim+1.
c
c  set end knots according to knotyp.
c
      hbeg = x(2) - x(1)
      hend = x(n) - x(n-1)
      if (knotyp.eq.1 )  then
c          extrapolate.
         t(2) = x(1) - hbeg
         t(ndim+3) = x(n) + hend
      else if ( knotyp.eq.2 )  then
c          periodic.
         t(2) = x(1) - hend
         t(ndim+3) = x(n) + hbeg
      else
c          quadruple end knots.
         t(2) = x(1)
         t(ndim+3) = x(n)
      endif
      t(1) = t(2)
      t(ndim+4) = t(ndim+3)
c
c  terminate.
c
      return
c------------- last line of dpchkt follows -----------------------------
      end
