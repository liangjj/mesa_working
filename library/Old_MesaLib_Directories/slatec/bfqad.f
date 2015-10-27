*deck bfqad
      subroutine bfqad (f, t, bcoef, n, k, id, x1, x2, tol, quad, ierr,
     +   work)
c***begin prologue  bfqad
c***purpose  compute the integral of a product of a function and a
c            derivative of a b-spline.
c***library   slatec
c***category  h2a2a1, e3, k6
c***type      single precision (bfqad-s, dbfqad-d)
c***keywords  integral of b-spline, quadrature
c***author  amos, d. e., (snla)
c***description
c
c     abstract
c         bfqad computes the integral on (x1,x2) of a product of a
c         function f and the id-th derivative of a k-th order b-spline,
c         using the b-representation (t,bcoef,n,k).  (x1,x2) must be
c         a subinterval of t(k) .le. x .le. t(n+1).  an integration
c         routine bsgq8 (a modification
c         of gaus8), integrates the product on sub-
c         intervals of (x1,x2) formed by included (distinct) knots.
c
c     description of arguments
c         input
c           f      - external function of one argument for the
c                    integrand bf(x)=f(x)*bvalu(t,bcoef,n,k,id,x,inbv,
c                    work)
c           t      - knot array of length n+k
c           bcoef  - coefficient array of length n
c           n      - length of coefficient array
c           k      - order of b-spline, k .ge. 1
c           id     - order of the spline derivative, 0 .le. id .le. k-1
c                    id=0 gives the spline function
c           x1,x2  - end points of quadrature interval in
c                    t(k) .le. x .le. t(n+1)
c           tol    - desired accuracy for the quadrature, suggest
c                    10.*stol .lt. tol .le. 0.1 where stol is the single
c                    precision unit roundoff for the machine = r1mach(4)
c
c         output
c           quad   - integral of bf(x) on (x1,x2)
c           ierr   - a status code
c                    ierr=1  normal return
c                         2  some quadrature on (x1,x2) does not meet
c                            the requested tolerance.
c           work   - work vector of length 3*k
c
c     error conditions
c         x1 or x2 not in t(k) .le. x .le. t(n+1) is a fatal error.
c         tol not greater than the single precision unit roundoff or
c         less than 0.1 is a fatal error.
c         some quadrature fails to meet the requested tolerance.
c
c***references  d. e. amos, quadrature subroutines for splines and
c                 b-splines, report sand79-1825, sandia laboratories,
c                 december 1979.
c***routines called  bsgq8, intrv, r1mach, xermsg
c***revision history  (yymmdd)
c   800901  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  bfqad
c
c
      integer id, ierr, iflg, ilo, il1, il2, k, left, mflag, n, npk, np1
      real a,aa,ans,b,bb,bcoef,q,quad,t,ta,tb,tol,work,wtol, x1,
     1 x2
      real r1mach, f
      dimension t(*), bcoef(*), work(*)
      external f
c***first executable statement  bfqad
      ierr = 1
      quad = 0.0e0
      if(k.lt.1) go to 100
      if(n.lt.k) go to 105
      if(id.lt.0 .or. id.ge.k) go to 110
      wtol = r1mach(4)
      if (tol.lt.wtol .or. tol.gt.0.1e0) go to 30
      aa = min(x1,x2)
      bb = max(x1,x2)
      if (aa.lt.t(k)) go to 20
      np1 = n + 1
      if (bb.gt.t(np1)) go to 20
      if (aa.eq.bb) return
      npk = n + k
c
      ilo = 1
      call intrv(t, npk, aa, ilo, il1, mflag)
      call intrv(t, npk, bb, ilo, il2, mflag)
      if (il2.ge.np1) il2 = n
      inbv = 1
      q = 0.0e0
      do 10 left=il1,il2
        ta = t(left)
        tb = t(left+1)
        if (ta.eq.tb) go to 10
        a = max(aa,ta)
        b = min(bb,tb)
        call bsgq8(f,t,bcoef,n,k,id,a,b,inbv,tol,ans,iflg,work)
        if (iflg.gt.1) ierr = 2
        q = q + ans
   10 continue
      if (x1.gt.x2) q = -q
      quad = q
      return
c
c
   20 continue
      call xermsg ('slatec', 'bfqad',
     +   'x1 or x2 or both do not satisfy t(k).le.x.le.t(n+1)', 2, 1)
      return
   30 continue
      call xermsg ('slatec', 'bfqad',
     +   'tol is less than the single precision tolerance or ' //
     +   'greater than 0.1', 2, 1)
      return
  100 continue
      call xermsg ('slatec', 'bfqad', 'k does not satisfy k.ge.1', 2,
     +   1)
      return
  105 continue
      call xermsg ('slatec', 'bfqad', 'n does not satisfy n.ge.k', 2,
     +   1)
      return
  110 continue
      call xermsg ('slatec', 'bfqad',
     +   'id does not satisfy 0 .le. id .lt. k', 2, 1)
      return
      end
