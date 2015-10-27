*deck dpfqad
      subroutine dpfqad (f, ldc, c, xi, lxi, k, id, x1, x2, tol, quad,
     +   ierr)
c***begin prologue  dpfqad
c***purpose  compute the integral on (x1,x2) of a product of a
c            function f and the id-th derivative of a b-spline,
c            (pp-representation).
c***library   slatec
c***category  h2a2a1, e3, k6
c***type      double precision (pfqad-s, dpfqad-d)
c***keywords  b-spline, data fitting, interpolation, quadrature, splines
c***author  amos, d. e., (snla)
c***description
c
c     abstract    **** a double precision routine ****
c         dpfqad computes the integral on (x1,x2) of a product of a
c         function f and the id-th derivative of a b-spline, using the
c         pp-representation (c,xi,lxi,k).  (x1,x2) is normally a sub-
c         interval of xi(1) .le. x .le. xi(lxi+1).  an integration
c         routine, dppgq8 (a modification of gaus8), integrates the
c         product on subintervals of (x1,x2) formed by the included
c         break points.  integration outside of (xi(1),xi(lxi+1)) is
c         permitted provided f is defined.
c
c         the maximum number of significant digits obtainable in
c         dbsqad is the smaller of 18 and the number of digits
c         carried in double precision arithmetic.
c
c     description of arguments
c         input      f,c,xi,x1,x2,tol are double precision
c           f      - external function of one argument for the
c                    integrand pf(x)=f(x)*dppval(ldc,c,xi,lxi,k,id,x,
c                    inppv)
c           ldc    - leading dimension of matrix c, ldc .ge. k
c           c(i,j) - right taylor derivatives at xi(j), i=1,k , j=1,lxi
c           xi(*)  - break point array of length lxi+1
c           lxi    - number of polynomial pieces
c           k      - order of b-spline, k .ge. 1
c           id     - order of the spline derivative, 0 .le. id .le. k-1
c                    id=0 gives the spline function
c           x1,x2  - end points of quadrature interval, normally in
c                    xi(1) .le. x .le. xi(lxi+1)
c           tol    - desired accuracy for the quadrature, suggest
c                    10.*dtol .lt. tol .le. 0.1 where dtol is the
c                    maximum of 1.0d-18 and double precision unit
c                    roundoff for the machine = d1mach(4)
c
c         output     quad is double precision
c           quad   - integral of pf(x) on (x1,x2)
c           ierr   - a status code
c                    ierr=1 normal return
c                         2 some quadrature does not meet the
c                           requested tolerance
c
c     error conditions
c         improper input is a fatal error.
c         some quadrature does not meet the requested tolerance.
c
c***references  d. e. amos, quadrature subroutines for splines and
c                 b-splines, report sand79-1825, sandia laboratories,
c                 december 1979.
c***routines called  d1mach, dintrv, dppgq8, xermsg
c***revision history  (yymmdd)
c   800901  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  dpfqad
c
      integer id,ierr,iflg,ilo,il1,il2,inppv,k,ldc,left,lxi,mf1,mf2
      double precision a,aa,ans,b,bb,c,q,quad,ta,tb,tol,wtol,xi,x1,x2
      double precision d1mach, f
      dimension xi(*), c(ldc,*)
      external f
c
c***first executable statement  dpfqad
      ierr = 1
      quad = 0.0d0
      if(k.lt.1) go to 100
      if(ldc.lt.k) go to 105
      if(id.lt.0 .or. id.ge.k) go to 110
      if(lxi.lt.1) go to 115
      wtol = d1mach(4)
      wtol = max(wtol,1.0d-18)
      if (tol.lt.wtol .or. tol.gt.0.1d0) go to 20
      aa = min(x1,x2)
      bb = max(x1,x2)
      if (aa.eq.bb) return
      ilo = 1
      call dintrv(xi, lxi, aa, ilo, il1, mf1)
      call dintrv(xi, lxi, bb, ilo, il2, mf2)
      q = 0.0d0
      inppv = 1
      do 10 left=il1,il2
        ta = xi(left)
        a = max(aa,ta)
        if (left.eq.1) a = aa
        tb = bb
        if (left.lt.lxi) tb = xi(left+1)
        b = min(bb,tb)
        call dppgq8(f,ldc,c,xi,lxi,k,id,a,b,inppv,tol,ans,iflg)
        if (iflg.gt.1) ierr = 2
        q = q + ans
   10 continue
      if (x1.gt.x2) q = -q
      quad = q
      return
c
   20 continue
      call xermsg ('slatec', 'dpfqad',
     +   'tol is less dtol or greater than 0.1', 2, 1)
      return
  100 continue
      call xermsg ('slatec', 'dpfqad', 'k does not satisfy k.ge.1', 2,
     +   1)
      return
  105 continue
      call xermsg ('slatec', 'dpfqad', 'ldc does not satisfy ldc.ge.k',
     +   2, 1)
      return
  110 continue
      call xermsg ('slatec', 'dpfqad',
     +   'id does not satisfy 0.le.id.lt.k', 2, 1)
      return
  115 continue
      call xermsg ('slatec', 'dpfqad', 'lxi does not satisfy lxi.ge.1',
     +   2, 1)
      return
      end
