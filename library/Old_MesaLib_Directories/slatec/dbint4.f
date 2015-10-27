*deck dbint4
      subroutine dbint4 (x, y, ndata, ibcl, ibcr, fbcl, fbcr, kntopt, t,
     +   bcoef, n, k, w)
c***begin prologue  dbint4
c***purpose  compute the b-representation of a cubic spline
c            which interpolates given data.
c***library   slatec
c***category  e1a
c***type      double precision (bint4-s, dbint4-d)
c***keywords  b-spline, cubic splines, data fitting, interpolation
c***author  amos, d. e., (snla)
c***description
c
c     abstract    **** a double precision routine ****
c
c         dbint4 computes the b representation (t,bcoef,n,k) of a
c         cubic spline (k=4) which interpolates data (x(i),y(i)),
c         i=1,ndata.  parameters ibcl, ibcr, fbcl, fbcr allow the
c         specification of the spline first or second derivative at
c         both x(1) and x(ndata).  when this data is not specified
c         by the problem, it is common practice to use a natural
c         spline by setting second derivatives at x(1) and x(ndata)
c         to zero (ibcl=ibcr=2,fbcl=fbcr=0.0).  the spline is defined
c         on t(4) .le. x .le. t(n+1) with (ordered) interior knots at
c         x(i) values where n=ndata+2.  the knots t(1),t(2),t(3) lie to
c         the left of t(4)=x(1) and the knots t(n+2), t(n+3), t(n+4)
c         lie to the right of t(n+1)=x(ndata) in increasing order.  if
c         no extrapolation outside (x(1),x(ndata)) is anticipated, the
c         knots t(1)=t(2)=t(3)=t(4)=x(1) and t(n+2)=t(n+3)=t(n+4)=
c         t(n+1)=x(ndata) can be specified by kntopt=1.  kntopt=2
c         selects a knot placement for t(1), t(2), t(3) to make the
c         first 7 knots symmetric about t(4)=x(1) and similarly for
c         t(n+2), t(n+3), t(n+4) about t(n+1)=x(ndata).  kntopt=3
c         allows the user to make his own selection, in increasing
c         order, for t(1), t(2), t(3) to the left of x(1) and t(n+2),
c         t(n+3), t(n+4) to the right of x(ndata) in the work array
c         w(1) through w(6).  in any case, the interpolation on
c         t(4) .le. x .le. t(n+1) by using function dbvalu is unique
c         for given boundary conditions.
c
c     description of arguments
c
c         input      x,y,fbcl,fbcr,w are double precision
c           x      - x vector of abscissae of length ndata, distinct
c                    and in increasing order
c           y      - y vector of ordinates of length ndata
c           ndata  - number of data points, ndata .ge. 2
c           ibcl   - selection parameter for left boundary condition
c                    ibcl = 1 constrain the first derivative at
c                             x(1) to fbcl
c                         = 2 constrain the second derivative at
c                             x(1) to fbcl
c           ibcr   - selection parameter for right boundary condition
c                    ibcr = 1 constrain first derivative at
c                             x(ndata) to fbcr
c                    ibcr = 2 constrain second derivative at
c                             x(ndata) to fbcr
c           fbcl   - left boundary values governed by ibcl
c           fbcr   - right boundary values governed by ibcr
c           kntopt - knot selection parameter
c                    kntopt = 1 sets knot multiplicity at t(4) and
c                               t(n+1) to 4
c                           = 2 sets a symmetric placement of knots
c                               about t(4) and t(n+1)
c                           = 3 sets t(i)=w(i) and t(n+1+i)=w(3+i),i=1,3
c                               where w(i),i=1,6 is supplied by the user
c           w      - work array of dimension at least 5*(ndata+2)
c                    if kntopt=3, then w(1),w(2),w(3) are knot values to
c                    the left of x(1) and w(4),w(5),w(6) are knot
c                    values to the right of x(ndata) in increasing
c                    order to be supplied by the user
c
c         output     t,bcoef are double precision
c           t      - knot array of length n+4
c           bcoef  - b spline coefficient array of length n
c           n      - number of coefficients, n=ndata+2
c           k      - order of spline, k=4
c
c     error conditions
c         improper  input is a fatal error
c         singular system of equations is a fatal error
c
c***references  d. e. amos, computation with splines and b-splines,
c                 report sand78-1968, sandia laboratories, march 1979.
c               carl de boor, package for calculating with b-splines,
c                 siam journal on numerical analysis 14, 3 (june 1977),
c                 pp. 441-472.
c               carl de boor, a practical guide to splines, applied
c                 mathematics series 27, springer-verlag, new york,
c                 1978.
c***routines called  d1mach, dbnfac, dbnslv, dbspvd, xermsg
c***revision history  (yymmdd)
c   800901  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  dbint4
c
      integer i, ibcl, ibcr, iflag, ilb, ileft, it, iub, iw, iwp, j,
     1 jw, k, kntopt, n, ndata, ndm, np, nwrow
      double precision bcoef,fbcl,fbcr,t,tol,txn,tx1,vnikx,w,wdtol,
     1 work,x,xl,y
      double precision d1mach
      dimension x(*), y(*), t(*), bcoef(*), w(5,*), vnikx(4,4), work(15)
c***first executable statement  dbint4
      wdtol = d1mach(4)
      tol = sqrt(wdtol)
      if (ndata.lt.2) go to 200
      ndm = ndata - 1
      do 10 i=1,ndm
        if (x(i).ge.x(i+1)) go to 210
   10 continue
      if (ibcl.lt.1 .or. ibcl.gt.2) go to 220
      if (ibcr.lt.1 .or. ibcr.gt.2) go to 230
      if (kntopt.lt.1 .or. kntopt.gt.3) go to 240
      k = 4
      n = ndata + 2
      np = n + 1
      do 20 i=1,ndata
        t(i+3) = x(i)
   20 continue
      go to (30, 50, 90), kntopt
c     set up knot array with multiplicity 4 at x(1) and x(ndata)
   30 continue
      do 40 i=1,3
        t(4-i) = x(1)
        t(np+i) = x(ndata)
   40 continue
      go to 110
c     set up knot array with symmetric placement about end points
   50 continue
      if (ndata.gt.3) go to 70
      xl = (x(ndata)-x(1))/3.0d0
      do 60 i=1,3
        t(4-i) = t(5-i) - xl
        t(np+i) = t(np+i-1) + xl
   60 continue
      go to 110
   70 continue
      tx1 = x(1) + x(1)
      txn = x(ndata) + x(ndata)
      do 80 i=1,3
        t(4-i) = tx1 - x(i+1)
        t(np+i) = txn - x(ndata-i)
   80 continue
      go to 110
c     set up knot array less than x(1) and greater than x(ndata) to be
c     supplied by user in work locations w(1) through w(6) when kntopt=3
   90 continue
      do 100 i=1,3
        t(4-i) = w(4-i,1)
        jw = max(1,i-1)
        iw = mod(i+2,5)+1
        t(np+i) = w(iw,jw)
        if (t(4-i).gt.t(5-i)) go to 250
        if (t(np+i).lt.t(np+i-1)) go to 250
  100 continue
  110 continue
c
      do 130 i=1,5
        do 120 j=1,n
          w(i,j) = 0.0d0
  120   continue
  130 continue
c     set up left interpolation point and left boundary condition for
c     right limits
      it = ibcl + 1
      call dbspvd(t, k, it, x(1), k, 4, vnikx, work)
      iw = 0
      if (abs(vnikx(3,1)).lt.tol) iw = 1
      do 140 j=1,3
        w(j+1,4-j) = vnikx(4-j,it)
        w(j,4-j) = vnikx(4-j,1)
  140 continue
      bcoef(1) = y(1)
      bcoef(2) = fbcl
c     set up interpolation equations for points i=2 to i=ndata-1
      ileft = 4
      if (ndm.lt.2) go to 170
      do 160 i=2,ndm
        ileft = ileft + 1
        call dbspvd(t, k, 1, x(i), ileft, 4, vnikx, work)
        do 150 j=1,3
          w(j+1,3+i-j) = vnikx(4-j,1)
  150   continue
        bcoef(i+1) = y(i)
  160 continue
c     set up right interpolation point and right boundary condition for
c     left limits(ileft is associated with t(n)=x(ndata-1))
  170 continue
      it = ibcr + 1
      call dbspvd(t, k, it, x(ndata), ileft, 4, vnikx, work)
      jw = 0
      if (abs(vnikx(2,1)).lt.tol) jw = 1
      do 180 j=1,3
        w(j+1,3+ndata-j) = vnikx(5-j,it)
        w(j+2,3+ndata-j) = vnikx(5-j,1)
  180 continue
      bcoef(n-1) = fbcr
      bcoef(n) = y(ndata)
c     solve system of equations
      ilb = 2 - jw
      iub = 2 - iw
      nwrow = 5
      iwp = iw + 1
      call dbnfac(w(iwp,1), nwrow, n, ilb, iub, iflag)
      if (iflag.eq.2) go to 190
      call dbnslv(w(iwp,1), nwrow, n, ilb, iub, bcoef)
      return
c
c
  190 continue
      call xermsg ('slatec', 'dbint4',
     +   'the system of equations is singular', 2, 1)
      return
  200 continue
      call xermsg ('slatec', 'dbint4', 'ndata is less than 2', 2, 1)
      return
  210 continue
      call xermsg ('slatec', 'dbint4',
     +   'x values are not distinct or not ordered', 2, 1)
      return
  220 continue
      call xermsg ('slatec', 'dbint4', 'ibcl is not 1 or 2', 2, 1)
      return
  230 continue
      call xermsg ('slatec', 'dbint4', 'ibcr is not 1 or 2', 2, 1)
      return
  240 continue
      call xermsg ('slatec', 'dbint4', 'kntopt is not 1, 2, or 3', 2,
     +   1)
      return
  250 continue
      call xermsg ('slatec', 'dbint4',
     +   'knot input through w array is not ordered properly', 2, 1)
      return
      end
