*deck dgamrn
      double precision function dgamrn (x)
c***begin prologue  dgamrn
c***subsidiary
c***purpose  subsidiary to dbskin
c***library   slatec
c***type      double precision (gamrn-s, dgamrn-d)
c***author  amos, d. e., (snla)
c***description
c
c     abstract   * a double precision routine *
c         dgamrn computes the gamma function ratio gamma(x)/gamma(x+0.5)
c         for real x.gt.0. if x.ge.xmin, an asymptotic expansion is
c         evaluated. if x.lt.xmin, an integer is added to x to form a
c         new value of x.ge.xmin and the asymptotic expansion is eval-
c         uated for this new value of x. successive application of the
c         recurrence relation
c
c                      w(x)=w(x+1)*(1+0.5/x)
c
c         reduces the argument to its original value. xmin and comp-
c         utational tolerances are computed as a function of the number
c         of digits carried in a word by calls to i1mach and d1mach.
c         however, the computational accuracy is limited to the max-
c         imum of unit roundoff (=d1mach(4)) and 1.0d-18 since critical
c         constants are given to only 18 digits.
c
c         input      x is double precision
c           x      - argument, x.gt.0.0d0
c
c         output      dgamrn is double precision
c           dgamrn  - ratio  gamma(x)/gamma(x+0.5)
c
c***see also  dbskin
c***references  y. l. luke, the special functions and their
c                 approximations, vol. 1, math in sci. and
c                 eng. series 53, academic press, new york, 1969,
c                 pp. 34-35.
c***routines called  d1mach, i1mach
c***revision history  (yymmdd)
c   820601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890911  removed unnecessary intrinsics.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910722  updated author section.  (als)
c   920520  added references section.  (wrb)
c***end prologue  dgamrn
      integer i, i1m11, k, mx, nx
      integer i1mach
      double precision fln, gr, rln, s, tol, trm, x, xdmy, xinc, xm,
     * xmin, xp, xsq
      double precision d1mach
      dimension gr(12)
      save gr
c
      data gr(1), gr(2), gr(3), gr(4), gr(5), gr(6), gr(7), gr(8),
     * gr(9), gr(10), gr(11), gr(12) /1.00000000000000000d+00,
     * -1.56250000000000000d-02,2.56347656250000000d-03,
     * -1.27983093261718750d-03,1.34351104497909546d-03,
     * -2.43289663922041655d-03,6.75423753364157164d-03,
     * -2.66369606131178216d-02,1.41527455519564332d-01,
     * -9.74384543032201613d-01,8.43686251229783675d+00,
     * -8.97258321640552515d+01/
c
c***first executable statement  dgamrn
      nx = int(x)
      tol = max(d1mach(4),1.0d-18)
      i1m11 = i1mach(14)
      rln = d1mach(5)*i1m11
      fln = min(rln,20.0d0)
      fln = max(fln,3.0d0)
      fln = fln - 3.0d0
      xm = 2.0d0 + fln*(0.2366d0+0.01723d0*fln)
      mx = int(xm) + 1
      xmin = mx
      xdmy = x - 0.25d0
      xinc = 0.0d0
      if (x.ge.xmin) go to 10
      xinc = xmin - nx
      xdmy = xdmy + xinc
   10 continue
      s = 1.0d0
      if (xdmy*tol.gt.1.0d0) go to 30
      xsq = 1.0d0/(xdmy*xdmy)
      xp = xsq
      do 20 k=2,12
        trm = gr(k)*xp
        if (abs(trm).lt.tol) go to 30
        s = s + trm
        xp = xp*xsq
   20 continue
   30 continue
      s = s/sqrt(xdmy)
      if (xinc.ne.0.0d0) go to 40
      dgamrn = s
      return
   40 continue
      nx = int(xinc)
      xp = 0.0d0
      do 50 i=1,nx
        s = s*(1.0d0+0.5d0/(x+xp))
        xp = xp + 1.0d0
   50 continue
      dgamrn = s
      return
      end
