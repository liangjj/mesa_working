*deck bsqad
      subroutine bsqad (t, bcoef, n, k, x1, x2, bquad, work)
c***begin prologue  bsqad
c***purpose  compute the integral of a k-th order b-spline using the
c            b-representation.
c***library   slatec
c***category  h2a2a1, e3, k6
c***type      single precision (bsqad-s, dbsqad-d)
c***keywords  integral of b-splines, quadrature
c***author  amos, d. e., (snla)
c***description
c
c     abstract
c         bsqad computes the integral on (x1,x2) of a k-th order
c         b-spline using the b-representation (t,bcoef,n,k).  orders
c         k as high as 20 are permitted by applying a 2, 6, or 10
c         point gauss formula on subintervals of (x1,x2) which are
c         formed by included (distinct) knots.
c
c         if orders k greater than 20 are needed, use bfqad with
c         f(x) = 1.
c
c     description of arguments
c         input
c           t      - knot array of length n+k
c           bcoef  - b-spline coefficient array of length n
c           n      - length of coefficient array
c           k      - order of b-spline, 1 .le. k .le. 20
c           x1,x2  - end points of quadrature interval in
c                    t(k) .le. x .le. t(n+1)
c
c         output
c           bquad  - integral of the b-spline over (x1,x2)
c           work   - work vector of length 3*k
c
c     error conditions
c         improper input is a fatal error
c
c***references  d. e. amos, quadrature subroutines for splines and
c                 b-splines, report sand79-1825, sandia laboratories,
c                 december 1979.
c***routines called  bvalu, intrv, xermsg
c***revision history  (yymmdd)
c   800901  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  bsqad
c
      integer i,il1,il2,ilo,inbv, jf,k,left,m,mf,mflag,n, npk, np1
      real a, aa, b, bb, bcoef, bma, bpa, bquad, c1, gpts, gwts, gx, q,
     1 sum, t, ta, tb, work, x1, x2, y1, y2
      real bvalu
      dimension t(*), bcoef(*), gpts(9), gwts(9), sum(5), work(*)
c
      save gpts, gwts
      data gpts(1), gpts(2), gpts(3), gpts(4), gpts(5), gpts(6),
     1     gpts(7), gpts(8), gpts(9)/
     2     5.77350269189625764e-01,     2.38619186083196909e-01,
     3     6.61209386466264514e-01,     9.32469514203152028e-01,
     4     1.48874338981631211e-01,     4.33395394129247191e-01,
     5     6.79409568299024406e-01,     8.65063366688984511e-01,
     6     9.73906528517171720e-01/
      data gwts(1), gwts(2), gwts(3), gwts(4), gwts(5), gwts(6),
     1     gwts(7), gwts(8), gwts(9)/
     2     1.00000000000000000e+00,     4.67913934572691047e-01,
     3     3.60761573048138608e-01,     1.71324492379170345e-01,
     4     2.95524224714752870e-01,     2.69266719309996355e-01,
     5     2.19086362515982044e-01,     1.49451349150580593e-01,
     6     6.66713443086881376e-02/
c
c***first executable statement  bsqad
      bquad = 0.0e0
      if(k.lt.1 .or. k.gt.20) go to 65
      if(n.lt.k) go to 70
      aa = min(x1,x2)
      bb = max(x1,x2)
      if (aa.lt.t(k)) go to 60
      np1 = n + 1
      if (bb.gt.t(np1)) go to 60
      if (aa.eq.bb) return
      npk = n + k
c     selection of 2, 6, or 10 point gauss formula
      jf = 0
      mf = 1
      if (k.le.4) go to 10
      jf = 1
      mf = 3
      if (k.le.12) go to 10
      jf = 4
      mf = 5
   10 continue
c
      do 20 i=1,mf
        sum(i) = 0.0e0
   20 continue
      ilo = 1
      inbv = 1
      call intrv(t, npk, aa, ilo, il1, mflag)
      call intrv(t, npk, bb, ilo, il2, mflag)
      if (il2.ge.np1) il2 = n
      do 40 left=il1,il2
        ta = t(left)
        tb = t(left+1)
        if (ta.eq.tb) go to 40
        a = max(aa,ta)
        b = min(bb,tb)
        bma = 0.5e0*(b-a)
        bpa = 0.5e0*(b+a)
        do 30 m=1,mf
          c1 = bma*gpts(jf+m)
          gx = -c1 + bpa
          y2 = bvalu(t,bcoef,n,k,0,gx,inbv,work)
          gx = c1 + bpa
          y1 = bvalu(t,bcoef,n,k,0,gx,inbv,work)
          sum(m) = sum(m) + (y1+y2)*bma
   30   continue
   40 continue
      q = 0.0e0
      do 50 m=1,mf
        q = q + gwts(jf+m)*sum(m)
   50 continue
      if (x1.gt.x2) q = -q
      bquad = q
      return
c
c
   60 continue
      call xermsg ('slatec', 'bsqad',
     +   'x1 or x2 or both do not satisfy t(k).le.x.le.t(n+1)', 2, 1)
      return
   65 continue
      call xermsg ('slatec', 'bsqad', 'k does not satisfy 1.le.k.le.20'
     +   , 2, 1)
      return
   70 continue
      call xermsg ('slatec', 'bsqad', 'n does not satisfy n.ge.k', 2,
     +   1)
      return
      end
