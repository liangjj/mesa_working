*deck bsgq8
      subroutine bsgq8 (fun, xt, bc, n, kk, id, a, b, inbv, err, ans,
     +   ierr, work)
c***begin prologue  bsgq8
c***subsidiary
c***purpose  subsidiary to bfqad
c***library   slatec
c***type      single precision (bsgq8-s, dbsgq8-d)
c***author  jones, r. e., (snla)
c***description
c
c     abstract
c        bsgq8, a modification of gaus8, integrates the
c        product of fun(x) by the id-th derivative of a spline
c        bvalu(xt,bc,n,kk,id,x,inbv,work)  between limits a and b.
c
c     description of arguments
c
c        input--
c        fun - name of external function of one argument which
c              multiplies bvalu.
c        xt  - knot array for bvalu
c        bc  - b-coefficient array for bvalu
c        n   - number of b-coefficients for bvalu
c        kk  - order of the spline, kk.ge.1
c        id  - order of the spline derivative, 0.le.id.le.kk-1
c        a   - lower limit of integral
c        b   - upper limit of integral (may be less than a)
c        inbv- initialization parameter for bvalu
c        err - is a requested pseudorelative error tolerance.  normally
c              pick a value of abs(err).lt.1e-3.  ans will normally
c              have no more error than abs(err) times the integral of
c              the absolute value of fun(x)*bvalu(xt,bc,n,kk,x,id,
c              inbv,work).
c
c
c        output--
c        err - will be an estimate of the absolute error in ans if the
c              input value of err was negative.  (err is unchanged if
c              the input value of err was nonnegative.)  the estimated
c              error is solely for information to the user and should
c              not be used as a correction to the computed integral.
c        ans - computed value of integral
c        ierr- a status code
c            --normal codes
c               1 ans most likely meets requested error tolerance,
c                 or a=b.
c              -1 a and b are too nearly equal to allow normal
c                 integration.  ans is set to zero.
c            --abnormal code
c               2 ans probably does not meet requested error tolerance.
c        work- work vector of length 3*k for bvalu
c
c***see also  bfqad
c***routines called  bvalu, i1mach, r1mach, xermsg
c***revision history  (yymmdd)
c   800901  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c   900328  added type section.  (wrb)
c   910408  updated the author section.  (wrb)
c***end prologue  bsgq8
c
      integer id, ierr, inbv, k, kk, kml, kmx, l, lmn, lmx, lr, mxl,
     1 n, nbits, nib, nlmn, nlmx
      integer i1mach
      real a, aa, ae, anib, ans, area, b, bc, c, ce, ee, ef, eps, err,
     1 est,gl,glr,gr,hh,sq2,tol,vl,vr,work,w1, w2, w3, w4, xt, x1,
     2 x2, x3, x4, x, h
      real r1mach, bvalu, g8, fun
      dimension xt(*), bc(*)
      dimension aa(30), hh(30), lr(30), vl(30), gr(30)
      save x1, x2, x3, x4, w1, w2, w3, w4, sq2, nlmn, kmx, kml
      data x1, x2, x3, x4/
     1     1.83434642495649805e-01,     5.25532409916328986e-01,
     2     7.96666477413626740e-01,     9.60289856497536232e-01/
      data w1, w2, w3, w4/
     1     3.62683783378361983e-01,     3.13706645877887287e-01,
     2     2.22381034453374471e-01,     1.01228536290376259e-01/
      data sq2/1.41421356e0/
      data nlmn/1/,kmx/5000/,kml/6/
      g8(x,h)=h*((w1*(fun(x-x1*h)*bvalu(xt,bc,n,kk,id,x-x1*h,inbv,work)+
     1                fun(x+x1*h)*bvalu(xt,bc,n,kk,id,x+x1*h,inbv,work))
     2           +w2*(fun(x-x2*h)*bvalu(xt,bc,n,kk,id,x-x2*h,inbv,work)+
     3              fun(x+x2*h)*bvalu(xt,bc,n,kk,id,x+x2*h,inbv,work)))
     4          +(w3*(fun(x-x3*h)*bvalu(xt,bc,n,kk,id,x-x3*h,inbv,work)+
     5                fun(x+x3*h)*bvalu(xt,bc,n,kk,id,x+x3*h,inbv,work))
     6           +w4*(fun(x-x4*h)*bvalu(xt,bc,n,kk,id,x-x4*h,inbv,work)+
     7             fun(x+x4*h)*bvalu(xt,bc,n,kk,id,x+x4*h,inbv,work))))
c
c     initialize
c
c***first executable statement  bsgq8
      k = i1mach(11)
      anib = r1mach(5)*k/0.30102000e0
      nbits = int(anib)
      nlmx = (nbits*5)/8
      ans = 0.0e0
      ierr = 1
      ce = 0.0e0
      if (a.eq.b) go to 140
      lmx = nlmx
      lmn = nlmn
      if (b.eq.0.0e0) go to 10
      if (sign(1.0e0,b)*a.le.0.0e0) go to 10
      c = abs(1.0e0-a/b)
      if (c.gt.0.1e0) go to 10
      if (c.le.0.0e0) go to 140
      anib = 0.5e0 - log(c)/0.69314718e0
      nib = int(anib)
      lmx = min(nlmx,nbits-nib-7)
      if (lmx.lt.1) go to 130
      lmn = min(lmn,lmx)
   10 tol = max(abs(err),2.0e0**(5-nbits))/2.0e0
      if (err.eq.0.0e0) tol = sqrt(r1mach(4))
      eps = tol
      hh(1) = (b-a)/4.0e0
      aa(1) = a
      lr(1) = 1
      l = 1
      est = g8(aa(l)+2.0e0*hh(l),2.0e0*hh(l))
      k = 8
      area = abs(est)
      ef = 0.5e0
      mxl = 0
c
c     compute refined estimates, estimate the error, etc.
c
   20 gl = g8(aa(l)+hh(l),hh(l))
      gr(l) = g8(aa(l)+3.0e0*hh(l),hh(l))
      k = k + 16
      area = area + (abs(gl)+abs(gr(l))-abs(est))
      glr = gl + gr(l)
      ee = abs(est-glr)*ef
      ae = max(eps*area,tol*abs(glr))
      if (ee-ae) 40, 40, 50
   30 mxl = 1
   40 ce = ce + (est-glr)
      if (lr(l)) 60, 60, 80
c
c     consider the left half of this level
c
   50 if (k.gt.kmx) lmx = kml
      if (l.ge.lmx) go to 30
      l = l + 1
      eps = eps*0.5e0
      ef = ef/sq2
      hh(l) = hh(l-1)*0.5e0
      lr(l) = -1
      aa(l) = aa(l-1)
      est = gl
      go to 20
c
c     proceed to right half at this level
c
   60 vl(l) = glr
   70 est = gr(l-1)
      lr(l) = 1
      aa(l) = aa(l) + 4.0e0*hh(l)
      go to 20
c
c     return one level
c
   80 vr = glr
   90 if (l.le.1) go to 120
      l = l - 1
      eps = eps*2.0e0
      ef = ef*sq2
      if (lr(l)) 100, 100, 110
  100 vl(l) = vl(l+1) + vr
      go to 70
  110 vr = vl(l+1) + vr
      go to 90
c
c      exit
c
  120 ans = vr
      if ((mxl.eq.0) .or. (abs(ce).le.2.0e0*tol*area)) go to 140
      ierr = 2
      call xermsg ('slatec', 'bsgq8',
     +   'ans is probably insufficiently accurate.', 3, 1)
      go to 140
  130 ierr = -1
      call xermsg ('slatec', 'bsgq8',
     +   'a and b are too nearly equal to allow normal integration. ' //
     +   ' ans is set to zero and ierr to -1.', 1, -1)
  140 continue
      if (err.lt.0.0e0) err = ce
      return
      end
