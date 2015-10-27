*deck dbsgq8
      subroutine dbsgq8 (fun, xt, bc, n, kk, id, a, b, inbv, err, ans,
     +   ierr, work)
c***begin prologue  dbsgq8
c***subsidiary
c***purpose  subsidiary to dbfqad
c***library   slatec
c***type      double precision (bsgq8-s, dbsgq8-d)
c***author  jones, r. e., (snla)
c***description
c
c     abstract    **** a double precision routine ****
c
c        dbsgq8, a modification of gaus8, integrates the
c        product of fun(x) by the id-th derivative of a spline
c        dbvalu(xt,bc,n,kk,id,x,inbv,work)  between limits a and b.
c
c     description of arguments
c
c        input-- fun,xt,bc,a,b,err are double precision
c        fun - name of external function of one argument which
c              multiplies dbvalu.
c        xt  - knot array for dbvalu
c        bc  - b-coefficient array for dbvalu
c        n   - number of b-coefficients for dbvalu
c        kk  - order of the spline, kk.ge.1
c        id  - order of the spline derivative, 0.le.id.le.kk-1
c        a   - lower limit of integral
c        b   - upper limit of integral (may be less than a)
c        inbv- initialization parameter for dbvalu
c        err - is a requested pseudorelative error tolerance.  normally
c              pick a value of abs(err).lt.1d-3.  ans will normally
c              have no more error than abs(err) times the integral of
c              the absolute value of fun(x)*dbvalu(xt,bc,n,kk,x,id,
c              inbv,work).
c
c
c        output-- err,ans,work are double precision
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
c        work- work vector of length 3*k for dbvalu
c
c***see also  dbfqad
c***routines called  d1mach, dbvalu, i1mach, xermsg
c***revision history  (yymmdd)
c   800901  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890911  removed unnecessary intrinsics.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c   900328  added type section.  (wrb)
c   910408  updated the author section.  (wrb)
c***end prologue  dbsgq8
c
      integer id, ierr, inbv, k, kk, kml, kmx, l, lmn, lmx, lr, mxl,
     1 n, nbits, nib, nlmn, nlmx
      integer i1mach
      double precision a,aa,ae,anib,ans,area,b,bc,c,ce,ee,ef,eps,err,
     1 est,gl,glr,gr,hh,sq2,tol,vl,vr,work,w1, w2, w3, w4, xt, x1,
     2 x2, x3, x4, x, h
      double precision d1mach, dbvalu, g8, fun
      dimension xt(*), bc(*), work(*)
      dimension aa(60), hh(60), lr(60), vl(60), gr(60)
      save x1, x2, x3, x4, w1, w2, w3, w4, sq2, nlmn, kmx, kml
      data x1, x2, x3, x4/
     1     1.83434642495649805d-01,     5.25532409916328986d-01,
     2     7.96666477413626740d-01,     9.60289856497536232d-01/
      data w1, w2, w3, w4/
     1     3.62683783378361983d-01,     3.13706645877887287d-01,
     2     2.22381034453374471d-01,     1.01228536290376259d-01/
      data sq2/1.41421356d0/
      data nlmn/1/,kmx/5000/,kml/6/
      g8(x,h)=h*((w1*(fun(x-x1*h)*dbvalu(xt,bc,n,kk,id,x-x1*h,inbv,work)
     1            + fun(x+x1*h)*dbvalu(xt,bc,n,kk,id,x+x1*h,inbv,work))
     2          +w2*(fun(x-x2*h)*dbvalu(xt,bc,n,kk,id,x-x2*h,inbv,work)+
     3              fun(x+x2*h)*dbvalu(xt,bc,n,kk,id,x+x2*h,inbv,work)))
     4         +(w3*(fun(x-x3*h)*dbvalu(xt,bc,n,kk,id,x-x3*h,inbv,work)+
     5               fun(x+x3*h)*dbvalu(xt,bc,n,kk,id,x+x3*h,inbv,work))
     6          +w4*(fun(x-x4*h)*dbvalu(xt,bc,n,kk,id,x-x4*h,inbv,work)+
     7             fun(x+x4*h)*dbvalu(xt,bc,n,kk,id,x+x4*h,inbv,work))))
c
c     initialize
c
c***first executable statement  dbsgq8
      k = i1mach(14)
      anib = d1mach(5)*k/0.30102000d0
      nbits = int(anib)
      nlmx = min((nbits*5)/8,60)
      ans = 0.0d0
      ierr = 1
      ce = 0.0d0
      if (a.eq.b) go to 140
      lmx = nlmx
      lmn = nlmn
      if (b.eq.0.0d0) go to 10
      if (sign(1.0d0,b)*a.le.0.0d0) go to 10
      c = abs(1.0d0-a/b)
      if (c.gt.0.1d0) go to 10
      if (c.le.0.0d0) go to 140
      anib = 0.5d0 - log(c)/0.69314718d0
      nib = int(anib)
      lmx = min(nlmx,nbits-nib-7)
      if (lmx.lt.1) go to 130
      lmn = min(lmn,lmx)
   10 tol = max(abs(err),2.0d0**(5-nbits))/2.0d0
      if (err.eq.0.0d0) tol = sqrt(d1mach(4))
      eps = tol
      hh(1) = (b-a)/4.0d0
      aa(1) = a
      lr(1) = 1
      l = 1
      est = g8(aa(l)+2.0d0*hh(l),2.0d0*hh(l))
      k = 8
      area = abs(est)
      ef = 0.5d0
      mxl = 0
c
c     compute refined estimates, estimate the error, etc.
c
   20 gl = g8(aa(l)+hh(l),hh(l))
      gr(l) = g8(aa(l)+3.0d0*hh(l),hh(l))
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
      eps = eps*0.5d0
      ef = ef/sq2
      hh(l) = hh(l-1)*0.5d0
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
      aa(l) = aa(l) + 4.0d0*hh(l)
      go to 20
c
c     return one level
c
   80 vr = glr
   90 if (l.le.1) go to 120
      l = l - 1
      eps = eps*2.0d0
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
      if ((mxl.eq.0) .or. (abs(ce).le.2.0d0*tol*area)) go to 140
      ierr = 2
      call xermsg ('slatec', 'dbsgq8',
     +   'ans is probably insufficiently accurate.', 3, 1)
      go to 140
  130 ierr = -1
      call xermsg ('slatec', 'dbsgq8',
     +   'a and b are too nearly equal to allow normal integration. ' //
     +   ' ans is set to zero and ierr to -1.', 1, -1)
  140 continue
      if (err.lt.0.0d0) err = ce
      return
      end
