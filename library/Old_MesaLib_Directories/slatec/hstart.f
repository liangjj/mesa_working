*deck hstart
      subroutine hstart (f, neq, a, b, y, yprime, etol, morder, small,
     +   big, spy, pv, yp, sf, rpar, ipar, h)
c***begin prologue  hstart
c***subsidiary
c***purpose  subsidiary to deabm, debdf and derkf
c***library   slatec
c***type      single precision (hstart-s, dhstrt-d)
c***author  watts, h. a., (snla)
c***description
c
c   hstart computes a starting step size to be used in solving initial
c   value problems in ordinary differential equations.
c **********************************************************************
c  abstract
c
c     subroutine hstart computes a starting step size to be used by an
c     initial value method in solving ordinary differential equations.
c     it is based on an estimate of the local lipschitz constant for the
c     differential equation (lower bound on a norm of the jacobian),
c     a bound on the differential equation (first derivative), and
c     a bound on the partial derivative of the equation with respect to
c     the independent variable.
c     (all approximated near the initial point a.)
c
c     subroutine hstart uses a function subprogram hvnrm for computing
c     a vector norm.  the maximum norm is presently utilized though it
c     can easily be replaced by any other vector norm.  it is presumed
c     that any replacement norm routine would be carefully coded to
c     prevent unnecessary underflows or overflows from occurring, and
c     also, would not alter the vector or number of components.
c
c **********************************************************************
c  on input you must provide the following
c
c      f -- this is a subroutine of the form
c                               f(x,u,uprime,rpar,ipar)
c             which defines the system of first order differential
c             equations to be solved.  for the given values of x and the
c             vector  u(*)=(u(1),u(2),...,u(neq)) , the subroutine must
c             evaluate the neq components of the system of differential
c             equations  du/dx=f(x,u)  and store the derivatives in the
c             array uprime(*), that is,  uprime(i) = * du(i)/dx *  for
c             equations i=1,...,neq.
c
c             subroutine f must not alter x or u(*).  you must declare
c             the name f in an external statement in your program that
c             calls hstart.  you must dimension u and uprime in f.
c
c             rpar and ipar are real and integer parameter arrays which
c             you can use for communication between your program and
c             subroutine f.  they are not used or altered by hstart.  if
c             you do not need rpar or ipar, ignore these parameters by
c             treating them as dummy arguments.  if you do choose to use
c             them, dimension them in your program and in f as arrays
c             of appropriate length.
c
c      neq -- this is the number of (first order) differential equations
c             to be integrated.
c
c      a -- this is the initial point of integration.
c
c      b -- this is a value of the independent variable used to define
c             the direction of integration.  a reasonable choice is to
c             set  b  to the first point at which a solution is desired.
c             you can also use  b, if necessary, to restrict the length
c             of the first integration step because the algorithm will
c             not compute a starting step length which is bigger than
c             abs(b-a), unless  b  has been chosen too close to  a.
c             (it is presumed that hstart has been called with  b
c             different from  a  on the machine being used.  also see
c             the discussion about the parameter  small.)
c
c      y(*) -- this is the vector of initial values of the neq solution
c             components at the initial point  a.
c
c      yprime(*) -- this is the vector of derivatives of the neq
c             solution components at the initial point  a.
c             (defined by the differential equations in subroutine f)
c
c      etol -- this is the vector of error tolerances corresponding to
c             the neq solution components.  it is assumed that all
c             elements are positive.  following the first integration
c             step, the tolerances are expected to be used by the
c             integrator in an error test which roughly requires that
c                        abs(local error) .le. etol
c             for each vector component.
c
c      morder -- this is the order of the formula which will be used by
c             the initial value method for taking the first integration
c             step.
c
c      small -- this is a small positive machine dependent constant
c             which is used for protecting against computations with
c             numbers which are too small relative to the precision of
c             floating point arithmetic.  small  should be set to
c             (approximately) the smallest positive real number such
c             that  (1.+small) .gt. 1.  on the machine being used. the
c             quantity  small**(3/8)  is used in computing increments of
c             variables for approximating derivatives by differences.
c             also the algorithm will not compute a starting step length
c             which is smaller than  100*small*abs(a).
c
c      big -- this is a large positive machine dependent constant which
c             is used for preventing machine overflows.  a reasonable
c             choice is to set big to (approximately) the square root of
c             the largest real number which can be held in the machine.
c
c      spy(*),pv(*),yp(*),sf(*) -- these are real work arrays of length
c             neq which provide the routine with needed storage space.
c
c      rpar,ipar -- these are parameter arrays, of real and integer
c             type, respectively, which can be used for communication
c             between your program and the f subroutine.  they are not
c             used or altered by hstart.
c
c **********************************************************************
c  on output  (after the return from hstart),
c
c      h -- is an appropriate starting step size to be attempted by the
c             differential equation method.
c
c           all parameters in the call list remain unchanged except for
c           the working arrays spy(*),pv(*),yp(*) and sf(*).
c
c **********************************************************************
c
c***see also  deabm, debdf, derkf
c***routines called  hvnrm
c***revision history  (yymmdd)
c   800501  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891024  changed references from vnorm to hvnrm.  (wrb)
c   891024  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910722  updated author section.  (als)
c***end prologue  hstart
c
      dimension y(*),yprime(*),etol(*),spy(*),pv(*),yp(*),sf(*),
     1   rpar(*),ipar(*)
      external f
c
c.......................................................................
c
c***first executable statement  hstart
      dx = b - a
      absdx = abs(dx)
      relper = small**0.375
      ynorm = hvnrm(y,neq)
c
c.......................................................................
c
c     compute a weighted approximate bound (dfdxb) on the partial
c     derivative of the equation with respect to the
c     independent variable. protect against an overflow. also
c     compute a weighted bound (fbnd) on the first derivative locally.
c
      da = sign(max(min(relper*abs(a),absdx),100.*small*abs(a)),dx)
      if (da .eq. 0.) da = relper*dx
      call f(a+da,y,sf,rpar,ipar)
c
      if (morder .eq. 1) go to 20
      power = 2./(morder+1)
      do 10 j=1,neq
        wtj = etol(j)**power
        spy(j) = sf(j)/wtj
        yp(j) = yprime(j)/wtj
   10   pv(j) = spy(j) - yp(j)
      go to 40
c
   20 do 30 j=1,neq
        spy(j) = sf(j)/etol(j)
        yp(j) = yprime(j)/etol(j)
   30   pv(j) = spy(j) - yp(j)
c
   40 delf = hvnrm(pv,neq)
      dfdxb = big
      if (delf .lt. big*abs(da)) dfdxb = delf/abs(da)
      ypnorm = hvnrm(yp,neq)
      fbnd = max(hvnrm(spy,neq),ypnorm)
c
c.......................................................................
c
c     compute an estimate (dfdub) of the local lipschitz constant for
c     the system of differential equations. this also represents an
c     estimate of the norm of the jacobian locally.
c     three iterations (two when neq=1) are used to estimate the
c     lipschitz constant by numerical differences. the first
c     perturbation vector is based on the initial derivatives and
c     direction of integration. the second perturbation vector is
c     formed using another evaluation of the differential equation.
c     the third perturbation vector is formed using perturbations based
c     only on the initial values. components that are zero are always
c     changed to non-zero values (except on the first iteration). when
c     information is available, care is taken to ensure that components
c     of the perturbation vector have signs which are consistent with
c     the slopes of local solution curves.
c     also choose the largest bound (fbnd) for the first derivative.
c     no attempt is made to keep the perturbation vector size constant.
c
      if (ypnorm .eq. 0.) go to 60
c                       use initial derivatives for first perturbation
      icase = 1
      do 50 j=1,neq
        spy(j) = yprime(j)
   50   yp(j) = yprime(j)
      go to 80
c                       cannot have a null perturbation vector
   60 icase = 2
      do 70 j=1,neq
        spy(j) = yprime(j)
   70   yp(j) = etol(j)
c
   80 dfdub = 0.
      lk = min(neq+1,3)
      do 260 k=1,lk
c                       set ypnorm and delx
        ypnorm = hvnrm(yp,neq)
        if (icase .eq. 1  .or.  icase .eq. 3) go to 90
        delx = sign(1.0,dx)
        go to 120
c                       try to enforce meaningful perturbation values
   90   delx = dx
        if (abs(delx)*ypnorm .ge. relper*ynorm) go to 100
        delxb = big
        if (relper*ynorm .lt. big*ypnorm) delxb = relper*ynorm/ypnorm
        delx = sign(delxb,dx)
  100   do 110 j=1,neq
          if (abs(delx*yp(j)) .gt. etol(j)) delx=sign(etol(j)/yp(j),dx)
  110     continue
c                       define perturbed vector of initial values
  120   do 130 j=1,neq
  130     pv(j) = y(j) + delx*yp(j)
        if (k .eq. 2) go to 150
c                       evaluate derivatives associated with perturbed
c                       vector  and  compute corresponding differences
        call f(a,pv,yp,rpar,ipar)
        do 140 j=1,neq
  140     pv(j) = yp(j) - yprime(j)
        go to 170
c                       use a shifted value of the independent variable
c                                             in computing one estimate
  150   call f(a+da,pv,yp,rpar,ipar)
        do 160 j=1,neq
  160     pv(j) = yp(j) - sf(j)
c                       choose largest bound on the weighted first
c                                                   derivative
  170   if (morder .eq. 1) go to 190
        do 180 j=1,neq
  180     yp(j) = yp(j)/etol(j)**power
        go to 210
  190   do 200 j=1,neq
  200     yp(j) = yp(j)/etol(j)
  210   fbnd = max(fbnd,hvnrm(yp,neq))
c                       compute bound on a local lipschitz constant
        delf = hvnrm(pv,neq)
        if (delf .eq. 0.) go to 220
        dely = abs(delx)*ypnorm
        if (delf .ge. big*dely) go to 270
        dfdub = max(dfdub,delf/dely)
c
  220   if (k .eq. lk) go to 280
c                       choose next perturbation vector
        do 250 j=1,neq
          if (k .eq. lk-1) go to 230
          icase = 3
          dy = abs(pv(j))
          if (dy .eq. 0.) dy = max(delf,etol(j))
          go to 240
  230     icase = 4
          dy = max(relper*abs(y(j)),etol(j))
  240     if (spy(j) .eq. 0.) spy(j) = yp(j)
          if (spy(j) .ne. 0.) dy = sign(dy,spy(j))
  250     yp(j) = dy
  260   continue
c
c                       protect against an overflow
  270 dfdub = big
c
c.......................................................................
c
c     compute a bound (ydpb) on the norm of the second derivative
c
  280 ydpb = dfdxb + dfdub*fbnd
c
c.......................................................................
c
c     compute a starting step size based on the above first and second
c     derivative information
c
c                       restrict the step length to be not bigger than
c                       abs(b-a).   (unless  b  is too close to  a)
      h = absdx
c
      if (ydpb .ne. 0.  .or.  fbnd .ne. 0.) go to 290
c
c                       both first derivative term (fbnd) and second
c                                    derivative term (ydpb) are zero
      go to 310
c
  290 if (ydpb .ne. 0.) go to 300
c
c                       only second derivative term (ydpb) is zero
      if (1.0 .lt. fbnd*absdx) h = 1./fbnd
      go to 310
c
c                       second derivative term (ydpb) is non-zero
  300 srydpb = sqrt(0.5*ydpb)
      if (1.0 .lt. srydpb*absdx) h = 1./srydpb
c
c                       further restrict the step length to be not
c                                                 bigger than  1/dfdub
  310 if (h*dfdub .gt. 1.) h = 1./dfdub
c
c                       finally, restrict the step length to be not
c                       smaller than  100*small*abs(a).  however, if
c                       a=0. and the computed h underflowed to zero,
c                       the algorithm returns  small*abs(b)  for the
c                                                       step length.
      h = max(h,100.*small*abs(a))
      if (h .eq. 0.) h = small*abs(b)
c
c                       now set direction of integration
      h = sign(h,dx)
c
      return
      end
