*deck dhstrt
      subroutine dhstrt (df, neq, a, b, y, yprime, etol, morder, small,
     +   big, spy, pv, yp, sf, rpar, ipar, h)
c***begin prologue  dhstrt
c***subsidiary
c***purpose  subsidiary to ddeabm, ddebdf and dderkf
c***library   slatec
c***type      double precision (hstart-s, dhstrt-d)
c***author  watts, h. a., (snla)
c***description
c
c   dhstrt computes a starting step size to be used in solving initial
c   value problems in ordinary differential equations.
c
c **********************************************************************
c  abstract
c
c     subroutine dhstrt computes a starting step size to be used by an
c     initial value method in solving ordinary differential equations.
c     it is based on an estimate of the local lipschitz constant for the
c     differential equation   (lower bound on a norm of the jacobian) ,
c     a bound on the differential equation  (first derivative) , and
c     a bound on the partial derivative of the equation with respect to
c     the independent variable.
c     (all approximated near the initial point a)
c
c     subroutine dhstrt uses a function subprogram dhvnrm for computing
c     a vector norm. the maximum norm is presently utilized though it
c     can easily be replaced by any other vector norm. it is presumed
c     that any replacement norm routine would be carefully coded to
c     prevent unnecessary underflows or overflows from occurring, and
c     also, would not alter the vector or number of components.
c
c **********************************************************************
c  on input you must provide the following
c
c      df -- this is a subroutine of the form
c                               df(x,u,uprime,rpar,ipar)
c             which defines the system of first order differential
c             equations to be solved. for the given values of x and the
c             vector  u(*)=(u(1),u(2),...,u(neq)) , the subroutine must
c             evaluate the neq components of the system of differential
c             equations  du/dx=df(x,u)  and store the derivatives in the
c             array uprime(*), that is,  uprime(i) = * du(i)/dx *  for
c             equations i=1,...,neq.
c
c             subroutine df must not alter x or u(*). you must declare
c             the name df in an external statement in your program that
c             calls dhstrt. you must dimension u and uprime in df.
c
c             rpar and ipar are double precision and integer parameter
c             arrays which you can use for communication between your
c             program and subroutine df. they are not used or altered by
c             dhstrt. if you do not need rpar or ipar, ignore these
c             parameters by treating them as dummy arguments. if you do
c             choose to use them, dimension them in your program and in
c             df as arrays of appropriate length.
c
c      neq -- this is the number of (first order) differential equations
c             to be integrated.
c
c      a -- this is the initial point of integration.
c
c      b -- this is a value of the independent variable used to define
c             the direction of integration. a reasonable choice is to
c             set  b  to the first point at which a solution is desired.
c             you can also use  b, if necessary, to restrict the length
c             of the first integration step because the algorithm will
c             not compute a starting step length which is bigger than
c             abs(b-a), unless  b  has been chosen too close to  a.
c             (it is presumed that dhstrt has been called with  b
c             different from  a  on the machine being used. also see the
c             discussion about the parameter  small.)
c
c      y(*) -- this is the vector of initial values of the neq solution
c             components at the initial point  a.
c
c      yprime(*) -- this is the vector of derivatives of the neq
c             solution components at the initial point  a.
c             (defined by the differential equations in subroutine df)
c
c      etol -- this is the vector of error tolerances corresponding to
c             the neq solution components. it is assumed that all
c             elements are positive. following the first integration
c             step, the tolerances are expected to be used by the
c             integrator in an error test which roughly requires that
c                        abs(local error)  .le.  etol
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
c             (approximately) the smallest positive double precision
c             number such that  (1.+small) .gt. 1.  on the machine being
c             used. the quantity  small**(3/8)  is used in computing
c             increments of variables for approximating derivatives by
c             differences.  also the algorithm will not compute a
c             starting step length which is smaller than
c             100*small*abs(a).
c
c      big -- this is a large positive machine dependent constant which
c             is used for preventing machine overflows. a reasonable
c             choice is to set big to (approximately) the square root of
c             the largest double precision number which can be held in
c             the machine.
c
c      spy(*),pv(*),yp(*),sf(*) -- these are double precision work
c             arrays of length neq which provide the routine with needed
c             storage space.
c
c      rpar,ipar -- these are parameter arrays, of double precision and
c             integer type, respectively, which can be used for
c             communication between your program and the df subroutine.
c             they are not used or altered by dhstrt.
c
c **********************************************************************
c  on output  (after the return from dhstrt),
c
c      h -- is an appropriate starting step size to be attempted by the
c             differential equation method.
c
c           all parameters in the call list remain unchanged except for
c           the working arrays spy(*),pv(*),yp(*), and sf(*).
c
c **********************************************************************
c
c***see also  ddeabm, ddebdf, dderkf
c***routines called  dhvnrm
c***revision history  (yymmdd)
c   820301  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   890911  removed unnecessary intrinsics.  (wrb)
c   891024  changed references from dvnorm to dhvnrm.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910722  updated author section.  (als)
c***end prologue  dhstrt
c
      integer ipar, j, k, lk, morder, neq
      double precision a, absdx, b, big, da, delf, dely,
     1      dfdub, dfdxb, dhvnrm,
     2      dx, dy, etol, fbnd, h, pv, relper, rpar, sf, small, spy,
     3      srydpb, tolexp, tolmin, tolp, tolsum, y, ydpb, yp, yprime
      dimension y(*),yprime(*),etol(*),spy(*),pv(*),yp(*),
     1          sf(*),rpar(*),ipar(*)
      external df
c
c     ..................................................................
c
c     begin block permitting ...exits to 160
c***first executable statement  dhstrt
         dx = b - a
         absdx = abs(dx)
         relper = small**0.375d0
c
c        ...............................................................
c
c             compute an approximate bound (dfdxb) on the partial
c             derivative of the equation with respect to the
c             independent variable. protect against an overflow.
c             also compute a bound (fbnd) on the first derivative
c             locally.
c
         da = sign(max(min(relper*abs(a),absdx),
     1                    100.0d0*small*abs(a)),dx)
         if (da .eq. 0.0d0) da = relper*dx
         call df(a+da,y,sf,rpar,ipar)
         do 10 j = 1, neq
            yp(j) = sf(j) - yprime(j)
   10    continue
         delf = dhvnrm(yp,neq)
         dfdxb = big
         if (delf .lt. big*abs(da)) dfdxb = delf/abs(da)
         fbnd = dhvnrm(sf,neq)
c
c        ...............................................................
c
c             compute an estimate (dfdub) of the local lipschitz
c             constant for the system of differential equations. this
c             also represents an estimate of the norm of the jacobian
c             locally.  three iterations (two when neq=1) are used to
c             estimate the lipschitz constant by numerical differences.
c             the first perturbation vector is based on the initial
c             derivatives and direction of integration. the second
c             perturbation vector is formed using another evaluation of
c             the differential equation.  the third perturbation vector
c             is formed using perturbations based only on the initial
c             values. components that are zero are always changed to
c             non-zero values (except on the first iteration). when
c             information is available, care is taken to ensure that
c             components of the perturbation vector have signs which are
c             consistent with the slopes of local solution curves.
c             also choose the largest bound (fbnd) for the first
c             derivative.
c
c                               perturbation vector size is held
c                               constant for all iterations. compute
c                               this change from the
c                                       size of the vector of initial
c                                       values.
         dely = relper*dhvnrm(y,neq)
         if (dely .eq. 0.0d0) dely = relper
         dely = sign(dely,dx)
         delf = dhvnrm(yprime,neq)
         fbnd = max(fbnd,delf)
         if (delf .eq. 0.0d0) go to 30
c           use initial derivatives for first perturbation
            do 20 j = 1, neq
               spy(j) = yprime(j)
               yp(j) = yprime(j)
   20       continue
         go to 50
   30    continue
c           cannot have a null perturbation vector
            do 40 j = 1, neq
               spy(j) = 0.0d0
               yp(j) = 1.0d0
   40       continue
            delf = dhvnrm(yp,neq)
   50    continue
c
         dfdub = 0.0d0
         lk = min(neq+1,3)
         do 140 k = 1, lk
c           define perturbed vector of initial values
            do 60 j = 1, neq
               pv(j) = y(j) + dely*(yp(j)/delf)
   60       continue
            if (k .eq. 2) go to 80
c              evaluate derivatives associated with perturbed
c              vector  and  compute corresponding differences
               call df(a,pv,yp,rpar,ipar)
               do 70 j = 1, neq
                  pv(j) = yp(j) - yprime(j)
   70          continue
            go to 100
   80       continue
c              use a shifted value of the independent variable
c                                    in computing one estimate
               call df(a+da,pv,yp,rpar,ipar)
               do 90 j = 1, neq
                  pv(j) = yp(j) - sf(j)
   90          continue
  100       continue
c           choose largest bounds on the first derivative
c                          and a local lipschitz constant
            fbnd = max(fbnd,dhvnrm(yp,neq))
            delf = dhvnrm(pv,neq)
c        ...exit
            if (delf .ge. big*abs(dely)) go to 150
            dfdub = max(dfdub,delf/abs(dely))
c     ......exit
            if (k .eq. lk) go to 160
c           choose next perturbation vector
            if (delf .eq. 0.0d0) delf = 1.0d0
            do 130 j = 1, neq
               if (k .eq. 2) go to 110
                  dy = abs(pv(j))
                  if (dy .eq. 0.0d0) dy = delf
               go to 120
  110          continue
                  dy = y(j)
                  if (dy .eq. 0.0d0) dy = dely/relper
  120          continue
               if (spy(j) .eq. 0.0d0) spy(j) = yp(j)
               if (spy(j) .ne. 0.0d0) dy = sign(dy,spy(j))
               yp(j) = dy
  130       continue
            delf = dhvnrm(yp,neq)
  140    continue
  150    continue
c
c        protect against an overflow
         dfdub = big
  160 continue
c
c     ..................................................................
c
c          compute a bound (ydpb) on the norm of the second derivative
c
      ydpb = dfdxb + dfdub*fbnd
c
c     ..................................................................
c
c          define the tolerance parameter upon which the starting step
c          size is to be based.  a value in the middle of the error
c          tolerance range is selected.
c
      tolmin = big
      tolsum = 0.0d0
      do 170 k = 1, neq
         tolexp = log10(etol(k))
         tolmin = min(tolmin,tolexp)
         tolsum = tolsum + tolexp
  170 continue
      tolp = 10.0d0**(0.5d0*(tolsum/neq + tolmin)/(morder+1))
c
c     ..................................................................
c
c          compute a starting step size based on the above first and
c          second derivative information
c
c                            restrict the step length to be not bigger
c                            than abs(b-a).   (unless  b  is too close
c                            to  a)
      h = absdx
c
      if (ydpb .ne. 0.0d0 .or. fbnd .ne. 0.0d0) go to 180
c
c        both first derivative term (fbnd) and second
c                     derivative term (ydpb) are zero
         if (tolp .lt. 1.0d0) h = absdx*tolp
      go to 200
  180 continue
c
      if (ydpb .ne. 0.0d0) go to 190
c
c        only second derivative term (ydpb) is zero
         if (tolp .lt. fbnd*absdx) h = tolp/fbnd
      go to 200
  190 continue
c
c        second derivative term (ydpb) is non-zero
         srydpb = sqrt(0.5d0*ydpb)
         if (tolp .lt. srydpb*absdx) h = tolp/srydpb
  200 continue
c
c     further restrict the step length to be not
c                               bigger than  1/dfdub
      if (h*dfdub .gt. 1.0d0) h = 1.0d0/dfdub
c
c     finally, restrict the step length to be not
c     smaller than  100*small*abs(a).  however, if
c     a=0. and the computed h underflowed to zero,
c     the algorithm returns  small*abs(b)  for the
c                                     step length.
      h = max(h,100.0d0*small*abs(a))
      if (h .eq. 0.0d0) h = small*abs(b)
c
c     now set direction of integration
      h = sign(h,dx)
c
      return
      end
