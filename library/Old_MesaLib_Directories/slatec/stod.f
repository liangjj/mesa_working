*deck stod
      subroutine stod (neq, y, yh, nyh, yh1, ewt, savf, acor, wm, iwm,
     +   f, jac, rpar, ipar)
c***begin prologue  stod
c***subsidiary
c***purpose  subsidiary to debdf
c***library   slatec
c***type      single precision (stod-s, dstod-d)
c***author  watts, h. a., (snla)
c***description
c
c   stod integrates a system of first order odes over one step in the
c   integrator package debdf.
c ----------------------------------------------------------------------
c stod  performs one step of the integration of an initial value
c problem for a system of ordinary differential equations.
c note.. stod  is independent of the value of the iteration method
c indicator miter, when this is .ne. 0, and hence is independent
c of the type of chord method used, or the jacobian structure.
c communication with stod  is done with the following variables..
c
c y      = an array of length .ge. n used as the y argument in
c          all calls to f and jac.
c neq    = integer array containing problem size in neq(1), and
c          passed as the neq argument in all calls to f and jac.
c yh     = an nyh by lmax array containing the dependent variables
c          and their approximate scaled derivatives, where
c          lmax = maxord + 1.  yh(i,j+1) contains the approximate
c          j-th derivative of y(i), scaled by h**j/factorial(j)
c          (j = 0,1,...,nq).  on entry for the first step, the first
c          two columns of yh must be set from the initial values.
c nyh    = a constant integer .ge. n, the first dimension of yh.
c yh1    = a one-dimensional array occupying the same space as yh.
c ewt    = an array of n elements with which the estimated local
c          errors in yh are compared.
c savf   = an array of working storage, of length n.
c acor   = a work array of length n, used for the accumulated
c          corrections.  on a successful return, acor(i) contains
c          the estimated one-step local error in y(i).
c wm,iwm = real and integer work arrays associated with matrix
c          operations in chord iteration (miter .ne. 0).
c pjac   = name of routine to evaluate and preprocess jacobian matrix
c          if a chord method is being used.
c slvs   = name of routine to solve linear system in chord iteration.
c h      = the step size to be attempted on the next step.
c          h is altered by the error control algorithm during the
c          problem.  h can be either positive or negative, but its
c          sign must remain constant throughout the problem.
c hmin   = the minimum absolute value of the step size h to be used.
c hmxi   = inverse of the maximum absolute value of h to be used.
c          hmxi = 0.0 is allowed and corresponds to an infinite hmax.
c          hmin and hmxi may be changed at any time, but will not
c          take effect until the next change of h is considered.
c tn     = the independent variable. tn is updated on each step taken.
c jstart = an integer used for input only, with the following
c          values and meanings..
c               0  perform the first step.
c           .gt.0  take a new step continuing from the last.
c              -1  take the next step with a new value of h, maxord,
c                    n, meth, miter, and/or matrix parameters.
c              -2  take the next step with a new value of h,
c                    but with other inputs unchanged.
c          on return, jstart is set to 1 to facilitate continuation.
c kflag  = a completion code with the following meanings..
c               0  the step was successful.
c              -1  the requested error could not be achieved.
c              -2  corrector convergence could not be achieved.
c          a return with kflag = -1 or -2 means either
c          abs(h) = hmin or 10 consecutive failures occurred.
c          on a return with kflag negative, the values of tn and
c          the yh array are as of the beginning of the last
c          step, and h is the last step size attempted.
c maxord = the maximum order of integration method to be allowed.
c meth/miter = the method flags.  see description in driver.
c n      = the number of first-order differential equations.
c ----------------------------------------------------------------------
c
c***see also  debdf
c***routines called  cfod, pjac, slvs, vnwrms
c***common blocks    debdf1
c***revision history  (yymmdd)
c   800901  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910722  updated author section.  (als)
c   920422  changed dimension statement.  (wrb)
c***end prologue  stod
      external f, jac
c
clll. optimize
      integer neq, nyh, iwm, i, i1, ialth, ier, iownd, iredo, iret,
     1   ipup, j, jb, jstart, kflag, l, lmax, m, maxord, meo, meth,
     2   miter, n, ncf, newq, nfe, nje, nq, nqnyh, nqu, nst, nstepj
      real y, yh, yh1, ewt, savf, acor, wm,
     1   rownd, conit, crate, el, elco, hold, rc, rmax, tesco,
     2   el0, h, hmin, hmxi, hu, tn, uround,
     3   dcon, ddn, del, delp, dsm, dup, exdn, exsm, exup,
     4   r, rh, rhdn, rhsm, rhup, told, vnwrms
      dimension         y(*), yh(nyh,*), yh1(*), ewt(*), savf(*),
     1   acor(*), wm(*), iwm(*), rpar(*), ipar(*)
      common /debdf1/ rownd, conit, crate, el(13), elco(13,12),
     1   hold, rc, rmax, tesco(3,12),
     2   el0, h, hmin, hmxi, hu, tn, uround, iownd(7), ksteps, iod(6),
     3   ialth, ipup, lmax, meo, nqnyh, nstepj,
     4   ier, jstart, kflag, l, meth, miter, maxord, n, nq, nst, nfe,
     5   nje, nqu
c
c
c***first executable statement  stod
      kflag = 0
      told = tn
      ncf = 0
      if (jstart .gt. 0) go to 200
      if (jstart .eq. -1) go to 100
      if (jstart .eq. -2) go to 160
c-----------------------------------------------------------------------
c on the first call, the order is set to 1, and other variables are
c initialized.  rmax is the maximum ratio by which h can be increased
c in a single step.  it is initially 1.e4 to compensate for the small
c initial h, but then is normally equal to 10.  if a failure
c occurs (in corrector convergence or error test), rmax is set at 2
c for the next increase.
c-----------------------------------------------------------------------
      lmax = maxord + 1
      nq = 1
      l = 2
      ialth = 2
      rmax = 10000.0e0
      rc = 0.0e0
      el0 = 1.0e0
      crate = 0.7e0
      delp = 0.0e0
      hold = h
      meo = meth
      nstepj = 0
      iret = 3
      go to 140
c-----------------------------------------------------------------------
c the following block handles preliminaries needed when jstart = -1.
c ipup is set to miter to force a matrix update.
c if an order increase is about to be considered (ialth = 1),
c ialth is reset to 2 to postpone consideration one more step.
c if the caller has changed meth, cfod  is called to reset
c the coefficients of the method.
c if the caller has changed maxord to a value less than the current
c order nq, nq is reduced to maxord, and a new h chosen accordingly.
c if h is to be changed, yh must be rescaled.
c if h or meth is being changed, ialth is reset to l = nq + 1
c to prevent further changes in h for that many steps.
c-----------------------------------------------------------------------
 100  ipup = miter
      lmax = maxord + 1
      if (ialth .eq. 1) ialth = 2
      if (meth .eq. meo) go to 110
      call cfod  (meth, elco, tesco)
      meo = meth
      if (nq .gt. maxord) go to 120
      ialth = l
      iret = 1
      go to 150
 110  if (nq .le. maxord) go to 160
 120  nq = maxord
      l = lmax
      do 125 i = 1,l
 125    el(i) = elco(i,nq)
      nqnyh = nq*nyh
      rc = rc*el(1)/el0
      el0 = el(1)
      conit = 0.5e0/(nq+2)
      ddn = vnwrms (n, savf, ewt)/tesco(1,l)
      exdn = 1.0e0/l
      rhdn = 1.0e0/(1.3e0*ddn**exdn + 0.0000013e0)
      rh = min(rhdn,1.0e0)
      iredo = 3
      if (h .eq. hold) go to 170
      rh = min(rh,abs(h/hold))
      h = hold
      go to 175
c-----------------------------------------------------------------------
c cfod  is called to get all the integration coefficients for the
c current meth.  then the el vector and related constants are reset
c whenever the order nq is changed, or at the start of the problem.
c-----------------------------------------------------------------------
 140  call cfod  (meth, elco, tesco)
 150  do 155 i = 1,l
 155    el(i) = elco(i,nq)
      nqnyh = nq*nyh
      rc = rc*el(1)/el0
      el0 = el(1)
      conit = 0.5e0/(nq+2)
      go to (160, 170, 200), iret
c-----------------------------------------------------------------------
c if h is being changed, the h ratio rh is checked against
c rmax, hmin, and hmxi, and the yh array rescaled.  ialth is set to
c l = nq + 1 to prevent a change of h for that many steps, unless
c forced by a convergence or error test failure.
c-----------------------------------------------------------------------
 160  if (h .eq. hold) go to 200
      rh = h/hold
      h = hold
      iredo = 3
      go to 175
 170  rh = max(rh,hmin/abs(h))
 175  rh = min(rh,rmax)
      rh = rh/max(1.0e0,abs(h)*hmxi*rh)
      r = 1.0e0
      do 180 j = 2,l
        r = r*rh
        do 180 i = 1,n
 180      yh(i,j) = yh(i,j)*r
      h = h*rh
      rc = rc*rh
      ialth = l
      if (iredo .eq. 0) go to 680
c-----------------------------------------------------------------------
c this section computes the predicted values by effectively
c multiplying the yh array by the pascal triangle matrix.
c rc is the ratio of new to old values of the coefficient  h*el(1).
c when rc differs from 1 by more than 30 percent, ipup is set to miter
c to force pjac to be called, if a jacobian is involved.
c in any case, pjac is called at least every 20-th step.
c-----------------------------------------------------------------------
 200  if (abs(rc-1.0e0) .gt. 0.3e0) ipup = miter
      if (nst .ge. nstepj+20) ipup = miter
      tn = tn + h
      i1 = nqnyh + 1
      do 215 jb = 1,nq
        i1 = i1 - nyh
        do 210 i = i1,nqnyh
 210      yh1(i) = yh1(i) + yh1(i+nyh)
 215    continue
      ksteps = ksteps + 1
c-----------------------------------------------------------------------
c up to 3 corrector iterations are taken.  a convergence test is
c made on the r.m.s. norm of each correction, weighted by the error
c weight vector ewt.  the sum of the corrections is accumulated in the
c vector acor(i).  the yh array is not altered in the corrector loop.
c-----------------------------------------------------------------------
 220  m = 0
      do 230 i = 1,n
 230    y(i) = yh(i,1)
      call f (tn, y, savf, rpar, ipar)
      nfe = nfe + 1
      if (ipup .le. 0) go to 250
c-----------------------------------------------------------------------
c if indicated, the matrix p = i - h*el(1)*j is reevaluated and
c preprocessed before starting the corrector iteration.  ipup is set
c to 0 as an indicator that this has been done.
c-----------------------------------------------------------------------
      ipup = 0
      rc = 1.0e0
      nstepj = nst
      crate = 0.7e0
      call pjac (neq, y, yh, nyh, ewt, acor, savf, wm, iwm, f, jac,
     1          rpar, ipar)
      if (ier .ne. 0) go to 430
 250  do 260 i = 1,n
 260    acor(i) = 0.0e0
 270  if (miter .ne. 0) go to 350
c-----------------------------------------------------------------------
c in the case of functional iteration, update y directly from
c the result of the last function evaluation.
c-----------------------------------------------------------------------
      do 290 i = 1,n
        savf(i) = h*savf(i) - yh(i,2)
 290    y(i) = savf(i) - acor(i)
      del = vnwrms (n, y, ewt)
      do 300 i = 1,n
        y(i) = yh(i,1) + el(1)*savf(i)
 300    acor(i) = savf(i)
      go to 400
c-----------------------------------------------------------------------
c in the case of the chord method, compute the corrector error,
c and solve the linear system with that as right-hand side and
c p as coefficient matrix.
c-----------------------------------------------------------------------
 350  do 360 i = 1,n
 360    y(i) = h*savf(i) - (yh(i,2) + acor(i))
      call slvs (wm, iwm, y, savf)
      if (ier .ne. 0) go to 410
      del = vnwrms (n, y, ewt)
      do 380 i = 1,n
        acor(i) = acor(i) + y(i)
 380    y(i) = yh(i,1) + el(1)*acor(i)
c-----------------------------------------------------------------------
c test for convergence.  if m.gt.0, an estimate of the convergence
c rate constant is stored in crate, and this is used in the test.
c-----------------------------------------------------------------------
 400  if (m .ne. 0) crate = max(0.2e0*crate,del/delp)
      dcon = del*min(1.0e0,1.5e0*crate)/(tesco(2,nq)*conit)
      if (dcon .le. 1.0e0) go to 450
      m = m + 1
      if (m .eq. 3) go to 410
      if (m .ge. 2 .and. del .gt. 2.0e0*delp) go to 410
      delp = del
      call f (tn, y, savf, rpar, ipar)
      nfe = nfe + 1
      go to 270
c-----------------------------------------------------------------------
c the corrector iteration failed to converge in 3 tries.
c if miter .ne. 0 and the jacobian is out of date, pjac is called for
c the next try.  otherwise the yh array is retracted to its values
c before prediction, and h is reduced, if possible.  if h cannot be
c reduced or 10 failures have occurred, exit with kflag = -2.
c-----------------------------------------------------------------------
 410  if (ipup .eq. 0) go to 430
      ipup = miter
      go to 220
 430  tn = told
      ncf = ncf + 1
      rmax = 2.0e0
      i1 = nqnyh + 1
      do 445 jb = 1,nq
        i1 = i1 - nyh
        do 440 i = i1,nqnyh
 440      yh1(i) = yh1(i) - yh1(i+nyh)
 445    continue
      if (abs(h) .le. hmin*1.00001e0) go to 670
      if (ncf .eq. 10) go to 670
      rh = 0.25e0
      ipup = miter
      iredo = 1
      go to 170
c-----------------------------------------------------------------------
c the corrector has converged.  ipup is set to -1 if miter .ne. 0,
c to signal that the jacobian involved may need updating later.
c the local error test is made and control passes to statement 500
c if it fails.
c-----------------------------------------------------------------------
 450  if (miter .ne. 0) ipup = -1
      if (m .eq. 0) dsm = del/tesco(2,nq)
      if (m .gt. 0) dsm = vnwrms (n, acor, ewt)/tesco(2,nq)
      if (dsm .gt. 1.0e0) go to 500
c-----------------------------------------------------------------------
c after a successful step, update the yh array.
c consider changing h if ialth = 1.  otherwise decrease ialth by 1.
c if ialth is then 1 and nq .lt. maxord, then acor is saved for
c use in a possible order increase on the next step.
c if a change in h is considered, an increase or decrease in order
c by one is considered also.  a change in h is made only if it is by a
c factor of at least 1.1.  if not, ialth is set to 3 to prevent
c testing for that many steps.
c-----------------------------------------------------------------------
      kflag = 0
      iredo = 0
      nst = nst + 1
      hu = h
      nqu = nq
      do 470 j = 1,l
        do 470 i = 1,n
 470      yh(i,j) = yh(i,j) + el(j)*acor(i)
      ialth = ialth - 1
      if (ialth .eq. 0) go to 520
      if (ialth .gt. 1) go to 690
      if (l .eq. lmax) go to 690
      do 490 i = 1,n
 490    yh(i,lmax) = acor(i)
      go to 690
c-----------------------------------------------------------------------
c the error test failed.  kflag keeps track of multiple failures.
c restore tn and the yh array to their previous values, and prepare
c to try the step again.  compute the optimum step size for this or
c one lower order.  after 2 or more failures, h is forced to decrease
c by a factor of 0.2 or less.
c-----------------------------------------------------------------------
 500  kflag = kflag - 1
      tn = told
      i1 = nqnyh + 1
      do 515 jb = 1,nq
        i1 = i1 - nyh
        do 510 i = i1,nqnyh
 510      yh1(i) = yh1(i) - yh1(i+nyh)
 515    continue
      rmax = 2.0e0
      if (abs(h) .le. hmin*1.00001e0) go to 660
      if (kflag .le. -3) go to 640
      iredo = 2
      rhup = 0.0e0
      go to 540
c-----------------------------------------------------------------------
c regardless of the success or failure of the step, factors
c rhdn, rhsm, and rhup are computed, by which h could be multiplied
c at order nq - 1, order nq, or order nq + 1, respectively.
c in the case of failure, rhup = 0.0 to avoid an order increase.
c the largest of these is determined and the new order chosen
c accordingly.  if the order is to be increased, we compute one
c additional scaled derivative.
c-----------------------------------------------------------------------
 520  rhup = 0.0e0
      if (l .eq. lmax) go to 540
      do 530 i = 1,n
 530    savf(i) = acor(i) - yh(i,lmax)
      dup = vnwrms (n, savf, ewt)/tesco(3,nq)
      exup = 1.0e0/(l+1)
      rhup = 1.0e0/(1.4e0*dup**exup + 0.0000014e0)
 540  exsm = 1.0e0/l
      rhsm = 1.0e0/(1.2e0*dsm**exsm + 0.0000012e0)
      rhdn = 0.0e0
      if (nq .eq. 1) go to 560
      ddn = vnwrms (n, yh(1,l), ewt)/tesco(1,nq)
      exdn = 1.0e0/nq
      rhdn = 1.0e0/(1.3e0*ddn**exdn + 0.0000013e0)
 560  if (rhsm .ge. rhup) go to 570
      if (rhup .gt. rhdn) go to 590
      go to 580
 570  if (rhsm .lt. rhdn) go to 580
      newq = nq
      rh = rhsm
      go to 620
 580  newq = nq - 1
      rh = rhdn
      if (kflag .lt. 0 .and. rh .gt. 1.0e0) rh = 1.0e0
      go to 620
 590  newq = l
      rh = rhup
      if (rh .lt. 1.1e0) go to 610
      r = el(l)/l
      do 600 i = 1,n
 600    yh(i,newq+1) = acor(i)*r
      go to 630
 610  ialth = 3
      go to 690
 620  if ((kflag .eq. 0) .and. (rh .lt. 1.1e0)) go to 610
      if (kflag .le. -2) rh = min(rh,0.2e0)
c-----------------------------------------------------------------------
c if there is a change of order, reset nq, l, and the coefficients.
c in any case h is reset according to rh and the yh array is rescaled.
c then exit from 680 if the step was ok, or redo the step otherwise.
c-----------------------------------------------------------------------
      if (newq .eq. nq) go to 170
 630  nq = newq
      l = nq + 1
      iret = 2
      go to 150
c-----------------------------------------------------------------------
c control reaches this section if 3 or more failures have occurred.
c if 10 failures have occurred, exit with kflag = -1.
c it is assumed that the derivatives that have accumulated in the
c yh array have errors of the wrong order.  hence the first
c derivative is recomputed, and the order is set to 1.  then
c h is reduced by a factor of 10, and the step is retried,
c until it succeeds or h reaches hmin.
c-----------------------------------------------------------------------
 640  if (kflag .eq. -10) go to 660
      rh = 0.1e0
      rh = max(hmin/abs(h),rh)
      h = h*rh
      do 645 i = 1,n
 645    y(i) = yh(i,1)
      call f (tn, y, savf, rpar, ipar)
      nfe = nfe + 1
      do 650 i = 1,n
 650    yh(i,2) = h*savf(i)
      ipup = miter
      ialth = 5
      if (nq .eq. 1) go to 200
      nq = 1
      l = 2
      iret = 3
      go to 150
c-----------------------------------------------------------------------
c all returns are made through this section.  h is saved in hold
c to allow the caller to change h on the next step.
c-----------------------------------------------------------------------
 660  kflag = -1
      go to 700
 670  kflag = -2
      go to 700
 680  rmax = 10.0e0
 690  r = 1.0e0/tesco(2,nqu)
      do 695 i = 1,n
 695    acor(i) = acor(i)*r
 700  hold = h
      jstart = 1
      return
c----------------------- end of subroutine stod  -----------------------
      end
