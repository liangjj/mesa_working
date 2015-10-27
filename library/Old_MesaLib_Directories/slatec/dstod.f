*deck dstod
      subroutine dstod (neq, y, yh, nyh, yh1, ewt, savf, acor, wm, iwm,
     +   df, djac, rpar, ipar)
c***begin prologue  dstod
c***subsidiary
c***purpose  subsidiary to ddebdf
c***library   slatec
c***type      double precision (stod-s, dstod-d)
c***author  watts, h. a., (snla)
c***description
c
c   dstod integrates a system of first order odes over one step in the
c   integrator package ddebdf.
c ----------------------------------------------------------------------
c dstod  performs one step of the integration of an initial value
c problem for a system of ordinary differential equations.
c note.. dstod  is independent of the value of the iteration method
c indicator miter, when this is .ne. 0, and hence is independent
c of the type of chord method used, or the jacobian structure.
c communication with dstod  is done with the following variables..
c
c y      = an array of length .ge. n used as the y argument in
c          all calls to df and djac.
c neq    = integer array containing problem size in neq(1), and
c          passed as the neq argument in all calls to df and djac.
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
c wm,iwm = double precision and integer work arrays associated with
c          matrix operations in chord iteration (miter .ne. 0).
c dpjac   = name of routine to evaluate and preprocess jacobian matrix
c          if a chord method is being used.
c dslvs   = name of routine to solve linear system in chord iteration.
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
c***see also  ddebdf
c***routines called  dcfod, dpjac, dslvs, dvnrms
c***common blocks    ddebd1
c***revision history  (yymmdd)
c   820301  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890911  removed unnecessary intrinsics.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910722  updated author section.  (als)
c   920422  changed dimension statement.  (wrb)
c***end prologue  dstod
c
      integer i, i1, ialth, ier, iod, iownd, ipar, ipup, iredo, iret,
     1      iwm, j, jb, jstart, kflag, ksteps, l, lmax, m, maxord,
     2      meo, meth, miter, n, ncf, neq, newq, nfe, nje, nq, nqnyh,
     3      nqu, nst, nstepj, nyh
      double precision acor, conit, crate, dcon, ddn,
     1      del, delp, dsm, dup, dvnrms, el, el0, elco,
     2      ewt, exdn, exsm, exup, h, hmin, hmxi, hold, hu, r, rc,
     3      rh, rhdn, rhsm, rhup, rmax, rownd, rpar, savf, tesco,
     4      tn, told, uround, wm, y, yh, yh1
      external df, djac
c
      dimension y(*),yh(nyh,*),yh1(*),ewt(*),savf(*),acor(*),wm(*),
     1          iwm(*),rpar(*),ipar(*)
      common /ddebd1/ rownd,conit,crate,el(13),elco(13,12),hold,rc,rmax,
     1                tesco(3,12),el0,h,hmin,hmxi,hu,tn,uround,iownd(7),
     2                ksteps,iod(6),ialth,ipup,lmax,meo,nqnyh,nstepj,
     3                ier,jstart,kflag,l,meth,miter,maxord,n,nq,nst,nfe,
     4                nje,nqu
c
c
c     begin block permitting ...exits to 690
c        begin block permitting ...exits to 60
c***first executable statement  dstod
            kflag = 0
            told = tn
            ncf = 0
            if (jstart .gt. 0) go to 160
            if (jstart .eq. -1) go to 10
               if (jstart .eq. -2) go to 90
c              ---------------------------------------------------------
c               on the first call, the order is set to 1, and other
c               variables are initialized.  rmax is the maximum ratio by
c               which h can be increased in a single step.  it is
c               initially 1.e4 to compensate for the small initial h,
c               but then is normally equal to 10.  if a failure occurs
c               (in corrector convergence or error test), rmax is set at
c               2 for the next increase.
c              ---------------------------------------------------------
               lmax = maxord + 1
               nq = 1
               l = 2
               ialth = 2
               rmax = 10000.0d0
               rc = 0.0d0
               el0 = 1.0d0
               crate = 0.7d0
               delp = 0.0d0
               hold = h
               meo = meth
               nstepj = 0
               iret = 3
            go to 50
   10       continue
c              begin block permitting ...exits to 30
c                 ------------------------------------------------------
c                  the following block handles preliminaries needed when
c                  jstart = -1.  ipup is set to miter to force a matrix
c                  update.  if an order increase is about to be
c                  considered (ialth = 1), ialth is reset to 2 to
c                  postpone consideration one more step.  if the caller
c                  has changed meth, dcfod  is called to reset the
c                  coefficients of the method.  if the caller has
c                  changed maxord to a value less than the current
c                  order nq, nq is reduced to maxord, and a new h chosen
c                  accordingly.  if h is to be changed, yh must be
c                  rescaled.  if h or meth is being changed, ialth is
c                  reset to l = nq + 1 to prevent further changes in h
c                  for that many steps.
c                 ------------------------------------------------------
                  ipup = miter
                  lmax = maxord + 1
                  if (ialth .eq. 1) ialth = 2
                  if (meth .eq. meo) go to 20
                     call dcfod(meth,elco,tesco)
                     meo = meth
c              ......exit
                     if (nq .gt. maxord) go to 30
                     ialth = l
                     iret = 1
c        ............exit
                     go to 60
   20             continue
                  if (nq .le. maxord) go to 90
   30          continue
               nq = maxord
               l = lmax
               do 40 i = 1, l
                  el(i) = elco(i,nq)
   40          continue
               nqnyh = nq*nyh
               rc = rc*el(1)/el0
               el0 = el(1)
               conit = 0.5d0/(nq+2)
               ddn = dvnrms(n,savf,ewt)/tesco(1,l)
               exdn = 1.0d0/l
               rhdn = 1.0d0/(1.3d0*ddn**exdn + 0.0000013d0)
               rh = min(rhdn,1.0d0)
               iredo = 3
               if (h .eq. hold) go to 660
               rh = min(rh,abs(h/hold))
               h = hold
               go to 100
   50       continue
c           ------------------------------------------------------------
c            dcfod  is called to get all the integration coefficients
c            for the current meth.  then the el vector and related
c            constants are reset whenever the order nq is changed, or at
c            the start of the problem.
c           ------------------------------------------------------------
            call dcfod(meth,elco,tesco)
   60    continue
   70    continue
c           begin block permitting ...exits to 680
               do 80 i = 1, l
                  el(i) = elco(i,nq)
   80          continue
               nqnyh = nq*nyh
               rc = rc*el(1)/el0
               el0 = el(1)
               conit = 0.5d0/(nq+2)
               go to (90,660,160), iret
c              ---------------------------------------------------------
c               if h is being changed, the h ratio rh is checked against
c               rmax, hmin, and hmxi, and the yh array rescaled.  ialth
c               is set to l = nq + 1 to prevent a change of h for that
c               many steps, unless forced by a convergence or error test
c               failure.
c              ---------------------------------------------------------
   90          continue
               if (h .eq. hold) go to 160
               rh = h/hold
               h = hold
               iredo = 3
  100          continue
  110          continue
                  rh = min(rh,rmax)
                  rh = rh/max(1.0d0,abs(h)*hmxi*rh)
                  r = 1.0d0
                  do 130 j = 2, l
                     r = r*rh
                     do 120 i = 1, n
                        yh(i,j) = yh(i,j)*r
  120                continue
  130             continue
                  h = h*rh
                  rc = rc*rh
                  ialth = l
                  if (iredo .ne. 0) go to 150
                     rmax = 10.0d0
                     r = 1.0d0/tesco(2,nqu)
                     do 140 i = 1, n
                        acor(i) = acor(i)*r
  140                continue
c     ...............exit
                     go to 690
  150             continue
c                 ------------------------------------------------------
c                  this section computes the predicted values by
c                  effectively multiplying the yh array by the pascal
c                  triangle matrix.  rc is the ratio of new to old
c                  values of the coefficient  h*el(1).  when rc differs
c                  from 1 by more than 30 percent, ipup is set to miter
c                  to force dpjac to be called, if a jacobian is
c                  involved.  in any case, dpjac is called at least
c                  every 20-th step.
c                 ------------------------------------------------------
  160             continue
  170             continue
c                    begin block permitting ...exits to 610
c                       begin block permitting ...exits to 490
                           if (abs(rc-1.0d0) .gt. 0.3d0) ipup = miter
                           if (nst .ge. nstepj + 20) ipup = miter
                           tn = tn + h
                           i1 = nqnyh + 1
                           do 190 jb = 1, nq
                              i1 = i1 - nyh
                              do 180 i = i1, nqnyh
                                 yh1(i) = yh1(i) + yh1(i+nyh)
  180                         continue
  190                      continue
                           ksteps = ksteps + 1
c                          ---------------------------------------------
c                           up to 3 corrector iterations are taken.  a
c                           convergence test is made on the r.m.s. norm
c                           of each correction, weighted by the error
c                           weight vector ewt.  the sum of the
c                           corrections is accumulated in the vector
c                           acor(i).  the yh array is not altered in the
c                           corrector loop.
c                          ---------------------------------------------
  200                      continue
                              m = 0
                              do 210 i = 1, n
                                 y(i) = yh(i,1)
  210                         continue
                              call df(tn,y,savf,rpar,ipar)
                              nfe = nfe + 1
                              if (ipup .le. 0) go to 220
c                                ---------------------------------------
c                                 if indicated, the matrix p = i -
c                                 h*el(1)*j is reevaluated and
c                                 preprocessed before starting the
c                                 corrector iteration.  ipup is set to 0
c                                 as an indicator that this has been
c                                 done.
c                                ---------------------------------------
                                 ipup = 0
                                 rc = 1.0d0
                                 nstepj = nst
                                 crate = 0.7d0
                                 call dpjac(neq,y,yh,nyh,ewt,acor,savf,
     1                                      wm,iwm,df,djac,rpar,ipar)
c                          ......exit
                                 if (ier .ne. 0) go to 440
  220                         continue
                              do 230 i = 1, n
                                 acor(i) = 0.0d0
  230                         continue
  240                         continue
                                 if (miter .ne. 0) go to 270
c                                   ------------------------------------
c                                    in the case of functional
c                                    iteration, update y directly from
c                                    the result of the last function
c                                    evaluation.
c                                   ------------------------------------
                                    do 250 i = 1, n
                                       savf(i) = h*savf(i) - yh(i,2)
                                       y(i) = savf(i) - acor(i)
  250                               continue
                                    del = dvnrms(n,y,ewt)
                                    do 260 i = 1, n
                                       y(i) = yh(i,1) + el(1)*savf(i)
                                       acor(i) = savf(i)
  260                               continue
                                 go to 300
  270                            continue
c                                   ------------------------------------
c                                    in the case of the chord method,
c                                    compute the corrector error, and
c                                    solve the linear system with that
c                                    as right-hand side and p as
c                                    coefficient matrix.
c                                   ------------------------------------
                                    do 280 i = 1, n
                                       y(i) = h*savf(i)
     1                                        - (yh(i,2) + acor(i))
  280                               continue
                                    call dslvs(wm,iwm,y,savf)
c                             ......exit
                                    if (ier .ne. 0) go to 430
                                    del = dvnrms(n,y,ewt)
                                    do 290 i = 1, n
                                       acor(i) = acor(i) + y(i)
                                       y(i) = yh(i,1) + el(1)*acor(i)
  290                               continue
  300                            continue
c                                ---------------------------------------
c                                 test for convergence.  if m.gt.0, an
c                                 estimate of the convergence rate
c                                 constant is stored in crate, and this
c                                 is used in the test.
c                                ---------------------------------------
                                 if (m .ne. 0)
     1                              crate = max(0.2d0*crate,del/delp)
                                 dcon = del*min(1.0d0,1.5d0*crate)
     1                                  /(tesco(2,nq)*conit)
                                 if (dcon .gt. 1.0d0) go to 420
c                                   ------------------------------------
c                                    the corrector has converged.  ipup
c                                    is set to -1 if miter .ne. 0, to
c                                    signal that the jacobian involved
c                                    may need updating later.  the local
c                                    error test is made and control
c                                    passes to statement 500 if it
c                                    fails.
c                                   ------------------------------------
                                    if (miter .ne. 0) ipup = -1
                                    if (m .eq. 0) dsm = del/tesco(2,nq)
                                    if (m .gt. 0)
     1                                 dsm = dvnrms(n,acor,ewt)
     2                                       /tesco(2,nq)
                                    if (dsm .gt. 1.0d0) go to 380
c                                      begin block
c                                      permitting ...exits to 360
c                                         ------------------------------
c                                          after a successful step,
c                                          update the yh array.
c                                          consider changing h if ialth
c                                          = 1.  otherwise decrease
c                                          ialth by 1.  if ialth is then
c                                          1 and nq .lt. maxord, then
c                                          acor is saved for use in a
c                                          possible order increase on
c                                          the next step.  if a change
c                                          in h is considered, an
c                                          increase or decrease in order
c                                          by one is considered also.  a
c                                          change in h is made only if
c                                          it is by a factor of at least
c                                          1.1.  if not, ialth is set to
c                                          3 to prevent testing for that
c                                          many steps.
c                                         ------------------------------
                                          kflag = 0
                                          iredo = 0
                                          nst = nst + 1
                                          hu = h
                                          nqu = nq
                                          do 320 j = 1, l
                                             do 310 i = 1, n
                                                yh(i,j) = yh(i,j)
     1                                                    + el(j)
     2                                                      *acor(i)
  310                                        continue
  320                                     continue
                                          ialth = ialth - 1
                                          if (ialth .ne. 0) go to 340
c                                            ---------------------------
c                                             regardless of the success
c                                             or failure of the step,
c                                             factors rhdn, rhsm, and
c                                             rhup are computed, by
c                                             which h could be
c                                             multiplied at order nq -
c                                             1, order nq, or order nq +
c                                             1, respectively.  in the
c                                             case of failure, rhup =
c                                             0.0 to avoid an order
c                                             increase.  the largest of
c                                             these is determined and
c                                             the new order chosen
c                                             accordingly.  if the order
c                                             is to be increased, we
c                                             compute one additional
c                                             scaled derivative.
c                                            ---------------------------
                                             rhup = 0.0d0
c                       .....................exit
                                             if (l .eq. lmax) go to 490
                                             do 330 i = 1, n
                                                savf(i) = acor(i)
     1                                                    - yh(i,lmax)
  330                                        continue
                                             dup = dvnrms(n,savf,ewt)
     1                                             /tesco(3,nq)
                                             exup = 1.0d0/(l+1)
                                             rhup = 1.0d0
     1                                              /(1.4d0*dup**exup
     2                                                + 0.0000014d0)
c                       .....................exit
                                             go to 490
  340                                     continue
c                                      ...exit
                                          if (ialth .gt. 1) go to 360
c                                      ...exit
                                          if (l .eq. lmax) go to 360
                                          do 350 i = 1, n
                                             yh(i,lmax) = acor(i)
  350                                     continue
  360                                  continue
                                       r = 1.0d0/tesco(2,nqu)
                                       do 370 i = 1, n
                                          acor(i) = acor(i)*r
  370                                  continue
c     .................................exit
                                       go to 690
  380                               continue
c                                   ------------------------------------
c                                    the error test failed.  kflag keeps
c                                    track of multiple failures.
c                                    restore tn and the yh array to
c                                    their previous values, and prepare
c                                    to try the step again.  compute the
c                                    optimum step size for this or one
c                                    lower order.  after 2 or more
c                                    failures, h is forced to decrease
c                                    by a factor of 0.2 or less.
c                                   ------------------------------------
                                    kflag = kflag - 1
                                    tn = told
                                    i1 = nqnyh + 1
                                    do 400 jb = 1, nq
                                       i1 = i1 - nyh
                                       do 390 i = i1, nqnyh
                                          yh1(i) = yh1(i) - yh1(i+nyh)
  390                                  continue
  400                               continue
                                    rmax = 2.0d0
                                    if (abs(h) .gt. hmin*1.00001d0)
     1                                 go to 410
c                                      ---------------------------------
c                                       all returns are made through
c                                       this section.  h is saved in
c                                       hold to allow the caller to
c                                       change h on the next step.
c                                      ---------------------------------
                                       kflag = -1
c     .................................exit
                                       go to 690
  410                               continue
c                    ...............exit
                                    if (kflag .le. -3) go to 610
                                    iredo = 2
                                    rhup = 0.0d0
c                       ............exit
                                    go to 490
  420                            continue
                                 m = m + 1
c                             ...exit
                                 if (m .eq. 3) go to 430
c                             ...exit
                                 if (m .ge. 2 .and. del .gt. 2.0d0*delp)
     1                              go to 430
                                 delp = del
                                 call df(tn,y,savf,rpar,ipar)
                                 nfe = nfe + 1
                              go to 240
  430                         continue
c                             ------------------------------------------
c                              the corrector iteration failed to
c                              converge in 3 tries.  if miter .ne. 0 and
c                              the jacobian is out of date, dpjac is
c                              called for the next try.  otherwise the
c                              yh array is retracted to its values
c                              before prediction, and h is reduced, if
c                              possible.  if h cannot be reduced or 10
c                              failures have occurred, exit with kflag =
c                              -2.
c                             ------------------------------------------
c                          ...exit
                              if (ipup .eq. 0) go to 440
                              ipup = miter
                           go to 200
  440                      continue
                           tn = told
                           ncf = ncf + 1
                           rmax = 2.0d0
                           i1 = nqnyh + 1
                           do 460 jb = 1, nq
                              i1 = i1 - nyh
                              do 450 i = i1, nqnyh
                                 yh1(i) = yh1(i) - yh1(i+nyh)
  450                         continue
  460                      continue
                           if (abs(h) .gt. hmin*1.00001d0) go to 470
                              kflag = -2
c     ........................exit
                              go to 690
  470                      continue
                           if (ncf .ne. 10) go to 480
                              kflag = -2
c     ........................exit
                              go to 690
  480                      continue
                           rh = 0.25d0
                           ipup = miter
                           iredo = 1
c                 .........exit
                           go to 650
  490                   continue
                        exsm = 1.0d0/l
                        rhsm = 1.0d0/(1.2d0*dsm**exsm + 0.0000012d0)
                        rhdn = 0.0d0
                        if (nq .eq. 1) go to 500
                           ddn = dvnrms(n,yh(1,l),ewt)/tesco(1,nq)
                           exdn = 1.0d0/nq
                           rhdn = 1.0d0/(1.3d0*ddn**exdn + 0.0000013d0)
  500                   continue
                        if (rhsm .ge. rhup) go to 550
                           if (rhup .le. rhdn) go to 540
                              newq = l
                              rh = rhup
                              if (rh .ge. 1.1d0) go to 520
                                 ialth = 3
                                 r = 1.0d0/tesco(2,nqu)
                                 do 510 i = 1, n
                                    acor(i) = acor(i)*r
  510                            continue
c     ...........................exit
                                 go to 690
  520                         continue
                              r = el(l)/l
                              do 530 i = 1, n
                                 yh(i,newq+1) = acor(i)*r
  530                         continue
                              nq = newq
                              l = nq + 1
                              iret = 2
c           ..................exit
                              go to 680
  540                      continue
                        go to 580
  550                   continue
                        if (rhsm .lt. rhdn) go to 580
                           newq = nq
                           rh = rhsm
                           if (kflag .eq. 0 .and. rh .lt. 1.1d0)
     1                        go to 560
                              if (kflag .le. -2) rh = min(rh,0.2d0)
c                             ------------------------------------------
c                              if there is a change of order, reset nq,
c                              l, and the coefficients.  in any case h
c                              is reset according to rh and the yh array
c                              is rescaled.  then exit from 680 if the
c                              step was ok, or redo the step otherwise.
c                             ------------------------------------------
c                 ............exit
                              if (newq .eq. nq) go to 650
                              nq = newq
                              l = nq + 1
                              iret = 2
c           ..................exit
                              go to 680
  560                      continue
                           ialth = 3
                           r = 1.0d0/tesco(2,nqu)
                           do 570 i = 1, n
                              acor(i) = acor(i)*r
  570                      continue
c     .....................exit
                           go to 690
  580                   continue
                        newq = nq - 1
                        rh = rhdn
                        if (kflag .lt. 0 .and. rh .gt. 1.0d0) rh = 1.0d0
                        if (kflag .eq. 0 .and. rh .lt. 1.1d0) go to 590
                           if (kflag .le. -2) rh = min(rh,0.2d0)
c                          ---------------------------------------------
c                           if there is a change of order, reset nq, l,
c                           and the coefficients.  in any case h is
c                           reset according to rh and the yh array is
c                           rescaled.  then exit from 680 if the step
c                           was ok, or redo the step otherwise.
c                          ---------------------------------------------
c                 .........exit
                           if (newq .eq. nq) go to 650
                           nq = newq
                           l = nq + 1
                           iret = 2
c           ...............exit
                           go to 680
  590                   continue
                        ialth = 3
                        r = 1.0d0/tesco(2,nqu)
                        do 600 i = 1, n
                           acor(i) = acor(i)*r
  600                   continue
c     ..................exit
                        go to 690
  610                continue
c                    ---------------------------------------------------
c                     control reaches this section if 3 or more failures
c                     have occurred.  if 10 failures have occurred, exit
c                     with kflag = -1.  it is assumed that the
c                     derivatives that have accumulated in the yh array
c                     have errors of the wrong order.  hence the first
c                     derivative is recomputed, and the order is set to
c                     1.  then h is reduced by a factor of 10, and the
c                     step is retried, until it succeeds or h reaches
c                     hmin.
c                    ---------------------------------------------------
                     if (kflag .ne. -10) go to 620
c                       ------------------------------------------------
c                        all returns are made through this section.  h
c                        is saved in hold to allow the caller to change
c                        h on the next step.
c                       ------------------------------------------------
                        kflag = -1
c     ..................exit
                        go to 690
  620                continue
                     rh = 0.1d0
                     rh = max(hmin/abs(h),rh)
                     h = h*rh
                     do 630 i = 1, n
                        y(i) = yh(i,1)
  630                continue
                     call df(tn,y,savf,rpar,ipar)
                     nfe = nfe + 1
                     do 640 i = 1, n
                        yh(i,2) = h*savf(i)
  640                continue
                     ipup = miter
                     ialth = 5
c              ......exit
                     if (nq .ne. 1) go to 670
                  go to 170
  650             continue
  660             continue
                  rh = max(rh,hmin/abs(h))
               go to 110
  670          continue
               nq = 1
               l = 2
               iret = 3
  680       continue
         go to 70
  690 continue
      hold = h
      jstart = 1
      return
c     ----------------------- end of subroutine dstod
c     -----------------------
      end
