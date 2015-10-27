*deck ddastp
      subroutine ddastp (x, y, yprime, neq, res, jac, h, wt, jstart,
     *   idid, rpar, ipar, phi, delta, e, wm, iwm, alpha, beta, gamma,
     *   psi, sigma, cj, cjold, hold, s, hmin, uround, iphase, jcalc, k,
     *   kold, ns, nonneg, ntemp)
c***begin prologue  ddastp
c***subsidiary
c***purpose  perform one step of the ddassl integration.
c***library   slatec (dassl)
c***type      double precision (sdastp-s, ddastp-d)
c***author  petzold, linda r., (llnl)
c***description
c-----------------------------------------------------------------------
c     ddastp solves a system of differential/
c     algebraic equations of the form
c     g(x,y,yprime) = 0,  for one step (normally
c     from x to x+h).
c
c     the methods used are modified divided
c     difference,fixed leading coefficient
c     forms of backward differentiation
c     formulas. the code adjusts the stepsize
c     and order to control the local error per
c     step.
c
c
c     the parameters represent
c     x  --        independent variable
c     y  --        solution vector at x
c     yprime --    derivative of solution vector
c                  after successful step
c     neq --       number of equations to be integrated
c     res --       external user-supplied subroutine
c                  to evaluate the residual.  the call is
c                  call res(x,y,yprime,delta,ires,rpar,ipar)
c                  x,y,yprime are input.  delta is output.
c                  on input, ires=0.  res should alter ires only
c                  if it encounters an illegal value of y or a
c                  stop condition.  set ires=-1 if an input value
c                  of y is illegal, and ddastp will try to solve
c                  the problem without getting ires = -1.  if
c                  ires=-2, ddastp returns control to the calling
c                  program with idid = -11.
c     jac --       external user-supplied routine to evaluate
c                  the iteration matrix (this is optional)
c                  the call is of the form
c                  call jac(x,y,yprime,pd,cj,rpar,ipar)
c                  pd is the matrix of partial derivatives,
c                  pd=dg/dy+cj*dg/dyprime
c     h --         appropriate step size for next step.
c                  normally determined by the code
c     wt --        vector of weights for error criterion.
c     jstart --    integer variable set 0 for
c                  first step, 1 otherwise.
c     idid --      completion code with the following meanings:
c                  idid= 1 -- the step was completed successfully
c                  idid=-6 -- the error test failed repeatedly
c                  idid=-7 -- the corrector could not converge
c                  idid=-8 -- the iteration matrix is singular
c                  idid=-9 -- the corrector could not converge.
c                             there were repeated error test
c                             failures on this step.
c                  idid=-10-- the corrector could not converge
c                             because ires was equal to minus one
c                  idid=-11-- ires equal to -2 was encountered,
c                             and control is being returned to
c                             the calling program
c     rpar,ipar -- real and integer parameter arrays that
c                  are used for communication between the
c                  calling program and external user routines
c                  they are not altered by ddastp
c     phi --       array of divided differences used by
c                  ddastp. the length is neq*(k+1),where
c                  k is the maximum order
c     delta,e --   work vectors for ddastp of length neq
c     wm,iwm --    real and integer arrays storing
c                  matrix information such as the matrix
c                  of partial derivatives,permutation
c                  vector, and various other information.
c
c     the other parameters are information
c     which is needed internally by ddastp to
c     continue from step to step.
c
c-----------------------------------------------------------------------
c***routines called  ddajac, ddanrm, ddaslv, ddatrp
c***revision history  (yymmdd)
c   830315  date written
c   901009  finished conversion to slatec 4.0 format (f.n.fritsch)
c   901019  merged changes made by c. ulrich with slatec 4.0 format.
c   901026  added explicit declarations for all variables and minor
c           cosmetic changes to prologue.  (fnf)
c***end prologue  ddastp
c
      integer  neq, jstart, idid, ipar(*), iwm(*), iphase, jcalc, k,
     *   kold, ns, nonneg, ntemp
      double precision
     *   x, y(*), yprime(*), h, wt(*), rpar(*), phi(neq,*), delta(*),
     *   e(*), wm(*), alpha(*), beta(*), gamma(*), psi(*), sigma(*), cj,
     *   cjold, hold, s, hmin, uround
      external  res, jac
c
      external  ddajac, ddanrm, ddaslv, ddatrp
      double precision  ddanrm
c
      integer  i, ier, ires, j, j1, kdiff, km1, knew, kp1, kp2, lctf,
     *   letf, lmxord, lnje, lnre, lnst, m, maxit, ncf, nef, nsf, nsp1
      double precision
     *   alpha0, alphas, cjlast, ck, delnrm, enorm, erk, erkm1,
     *   erkm2, erkp1, err, est, hnew, oldnrm, pnorm, r, rate, temp1,
     *   temp2, terk, terkm1, terkm2, terkp1, xold, xrate
      logical  convgd
c
      parameter (lmxord=3)
      parameter (lnst=11)
      parameter (lnre=12)
      parameter (lnje=13)
      parameter (letf=14)
      parameter (lctf=15)
c
      data maxit/4/
      data xrate/0.25d0/
c
c
c
c
c
c-----------------------------------------------------------------------
c     block 1.
c     initialize. on the first call,set
c     the order to 1 and initialize
c     other variables.
c-----------------------------------------------------------------------
c
c     initializations for all calls
c***first executable statement  ddastp
      idid=1
      xold=x
      ncf=0
      nsf=0
      nef=0
      if(jstart .ne. 0) go to 120
c
c     if this is the first step,perform
c     other initializations
      iwm(letf) = 0
      iwm(lctf) = 0
      k=1
      kold=0
      hold=0.0d0
      jstart=1
      psi(1)=h
      cjold = 1.0d0/h
      cj = cjold
      s = 100.d0
      jcalc = -1
      delnrm=1.0d0
      iphase = 0
      ns=0
120   continue
c
c
c
c
c
c-----------------------------------------------------------------------
c     block 2
c     compute coefficients of formulas for
c     this step.
c-----------------------------------------------------------------------
200   continue
      kp1=k+1
      kp2=k+2
      km1=k-1
      xold=x
      if(h.ne.hold.or.k .ne. kold) ns = 0
      ns=min(ns+1,kold+2)
      nsp1=ns+1
      if(kp1 .lt. ns)go to 230
c
      beta(1)=1.0d0
      alpha(1)=1.0d0
      temp1=h
      gamma(1)=0.0d0
      sigma(1)=1.0d0
      do 210 i=2,kp1
         temp2=psi(i-1)
         psi(i-1)=temp1
         beta(i)=beta(i-1)*psi(i-1)/temp2
         temp1=temp2+h
         alpha(i)=h/temp1
         sigma(i)=(i-1)*sigma(i-1)*alpha(i)
         gamma(i)=gamma(i-1)+alpha(i-1)/h
210      continue
      psi(kp1)=temp1
230   continue
c
c     compute alphas, alpha0
      alphas = 0.0d0
      alpha0 = 0.0d0
      do 240 i = 1,k
        alphas = alphas - 1.0d0/i
        alpha0 = alpha0 - alpha(i)
240     continue
c
c     compute leading coefficient cj
      cjlast = cj
      cj = -alphas/h
c
c     compute variable stepsize error coefficient ck
      ck = abs(alpha(kp1) + alphas - alpha0)
      ck = max(ck,alpha(kp1))
c
c     decide whether new jacobian is needed
      temp1 = (1.0d0 - xrate)/(1.0d0 + xrate)
      temp2 = 1.0d0/temp1
      if (cj/cjold .lt. temp1 .or. cj/cjold .gt. temp2) jcalc = -1
      if (cj .ne. cjlast) s = 100.d0
c
c     change phi to phi star
      if(kp1 .lt. nsp1) go to 280
      do 270 j=nsp1,kp1
         do 260 i=1,neq
260         phi(i,j)=beta(j)*phi(i,j)
270      continue
280   continue
c
c     update time
      x=x+h
c
c
c
c
c
c-----------------------------------------------------------------------
c     block 3
c     predict the solution and derivative,
c     and solve the corrector equation
c-----------------------------------------------------------------------
c
c     first,predict the solution and derivative
300   continue
      do 310 i=1,neq
         y(i)=phi(i,1)
310      yprime(i)=0.0d0
      do 330 j=2,kp1
         do 320 i=1,neq
            y(i)=y(i)+phi(i,j)
320         yprime(i)=yprime(i)+gamma(j)*phi(i,j)
330   continue
      pnorm = ddanrm (neq,y,wt,rpar,ipar)
c
c
c
c     solve the corrector equation using a
c     modified newton scheme.
      convgd= .true.
      m=0
      iwm(lnre)=iwm(lnre)+1
      ires = 0
      call res(x,y,yprime,delta,ires,rpar,ipar)
      if (ires .lt. 0) go to 380
c
c
c     if indicated,reevaluate the
c     iteration matrix pd = dg/dy + cj*dg/dyprime
c     (where g(x,y,yprime)=0). set
c     jcalc to 0 as an indicator that
c     this has been done.
      if(jcalc .ne. -1)go to 340
      iwm(lnje)=iwm(lnje)+1
      jcalc=0
      call ddajac(neq,x,y,yprime,delta,cj,h,
     * ier,wt,e,wm,iwm,res,ires,uround,jac,rpar,
     * ipar,ntemp)
      cjold=cj
      s = 100.d0
      if (ires .lt. 0) go to 380
      if(ier .ne. 0)go to 380
      nsf=0
c
c
c     initialize the error accumulation vector e.
340   continue
      do 345 i=1,neq
345      e(i)=0.0d0
c
c
c     corrector loop.
350   continue
c
c     multiply residual by temp1 to accelerate convergence
      temp1 = 2.0d0/(1.0d0 + cj/cjold)
      do 355 i = 1,neq
355     delta(i) = delta(i) * temp1
c
c     compute a new iterate (back-substitution).
c     store the correction in delta.
      call ddaslv(neq,delta,wm,iwm)
c
c     update y, e, and yprime
      do 360 i=1,neq
         y(i)=y(i)-delta(i)
         e(i)=e(i)-delta(i)
360      yprime(i)=yprime(i)-cj*delta(i)
c
c     test for convergence of the iteration
      delnrm=ddanrm(neq,delta,wt,rpar,ipar)
      if (delnrm .le. 100.d0*uround*pnorm) go to 375
      if (m .gt. 0) go to 365
         oldnrm = delnrm
         go to 367
365   rate = (delnrm/oldnrm)**(1.0d0/m)
      if (rate .gt. 0.90d0) go to 370
      s = rate/(1.0d0 - rate)
367   if (s*delnrm .le. 0.33d0) go to 375
c
c     the corrector has not yet converged.
c     update m and test whether the
c     maximum number of iterations have
c     been tried.
      m=m+1
      if(m.ge.maxit)go to 370
c
c     evaluate the residual
c     and go back to do another iteration
      iwm(lnre)=iwm(lnre)+1
      ires = 0
      call res(x,y,yprime,delta,ires,
     *  rpar,ipar)
      if (ires .lt. 0) go to 380
      go to 350
c
c
c     the corrector failed to converge in maxit
c     iterations. if the iteration matrix
c     is not current,re-do the step with
c     a new iteration matrix.
370   continue
      if(jcalc.eq.0)go to 380
      jcalc=-1
      go to 300
c
c
c     the iteration has converged.  if nonnegativity of solution is
c     required, set the solution nonnegative, if the perturbation
c     to do it is small enough.  if the change is too large, then
c     consider the corrector iteration to have failed.
375   if(nonneg .eq. 0) go to 390
      do 377 i = 1,neq
377      delta(i) = min(y(i),0.0d0)
      delnrm = ddanrm(neq,delta,wt,rpar,ipar)
      if(delnrm .gt. 0.33d0) go to 380
      do 378 i = 1,neq
378      e(i) = e(i) - delta(i)
      go to 390
c
c
c     exits from block 3
c     no convergence with current iteration
c     matrix,or singular iteration matrix
380   convgd= .false.
390   jcalc = 1
      if(.not.convgd)go to 600
c
c
c
c
c
c-----------------------------------------------------------------------
c     block 4
c     estimate the errors at orders k,k-1,k-2
c     as if constant stepsize was used. estimate
c     the local error at order k and test
c     whether the current step is successful.
c-----------------------------------------------------------------------
c
c     estimate errors at orders k,k-1,k-2
      enorm = ddanrm(neq,e,wt,rpar,ipar)
      erk = sigma(k+1)*enorm
      terk = (k+1)*erk
      est = erk
      knew=k
      if(k .eq. 1)go to 430
      do 405 i = 1,neq
405     delta(i) = phi(i,kp1) + e(i)
      erkm1=sigma(k)*ddanrm(neq,delta,wt,rpar,ipar)
      terkm1 = k*erkm1
      if(k .gt. 2)go to 410
      if(terkm1 .le. 0.5d0*terk)go to 420
      go to 430
410   continue
      do 415 i = 1,neq
415     delta(i) = phi(i,k) + delta(i)
      erkm2=sigma(k-1)*ddanrm(neq,delta,wt,rpar,ipar)
      terkm2 = (k-1)*erkm2
      if(max(terkm1,terkm2).gt.terk)go to 430
c     lower the order
420   continue
      knew=k-1
      est = erkm1
c
c
c     calculate the local error for the current step
c     to see if the step was successful
430   continue
      err = ck * enorm
      if(err .gt. 1.0d0)go to 600
c
c
c
c
c
c-----------------------------------------------------------------------
c     block 5
c     the step is successful. determine
c     the best order and stepsize for
c     the next step. update the differences
c     for the next step.
c-----------------------------------------------------------------------
      idid=1
      iwm(lnst)=iwm(lnst)+1
      kdiff=k-kold
      kold=k
      hold=h
c
c
c     estimate the error at order k+1 unless:
c        already decided to lower order, or
c        already using maximum order, or
c        stepsize not constant, or
c        order raised in previous step
      if(knew.eq.km1.or.k.eq.iwm(lmxord))iphase=1
      if(iphase .eq. 0)go to 545
      if(knew.eq.km1)go to 540
      if(k.eq.iwm(lmxord)) go to 550
      if(kp1.ge.ns.or.kdiff.eq.1)go to 550
      do 510 i=1,neq
510      delta(i)=e(i)-phi(i,kp2)
      erkp1 = (1.0d0/(k+2))*ddanrm(neq,delta,wt,rpar,ipar)
      terkp1 = (k+2)*erkp1
      if(k.gt.1)go to 520
      if(terkp1.ge.0.5d0*terk)go to 550
      go to 530
520   if(terkm1.le.min(terk,terkp1))go to 540
      if(terkp1.ge.terk.or.k.eq.iwm(lmxord))go to 550
c
c     raise order
530   k=kp1
      est = erkp1
      go to 550
c
c     lower order
540   k=km1
      est = erkm1
      go to 550
c
c     if iphase = 0, increase order by one and multiply stepsize by
c     factor two
545   k = kp1
      hnew = h*2.0d0
      h = hnew
      go to 575
c
c
c     determine the appropriate stepsize for
c     the next step.
550   hnew=h
      temp2=k+1
      r=(2.0d0*est+0.0001d0)**(-1.0d0/temp2)
      if(r .lt. 2.0d0) go to 555
      hnew = 2.0d0*h
      go to 560
555   if(r .gt. 1.0d0) go to 560
      r = max(0.5d0,min(0.9d0,r))
      hnew = h*r
560   h=hnew
c
c
c     update differences for next step
575   continue
      if(kold.eq.iwm(lmxord))go to 585
      do 580 i=1,neq
580      phi(i,kp2)=e(i)
585   continue
      do 590 i=1,neq
590      phi(i,kp1)=phi(i,kp1)+e(i)
      do 595 j1=2,kp1
         j=kp1-j1+1
         do 595 i=1,neq
595      phi(i,j)=phi(i,j)+phi(i,j+1)
      return
c
c
c
c
c
c-----------------------------------------------------------------------
c     block 6
c     the step is unsuccessful. restore x,psi,phi
c     determine appropriate stepsize for
c     continuing the integration, or exit with
c     an error flag if there have been many
c     failures.
c-----------------------------------------------------------------------
600   iphase = 1
c
c     restore x,phi,psi
      x=xold
      if(kp1.lt.nsp1)go to 630
      do 620 j=nsp1,kp1
         temp1=1.0d0/beta(j)
         do 610 i=1,neq
610         phi(i,j)=temp1*phi(i,j)
620      continue
630   continue
      do 640 i=2,kp1
640      psi(i-1)=psi(i)-h
c
c
c     test whether failure is due to corrector iteration
c     or error test
      if(convgd)go to 660
      iwm(lctf)=iwm(lctf)+1
c
c
c     the newton iteration failed to converge with
c     a current iteration matrix.  determine the cause
c     of the failure and take appropriate action.
      if(ier.eq.0)go to 650
c
c     the iteration matrix is singular. reduce
c     the stepsize by a factor of 4. if
c     this happens three times in a row on
c     the same step, return with an error flag
      nsf=nsf+1
      r = 0.25d0
      h=h*r
      if (nsf .lt. 3 .and. abs(h) .ge. hmin) go to 690
      idid=-8
      go to 675
c
c
c     the newton iteration failed to converge for a reason
c     other than a singular iteration matrix.  if ires = -2, then
c     return.  otherwise, reduce the stepsize and try again, unless
c     too many failures have occurred.
650   continue
      if (ires .gt. -2) go to 655
      idid = -11
      go to 675
655   ncf = ncf + 1
      r = 0.25d0
      h = h*r
      if (ncf .lt. 10 .and. abs(h) .ge. hmin) go to 690
      idid = -7
      if (ires .lt. 0) idid = -10
      if (nef .ge. 3) idid = -9
      go to 675
c
c
c     the newton scheme converged, and the cause
c     of the failure was the error estimate
c     exceeding the tolerance.
660   nef=nef+1
      iwm(letf)=iwm(letf)+1
      if (nef .gt. 1) go to 665
c
c     on first error test failure, keep current order or lower
c     order by one.  compute new stepsize based on differences
c     of the solution.
      k = knew
      temp2 = k + 1
      r = 0.90d0*(2.0d0*est+0.0001d0)**(-1.0d0/temp2)
      r = max(0.25d0,min(0.9d0,r))
      h = h*r
      if (abs(h) .ge. hmin) go to 690
      idid = -6
      go to 675
c
c     on second error test failure, use the current order or
c     decrease order by one.  reduce the stepsize by a factor of
c     four.
665   if (nef .gt. 2) go to 670
      k = knew
      h = 0.25d0*h
      if (abs(h) .ge. hmin) go to 690
      idid = -6
      go to 675
c
c     on third and subsequent error test failures, set the order to
c     one and reduce the stepsize by a factor of four.
670   k = 1
      h = 0.25d0*h
      if (abs(h) .ge. hmin) go to 690
      idid = -6
      go to 675
c
c
c
c
c     for all crashes, restore y to its last value,
c     interpolate to find yprime at last x, and return
675   continue
      call ddatrp(x,x,y,yprime,neq,k,phi,psi)
      return
c
c
c     go back and try this step again
690   go to 200
c
c------end of subroutine ddastp------
      end
