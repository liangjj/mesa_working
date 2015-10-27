*deck sdaini
      subroutine sdaini (x, y, yprime, neq, res, jac, h, wt, idid, rpar,
     *   ipar, phi, delta, e, wm, iwm, hmin, uround, nonneg, ntemp)
c***begin prologue  sdaini
c***subsidiary
c***purpose  initialization routine for sdassl.
c***library   slatec (dassl)
c***type      single precision (sdaini-s, ddaini-d)
c***author  petzold, linda r., (llnl)
c***description
c-----------------------------------------------------------------
c     sdaini takes one step of size h or smaller
c     with the backward euler method, to
c     find yprime.  x and y are updated to be consistent with the
c     new step.  a modified damped newton iteration is used to
c     solve the corrector iteration.
c
c     the initial guess for yprime is used in the
c     prediction, and in forming the iteration
c     matrix, but is not involved in the
c     error test. this may have trouble
c     converging if the initial guess is no
c     good, or if g(x,y,yprime) depends
c     nonlinearly on yprime.
c
c     the parameters represent:
c     x --         independent variable
c     y --         solution vector at x
c     yprime --    derivative of solution vector
c     neq --       number of equations
c     h --         stepsize. imder may use a stepsize
c                  smaller than h.
c     wt --        vector of weights for error
c                  criterion
c     idid --      completion code with the following meanings
c                  idid= 1 -- yprime was found successfully
c                  idid=-12 -- sdaini failed to find yprime
c     rpar,ipar -- real and integer parameter arrays
c                  that are not altered by sdaini
c     phi --       work space for sdaini
c     delta,e --   work space for sdaini
c     wm,iwm --    real and integer arrays storing
c                  matrix information
c
c-----------------------------------------------------------------
c***routines called  sdajac, sdanrm, sdaslv
c***revision history  (yymmdd)
c   830315  date written
c   901009  finished conversion to slatec 4.0 format (f.n.fritsch)
c   901019  merged changes made by c. ulrich with slatec 4.0 format.
c   901026  added explicit declarations for all variables and minor
c           cosmetic changes to prologue.  (fnf)
c   901030  minor corrections to declarations.  (fnf)
c***end prologue  sdaini
c
      integer  neq, idid, ipar(*), iwm(*), nonneg, ntemp
      real  x, y(*), yprime(*), h, wt(*), rpar(*), phi(neq,*), delta(*),
     *   e(*), wm(*), hmin, uround
      external  res, jac
c
      external  sdajac, sdanrm, sdaslv
      real  sdanrm
c
      integer  i, ier, ires, jcalc, lnje, lnre, m, maxit, mjac, ncf,
     *   nef, nsf
      real  cj, damp, delnrm, err, oldnrm, r, rate, s, xold, ynorm
      logical  convgd
c
      parameter (lnre=12)
      parameter (lnje=13)
c
      data maxit/10/,mjac/5/
      data damp/0.75e0/
c
c
c---------------------------------------------------
c     block 1.
c     initializations.
c---------------------------------------------------
c
c***first executable statement  sdaini
      idid=1
      nef=0
      ncf=0
      nsf=0
      xold=x
      ynorm=sdanrm(neq,y,wt,rpar,ipar)
c
c     save y and yprime in phi
      do 100 i=1,neq
         phi(i,1)=y(i)
100      phi(i,2)=yprime(i)
c
c
c----------------------------------------------------
c     block 2.
c     do one backward euler step.
c----------------------------------------------------
c
c     set up for start of corrector iteration
200   cj=1.0e0/h
      x=x+h
c
c     predict solution and derivative
      do 250 i=1,neq
250     y(i)=y(i)+h*yprime(i)
c
      jcalc=-1
      m=0
      convgd=.true.
c
c
c     corrector loop.
300   iwm(lnre)=iwm(lnre)+1
      ires=0
c
      call res(x,y,yprime,delta,ires,rpar,ipar)
      if (ires.lt.0) go to 430
c
c
c     evaluate the iteration matrix
      if (jcalc.ne.-1) go to 310
      iwm(lnje)=iwm(lnje)+1
      jcalc=0
      call sdajac(neq,x,y,yprime,delta,cj,h,
     *   ier,wt,e,wm,iwm,res,ires,
     *   uround,jac,rpar,ipar,ntemp)
c
      s=1000000.e0
      if (ires.lt.0) go to 430
      if (ier.ne.0) go to 430
      nsf=0
c
c
c
c     multiply residual by damping factor
310   continue
      do 320 i=1,neq
320      delta(i)=delta(i)*damp
c
c     compute a new iterate (back substitution)
c     store the correction in delta
c
      call sdaslv(neq,delta,wm,iwm)
c
c     update y and yprime
      do 330 i=1,neq
         y(i)=y(i)-delta(i)
330      yprime(i)=yprime(i)-cj*delta(i)
c
c     test for convergence of the iteration.
c
      delnrm=sdanrm(neq,delta,wt,rpar,ipar)
      if (delnrm.le.100.e0*uround*ynorm)
     *   go to 400
c
      if (m.gt.0) go to 340
         oldnrm=delnrm
         go to 350
c
340   rate=(delnrm/oldnrm)**(1.0e0/m)
      if (rate.gt.0.90e0) go to 430
      s=rate/(1.0e0-rate)
c
350   if (s*delnrm .le. 0.33e0) go to 400
c
c
c     the corrector has not yet converged. update
c     m and and test whether the maximum
c     number of iterations have been tried.
c     every mjac iterations, get a new
c     iteration matrix.
c
      m=m+1
      if (m.ge.maxit) go to 430
c
      if ((m/mjac)*mjac.eq.m) jcalc=-1
      go to 300
c
c
c     the iteration has converged.
c     check nonnegativity constraints
400   if (nonneg.eq.0) go to 450
      do 410 i=1,neq
410      delta(i)=min(y(i),0.0e0)
c
      delnrm=sdanrm(neq,delta,wt,rpar,ipar)
      if (delnrm.gt.0.33e0) go to 430
c
      do 420 i=1,neq
         y(i)=y(i)-delta(i)
420      yprime(i)=yprime(i)-cj*delta(i)
      go to 450
c
c
c     exits from corrector loop.
430   convgd=.false.
450   if (.not.convgd) go to 600
c
c
c
c-----------------------------------------------------
c     block 3.
c     the corrector iteration converged.
c     do error test.
c-----------------------------------------------------
c
      do 510 i=1,neq
510      e(i)=y(i)-phi(i,1)
      err=sdanrm(neq,e,wt,rpar,ipar)
c
      if (err.le.1.0e0) return
c
c
c
c--------------------------------------------------------
c     block 4.
c     the backward euler step failed. restore x, y
c     and yprime to their original values.
c     reduce stepsize and try again, if
c     possible.
c---------------------------------------------------------
c
600   continue
      x = xold
      do 610 i=1,neq
         y(i)=phi(i,1)
610      yprime(i)=phi(i,2)
c
      if (convgd) go to 640
      if (ier.eq.0) go to 620
         nsf=nsf+1
         h=h*0.25e0
         if (nsf.lt.3.and.abs(h).ge.hmin) go to 690
         idid=-12
         return
620   if (ires.gt.-2) go to 630
         idid=-12
         return
630   ncf=ncf+1
      h=h*0.25e0
      if (ncf.lt.10.and.abs(h).ge.hmin) go to 690
         idid=-12
         return
c
640   nef=nef+1
      r=0.90e0/(2.0e0*err+0.0001e0)
      r=max(0.1e0,min(0.5e0,r))
      h=h*r
      if (abs(h).ge.hmin.and.nef.lt.10) go to 690
         idid=-12
         return
690      go to 200
c
c-------------end of subroutine sdaini----------------------
      end
