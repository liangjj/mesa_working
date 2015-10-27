*deck derkfs
      subroutine derkfs (f, neq, t, y, tout, info, rtol, atol, idid, h,
     +   tolfac, yp, f1, f2, f3, f4, f5, ys, told, dtsign, u26, rer,
     +   init, ksteps, kop, iquit, stiff, nonstf, ntstep, nstifs, rpar,
     +   ipar)
c***begin prologue  derkfs
c***subsidiary
c***purpose  subsidiary to derkf
c***library   slatec
c***type      single precision (derkfs-s, drkfs-d)
c***author  watts, h. a., (snla)
c***description
c
c     fehlberg fourth-fifth order runge-kutta method
c **********************************************************************
c
c     derkfs integrates a system of first order ordinary differential
c     equations as described in the comments for derkf .
c
c     the arrays yp,f1,f2,f3,f4,f5,and ys  (of length at least neq)
c     appear in the call list for variable dimensioning purposes.
c
c     the variables h,tolfac,told,dtsign,u26,rer,init,ksteps,kop,iquit,
c     stiff,nonstf,ntstep, and nstifs are used internally by the code
c     and appear in the call list to eliminate local retention of
c     variables between calls. accordingly, these variables and the
c     array yp should not be altered.
c     items of possible interest are
c         h  - an appropriate step size to be used for the next step
c         tolfac - factor of change in the tolerances
c         yp - derivative of solution vector at t
c         ksteps - counter on the number of steps attempted
c
c **********************************************************************
c
c***see also  derkf
c***routines called  defehl, hstart, hvnrm, r1mach, xermsg
c***revision history  (yymmdd)
c   800501  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   891024  changed references from vnorm to hvnrm.  (wrb)
c   891024  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   900510  convert xerrwv calls to xermsg calls, replace gotos with
c           if-then-elses.  (rwc)
c   910722  updated author section.  (als)
c***end prologue  derkfs
c
      logical hfaild,output,stiff,nonstf
      character*8 xern1
      character*16 xern3, xern4
c
      dimension y(*),yp(*),f1(*),f2(*),f3(*),f4(*),f5(*),
     1          ys(*),info(15),rtol(*),atol(*),rpar(*),ipar(*)
c
      external f
c
c.......................................................................
c
c  a fifth order method will generally not be capable of delivering
c  accuracies near limiting precision on computers with long
c  wordlengths. to protect against limiting precision difficulties
c  arising from unreasonable accuracy requests, an appropriate
c  tolerance threshold remin is assigned for this method. this value
c  should not be changed across different machines.
c
      save remin, mxstep, mxkop
      data remin/1.e-12/
c
c.......................................................................
c
c  the expense of solving the problem is monitored by counting the
c  number of  steps attempted. when this exceeds  mxstep, the counter
c  is reset to zero and the user is informed about possible excessive
c  work.
c
      data mxstep/500/
c
c.......................................................................
c
c  inefficiency caused by too frequent output is monitored by counting
c  the number of step sizes which are severely shortened due solely to
c  the choice of output points. when the number of abuses exceed mxkop,
c  the counter is reset to zero and the user is informed about possible
c  misuse of the code.
c
      data mxkop/100/
c
c.......................................................................
c
c***first executable statement  derkfs
      if (info(1) .eq. 0) then
c
c on the first call , perform initialization --
c        define the machine unit roundoff quantity  u  by calling the
c        function routine  r1mach. the user must make sure that the
c        values set in r1mach are relevant to the computer being used.
c
         u = r1mach(4)
c                       -- set associated machine dependent parameters
         u26 = 26.*u
         rer = 2.*u+remin
c                       -- set termination flag
         iquit = 0
c                       -- set initialization indicator
         init = 0
c                       -- set counter for impact of output points
         kop = 0
c                       -- set counter for attempted steps
         ksteps = 0
c                       -- set indicators for stiffness detection
         stiff = .false.
         nonstf = .false.
c                       -- set step counters for stiffness detection
         ntstep = 0
         nstifs = 0
c                       -- reset info(1) for subsequent calls
         info(1) = 1
      endif
c
c.......................................................................
c
c        check validity of input parameters on each entry
c
      if (info(1) .ne. 0 .and. info(1) .ne. 1) then
         write (xern1, '(i8)') info(1)
         call xermsg ('slatec', 'derkfs',
     *      'in derkf, info(1) must be set to 0 ' //
     *      'for the start of a new problem, and must be set to 1 ' //
     *      'following an interrupted task.  you are attempting to ' //
     *      'continue the integration illegally by calling the code ' //
     *      'with  info(1) = ' // xern1, 3, 1)
         idid = -33
      endif
c
      if (info(2) .ne. 0 .and. info(2) .ne. 1) then
         write (xern1, '(i8)') info(2)
         call xermsg ('slatec', 'derkfs',
     *      'in derkf, info(2) must be 0 or 1 indicating scalar ' //
     *      'and vector error tolerances, respectively.  you have ' //
     *      'called the code with info(2) = ' // xern1, 4, 1)
         idid = -33
      endif
c
      if (info(3) .ne. 0 .and. info(3) .ne. 1) then
         write (xern1, '(i8)') info(3)
         call xermsg ('slatec', 'derkfs',
     *      'in derkf, info(3) must be 0 or 1 indicating the ' //
     *      'or intermediate-output mode of integration, ' //
     *      'respectively.  you have called the code ' //
     *      'with  info(3) = ' // xern1, 5, 1)
         idid = -33
      endif
c
      if (neq .lt. 1) then
         write (xern1, '(i8)') neq
         call xermsg ('slatec', 'derkfs',
     *      'in derkf, the number of equations neq must be a ' //
     *      'positive integer.  you have called the ' //
     *      'code with neq = ' // xern1, 6, 1)
         idid = -33
      endif
c
      nrtolp = 0
      natolp = 0
      do 10 k=1,neq
         if (nrtolp .eq. 0 .and. rtol(k) .lt. 0.d0) then
            write (xern1, '(i8)') k
            write (xern3, '(1pe15.6)') rtol(k)
            call xermsg ('slatec', 'derkfs',
     *         'in derkf, the relative error ' //
     *         'tolerances rtol must be non-negative.  you have ' //
     *         'called the code with  rtol(' // xern1 // ') = ' //
     *         xern3 // '.  in the case of vector error tolerances, ' //
     *         'no further checking of rtol components is done.', 7, 1)
            idid = -33
            nrtolp = 1
         endif
c
         if (natolp .eq. 0 .and. atol(k) .lt. 0.d0) then
            write (xern1, '(i8)') k
            write (xern3, '(1pe15.6)') atol(k)
            call xermsg ('slatec', 'derkfs',
     *         'in derkf, the absolute error ' //
     *         'tolerances atol must be non-negative.  you have ' //
     *         'called the code with  atol(' // xern1 // ') = ' //
     *         xern3 // '.  in the case of vector error tolerances, ' //
     *         'no further checking of atol components is done.', 8, 1)
            idid = -33
            natolp = 1
         endif
c
         if (info(2) .eq. 0) go to 20
         if (natolp.gt.0 .and. nrtolp.gt.0) go to 20
   10 continue
c
c
c     check some continuation possibilities
c
   20 if (init .ne. 0) then
         if (t .eq. tout) then
            write (xern3, '(1pe15.6)') t
            call xermsg ('slatec', 'derkfs',
     *         'in derkf, you have called the ' //
     *         'code with  t = tout = ' // xern3 // '$$this is not ' //
     *         'allowed on continuation calls.', 9, 1)
            idid=-33
         endif
c
         if (t .ne. told) then
            write (xern3, '(1pe15.6)') told
            write (xern4, '(1pe15.6)') t
            call xermsg ('slatec', 'derkfs',
     *         'in derkf, you have changed the ' //
     *         'value of t from ' // xern3 // ' to ' // xern4 //
     *         '$$this is not allowed on continuation calls.', 10, 1)
            idid=-33
         endif
c
         if (init .ne. 1) then
            if (dtsign*(tout-t) .lt. 0.d0) then
               write (xern3, '(1pe15.6)') tout
               call xermsg ('slatec', 'derkfs',
     *            'in derkf, by calling the code ' //
     *            'with tout = ' // xern3 // ' you are attempting ' //
     *            'to change the direction of integration.$$this is ' //
     *            'not allowed without restarting.', 11, 1)
               idid=-33
            endif
         endif
      endif
c
c     invalid input detected
c
      if (idid .eq. (-33)) then
         if (iquit .ne. (-33)) then
            iquit = -33
            goto 909
         else
            call xermsg ('slatec', 'derkfs',
     *         'in derkf, invalid input was ' //
     *         'detected on successive entries.  it is impossible ' //
     *         'to proceed because you have not corrected the ' //
     *         'problem, so execution is being terminated.', 12, 2)
            return
         endif
      endif
c
c.......................................................................
c
c     rtol = atol = 0. is allowed as valid input and interpreted as
c     asking for the most accurate solution possible. in this case,
c     the relative error tolerance rtol is reset to the smallest value
c     rer which is likely to be reasonable for this method and machine.
c
      do 50 k=1,neq
        if (rtol(k)+atol(k) .gt. 0.) go to 45
        rtol(k)=rer
        idid=-2
   45   if (info(2) .eq. 0) go to 55
   50   continue
c
   55 if (idid .ne. (-2)) go to 60
c
c                       rtol=atol=0 on input, so rtol was changed to a
c                                                small positive value
      tolfac=1.
      go to 909
c
c     branch on status of initialization indicator
c            init=0 means initial derivatives and starting step size
c                   not yet computed
c            init=1 means starting step size not yet computed
c            init=2 means no further initialization required
c
   60 if (init .eq. 0) go to 65
      if (init .eq. 1) go to 70
      go to 80
c
c.......................................................................
c
c     more initialization --
c                         -- evaluate initial derivatives
c
   65 init=1
      a=t
      call f(a,y,yp,rpar,ipar)
      if (t .eq. tout) go to 666
c
c                         -- set sign of integration direction  and
c                         -- estimate starting step size
c
   70 init=2
      dtsign=sign(1.,tout-t)
      u=r1mach(4)
      big=sqrt(r1mach(2))
      ute=u**0.375
      dy=ute*hvnrm(y,neq)
      if (dy .eq. 0.) dy=ute
      ktol=1
      do 75 k=1,neq
        if (info(2) .eq. 1)  ktol=k
        tol=rtol(ktol)*abs(y(k))+atol(ktol)
        if (tol .eq. 0.) tol=dy*rtol(ktol)
   75   f1(k)=tol
c
      call hstart (f,neq,t,tout,y,yp,f1,4,u,big,f2,f3,f4,f5,rpar,ipar,h)
c
c.......................................................................
c
c     set step size for integration in the direction from t to tout
c     and set output point indicator
c
   80 dt=tout-t
      h=sign(h,dt)
      output= .false.
c
c     test to see if derkf is being severely impacted by too many
c     output points
c
      if (abs(h) .ge. 2.*abs(dt)) kop=kop+1
      if (kop .le. mxkop) go to 85
c
c                       unnecessary frequency of output is restricting
c                                                 the step size choice
      idid=-5
      kop=0
      go to 909
c
   85 if (abs(dt) .gt. u26*abs(t)) go to 100
c
c     if too close to output point,extrapolate and return
c
      do 90 k=1,neq
   90   y(k)=y(k)+dt*yp(k)
      a=tout
      call f(a,y,yp,rpar,ipar)
      ksteps=ksteps+1
      go to 666
c
c **********************************************************************
c **********************************************************************
c     step by step integration
c
  100 hfaild= .false.
c
c     to protect against impossible accuracy requests, compute a
c     tolerance factor based on the requested error tolerance and a
c     level of accuracy achievable at limiting precision
c
      tolfac=0.
      ktol=1
      do 125 k=1,neq
        if (info(2) .eq. 1) ktol=k
        et=rtol(ktol)*abs(y(k))+atol(ktol)
        if (et .gt. 0.) go to 120
        tolfac=max(tolfac,rer/rtol(ktol))
        go to 125
  120   tolfac=max(tolfac,abs(y(k))*(rer/et))
  125   continue
      if (tolfac .le. 1.) go to 150
c
c                       requested error unattainable due to limited
c                                               precision available
      tolfac=2.*tolfac
      idid=-2
      go to 909
c
c     set smallest allowable step size
c
  150 hmin=u26*abs(t)
c
c     adjust step size if necessary to hit the output point --
c     look ahead two steps to avoid drastic changes in the step size and
c     thus lessen the impact of output points on the code.
c     stretch the step size by, at most, an amount equal to the
c     safety factor of 9/10.
c
      dt=tout-t
      if (abs(dt) .ge. 2.*abs(h)) go to 200
      if (abs(dt) .gt. abs(h)/0.9) go to 175
c
c     the next step, if successful, will complete the integration to
c     the output point
c
      output= .true.
      h=dt
      go to 200
c
  175 h=0.5*dt
c
c
c **********************************************************************
c     core integrator for taking a single step
c **********************************************************************
c     to avoid problems with zero crossings, relative error is measured
c     using the average of the magnitudes of the solution at the
c     beginning and end of a step.
c     the error estimate formula has been grouped to control loss of
c     significance.
c     local error estimates for a first order method using the same
c     step size as the fehlberg method are calculated as part of the
c     test for stiffness.
c     to distinguish the various arguments, h is not permitted
c     to become smaller than 26 units of roundoff in t.
c     practical limits on the change in the step size are enforced to
c     smooth the step size selection process and to avoid excessive
c     chattering on problems having discontinuities.
c     to prevent unnecessary failures, the code uses 9/10 the step size
c     it estimates will succeed.
c     after a step failure, the step size is not allowed to increase for
c     the next attempted step. this makes the code more efficient on
c     problems having discontinuities and more effective in general
c     since local extrapolation is being used and extra caution seems
c     warranted.
c.......................................................................
c
c     monitor number of steps attempted
c
  200 if (ksteps .le. mxstep) go to 222
c
c                       a significant amount of work has been expended
      idid=-1
      ksteps=0
      if (.not. stiff) go to 909
c
c                       problem appears to be stiff
      idid=-4
      stiff= .false.
      nonstf= .false.
      ntstep=0
      nstifs=0
      go to 909
c
c     advance an approximate solution over one step of length h
c
  222 call defehl(f,neq,t,y,h,yp,f1,f2,f3,f4,f5,ys,rpar,ipar)
      ksteps=ksteps+1
c
c.......................................................................
c
c     compute and test allowable tolerances versus local error
c     estimates.  note that relative error is measured with respect to
c     the average of the magnitudes of the solution at the beginning
c     and end of the step.
c     local error estimates for a special first order method are
c     calculated only when the stiffness detection is turned on.
c
      eeoet=0.
      estiff=0.
      ktol=1
      do 350 k=1,neq
        yavg=0.5*(abs(y(k))+abs(ys(k)))
        if (info(2) .eq. 1) ktol=k
        et=rtol(ktol)*yavg+atol(ktol)
        if (et .gt. 0.) go to 325
c
c                       pure relative error inappropriate when solution
c                                                              vanishes
        idid=-3
        go to 909
c
  325   ee=abs((-2090.*yp(k)+(21970.*f3(k)-15048.*f4(k)))+
     1                        (22528.*f2(k)-27360.*f5(k)))
        if (stiff .or. nonstf) go to 350
        es=abs(h*(0.055455*yp(k)-0.035493*f1(k)-0.036571*f2(k)+
     1            0.023107*f3(k)-0.009515*f4(k)+0.003017*f5(k)))
        estiff=max(estiff,es/et)
  350   eeoet=max(eeoet,ee/et)
c
      esttol=abs(h)*eeoet/752400.
c
      if (esttol .le. 1.) go to 500
c
c.......................................................................
c
c     unsuccessful step
c
      if (abs(h) .gt. hmin) go to 400
c
c                       requested error unattainable at smallest
c                                            allowable step size
      tolfac=1.69*esttol
      idid=-2
      go to 909
c
c                       reduce the step size , try again
c                       the decrease is limited to a factor of 1/10
c
  400 hfaild= .true.
      output= .false.
      s=0.1
      if (esttol .lt. 59049.) s=0.9/esttol**0.2
      h=sign(max(s*abs(h),hmin),h)
      go to 200
c
c.......................................................................
c
c     successful step
c                       store solution at t+h
c                       and evaluate derivatives there
c
  500 t=t+h
      do 525 k=1,neq
  525   y(k)=ys(k)
      a=t
      call f(a,y,yp,rpar,ipar)
c
c                       choose next step size
c                       the increase is limited to a factor of 5
c                       if step failure has just occurred, next
c                          step size is not allowed to increase
c
      s=5.
      if (esttol .gt. 1.889568e-4) s=0.9/esttol**0.2
      if (hfaild) s=min(s,1.)
      h=sign(max(s*abs(h),hmin),h)
c
c.......................................................................
c
c     check for stiffness (if not already detected)
c
c     in a sequence of 50 successful steps by the fehlberg method, 25
c     successful steps by the first order method indicates stiffness
c     and turns the test off. if 26 failures by the first order method
c     occur, the test is turned off until this sequence of 50 steps
c     by the fehlberg method is completed.
c
      if (stiff) go to 600
      ntstep=mod(ntstep+1,50)
      if (ntstep .eq. 1) nonstf= .false.
      if (nonstf) go to 600
      if (estiff .gt. 1.) go to 550
c
c                       successful step with first order method
      nstifs=nstifs+1
c                       turn test off after 25 indications of stiffness
      if (nstifs .eq. 25) stiff= .true.
      go to 600
c
c                       unsuccessful step with first order method
  550 if (ntstep-nstifs .le. 25) go to 600
c                       turn stiffness detection off for this block of
c                                                          fifty steps
      nonstf= .true.
c                       reset stiff step counter
      nstifs=0
c
c **********************************************************************
c     end of core integrator
c **********************************************************************
c
c
c     should we take another step
c
  600 if (output) go to 666
      if (info(3) .eq. 0) go to 100
c
c **********************************************************************
c **********************************************************************
c
c     integration successfully completed
c
c                 one-step mode
      idid=1
      told=t
      return
c
c                 interval mode
  666 idid=2
      t=tout
      told=t
      return
c
c     integration task interrupted
c
  909 info(1)=-1
      told=t
      if (idid .ne. (-2)) return
c
c                       the error tolerances are increased to values
c                               which are appropriate for continuing
      rtol(1)=tolfac*rtol(1)
      atol(1)=tolfac*atol(1)
      if (info(2) .eq. 0) return
      do 939 k=2,neq
        rtol(k)=tolfac*rtol(k)
  939   atol(k)=tolfac*atol(k)
      return
      end
