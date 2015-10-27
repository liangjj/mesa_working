*deck ddes
      subroutine ddes (df, neq, t, y, tout, info, rtol, atol, idid,
     +   ypout, yp, yy, wt, p, phi, alpha, beta, psi, v, w, sig, g, gi,
     +   h, eps, x, xold, hold, told, delsgn, tstop, twou, fouru, start,
     +   phase1, nornd, stiff, intout, ns, kord, kold, init, ksteps,
     +   kle4, iquit, kprev, ivc, iv, kgi, rpar, ipar)
c***begin prologue  ddes
c***subsidiary
c***purpose  subsidiary to ddeabm
c***library   slatec
c***type      double precision (des-s, ddes-d)
c***author  watts, h. a., (snla)
c***description
c
c   ddeabm merely allocates storage for ddes to relieve the user of the
c   inconvenience of a long call list.  consequently  ddes  is used as
c   described in the comments for  ddeabm .
c
c***see also  ddeabm
c***routines called  d1mach, dintp, dsteps, xermsg
c***revision history  (yymmdd)
c   820301  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   900510  convert xerrwv calls to xermsg calls, cvt gotos to
c           if-then-else.  (rwc)
c   910722  updated author section.  (als)
c***end prologue  ddes
c
      integer idid, info, init, ipar, iquit, iv, ivc, k, kgi, kle4,
     1      kold, kord, kprev, ksteps, l, ltol, maxnum, natolp, neq,
     2      nrtolp, ns
      double precision a, absdel, alpha, atol, beta, d1mach,
     1      del, delsgn, dt, eps, fouru, g, gi, h,
     2      ha, hold, p, phi, psi, rpar, rtol, sig, t, told, tout,
     3      tstop, twou, u, v, w, wt, x, xold, y, yp, ypout, yy
      logical stiff,crash,start,phase1,nornd,intout
c
      dimension y(*),yy(*),wt(*),phi(neq,16),p(*),yp(*),
     1  ypout(*),psi(12),alpha(12),beta(12),sig(13),v(12),w(12),g(13),
     2  gi(11),iv(10),info(15),rtol(*),atol(*),rpar(*),ipar(*)
      character*8 xern1
      character*16 xern3, xern4
c
      external df
c
c.......................................................................
c
c  the expense of solving the problem is monitored by counting the
c  number of  steps attempted. when this exceeds  maxnum, the counter
c  is reset to zero and the user is informed about possible excessive
c  work.
c
      save maxnum
      data maxnum/500/
c
c.......................................................................
c
c***first executable statement  ddes
      if (info(1) .eq. 0) then
c
c on the first call , perform initialization --
c        define the machine unit roundoff quantity  u  by calling the
c        function routine  d1mach. the user must make sure that the
c        values set in d1mach are relevant to the computer being used.
c
         u=d1mach(4)
c                       -- set associated machine dependent parameters
         twou=2.d0*u
         fouru=4.d0*u
c                       -- set termination flag
         iquit=0
c                       -- set initialization indicator
         init=0
c                       -- set counter for attempted steps
         ksteps=0
c                       -- set indicator for intermediate-output
         intout= .false.
c                       -- set indicator for stiffness detection
         stiff= .false.
c                       -- set step counter for stiffness detection
         kle4=0
c                       -- set indicators for steps code
         start= .true.
         phase1= .true.
         nornd= .true.
c                       -- reset info(1) for subsequent calls
         info(1)=1
      endif
c
c.......................................................................
c
c      check validity of input parameters on each entry
c
      if (info(1) .ne. 0  .and.  info(1) .ne. 1) then
         write (xern1, '(i8)') info(1)
         call xermsg ('slatec', 'ddes', 'in ddeabm, info(1) must be ' //
     *      'set to 0 for the start of a new problem, and must be ' //
     *      'set to 1 following an interrupted task.  you are ' //
     *      'attempting to continue the integration illegally by ' //
     *      'calling the code with info(1) = ' // xern1, 3, 1)
         idid=-33
      endif
c
      if (info(2) .ne. 0  .and.  info(2) .ne. 1) then
         write (xern1, '(i8)') info(2)
         call xermsg ('slatec', 'ddes', 'in ddeabm, info(2) must be ' //
     *      '0 or 1 indicating scalar and vector error tolerances, ' //
     *      'respectively.  you have called the code with info(2) = ' //
     *      xern1, 4, 1)
         idid=-33
      endif
c
      if (info(3) .ne. 0  .and.  info(3) .ne. 1) then
         write (xern1, '(i8)') info(3)
         call xermsg ('slatec', 'ddes', 'in ddeabm, info(3) must be ' //
     *      '0 or 1 indicating the interval or intermediate-output ' //
     *      'mode of integration, respectively.  you have called ' //
     *      'the code with  info(3) = ' // xern1, 5, 1)
         idid=-33
      endif
c
      if (info(4) .ne. 0  .and.  info(4) .ne. 1) then
         write (xern1, '(i8)') info(4)
         call xermsg ('slatec', 'ddes', 'in ddeabm, info(4) must be ' //
     *      '0 or 1 indicating whether or not the integration ' //
     *      'interval is to be restricted by a point tstop.  you ' //
     *      'have called the code with info(4) = ' // xern1, 14, 1)
         idid=-33
      endif
c
      if (neq .lt. 1) then
         write (xern1, '(i8)') neq
         call xermsg ('slatec', 'ddes', 'in ddeabm,  the number of ' //
     *      'equations neq must be a positive integer.  you have ' //
     *      'called the code with  neq = ' // xern1, 6, 1)
         idid=-33
      endif
c
      nrtolp = 0
      natolp = 0
      do 90 k=1,neq
         if (nrtolp .eq. 0 .and. rtol(k) .lt. 0.d0) then
            write (xern1, '(i8)') k
            write (xern3, '(1pe15.6)') rtol(k)
            call xermsg ('slatec', 'ddes', 'in ddeabm, the relative ' //
     *         'error tolerances rtol must be non-negative.  you ' //
     *         'have called the code with  rtol(' // xern1 // ') = ' //
     *         xern3 // '.  in the case of vector error tolerances, ' //
     *         'no further checking of rtol components is done.', 7, 1)
            idid = -33
            nrtolp = 1
         endif
c
         if (natolp .eq. 0 .and. atol(k) .lt. 0.d0) then
            write (xern1, '(i8)') k
            write (xern3, '(1pe15.6)') atol(k)
            call xermsg ('slatec', 'ddes', 'in ddeabm, the absolute ' //
     *         'error tolerances atol must be non-negative.  you ' //
     *         'have called the code with  atol(' // xern1 // ') = ' //
     *         xern3 // '.  in the case of vector error tolerances, ' //
     *         'no further checking of atol components is done.', 8, 1)
            idid = -33
            natolp = 1
         endif
c
         if (info(2) .eq. 0) go to 100
         if (natolp.gt.0 .and. nrtolp.gt.0) go to 100
   90 continue
c
  100 if (info(4) .eq. 1) then
         if (sign(1.d0,tout-t) .ne. sign(1.d0,tstop-t)
     1      .or. abs(tout-t) .gt. abs(tstop-t)) then
            write (xern3, '(1pe15.6)') tout
            write (xern4, '(1pe15.6)') tstop
            call xermsg ('slatec', 'ddes', 'in ddeabm, you have ' //
     *         'called the code with  tout = ' // xern3 // ' but ' //
     *         'you have also told the code (info(4) = 1) not to ' //
     *         'integrate past the point tstop = ' // xern4 //
     *         ' these instructions conflict.', 14, 1)
            idid=-33
         endif
      endif
c
c     check some continuation possibilities
c
      if (init .ne. 0) then
         if (t .eq. tout) then
            write (xern3, '(1pe15.6)') t
            call xermsg ('slatec', 'ddes', 'in ddeabm, you have ' //
     *         'called the code with  t = tout = ' // xern3 //
     *         '$$this is not allowed on continuation calls.', 9, 1)
            idid=-33
         endif
c
         if (t .ne. told) then
            write (xern3, '(1pe15.6)') told
            write (xern4, '(1pe15.6)') t
            call xermsg ('slatec', 'ddes', 'in ddeabm, you have ' //
     *         'changed the value of t from ' // xern3 // ' to ' //
     *         xern4 //'  this is not allowed on continuation calls.',
     *         10, 1)
            idid=-33
         endif
c
         if (init .ne. 1) then
            if (delsgn*(tout-t) .lt. 0.d0) then
               write (xern3, '(1pe15.6)') tout
               call xermsg ('slatec', 'ddes', 'in ddeabm, by ' //
     *            'calling the code with tout = ' // xern3 //
     *            ' you are attempting to change the direction of ' //
     *            'integration.$$this is not allowed without ' //
     *            'restarting.', 11, 1)
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
            info(1) = -1
         else
            call xermsg ('slatec', 'ddes', 'in ddeabm, invalid ' //
     *         'input was detected on successive entries.  it is ' //
     *         'impossible to proceed because you have not ' //
     *         'corrected the problem, so execution is being ' //
     *         'terminated.', 12, 2)
         endif
         return
      endif
c
c.......................................................................
c
c     rtol = atol = 0. is allowed as valid input and interpreted as
c     asking for the most accurate solution possible. in this case,
c     the relative error tolerance rtol is reset to the smallest value
c     fouru which is likely to be reasonable for this method and machine
c
      do 180 k=1,neq
        if (rtol(k)+atol(k) .gt. 0.d0) go to 170
        rtol(k)=fouru
        idid=-2
  170   if (info(2) .eq. 0) go to 190
  180   continue
c
  190 if (idid .ne. (-2)) go to 200
c                       rtol=atol=0 on input, so rtol is changed to a
c                                                small positive value
      info(1)=-1
      return
c
c     branch on status of initialization indicator
c            init=0 means initial derivatives and nominal step size
c                   and direction not yet set
c            init=1 means nominal step size and direction not yet set
c            init=2 means no further initialization required
c
  200 if (init .eq. 0) go to 210
      if (init .eq. 1) go to 220
      go to 240
c
c.......................................................................
c
c     more initialization --
c                         -- evaluate initial derivatives
c
  210 init=1
      a=t
      call df(a,y,yp,rpar,ipar)
      if (t .ne. tout) go to 220
      idid=2
      do 215 l = 1,neq
  215    ypout(l) = yp(l)
      told=t
      return
c
c                         -- set independent and dependent variables
c                                              x and yy(*) for steps
c                         -- set sign of integration direction
c                         -- initialize the step size
c
  220 init = 2
      x = t
      do 230 l = 1,neq
  230   yy(l) = y(l)
      delsgn = sign(1.0d0,tout-t)
      h = sign(max(fouru*abs(x),abs(tout-x)),tout-x)
c
c.......................................................................
c
c   on each call set information which determines the allowed interval
c   of integration before returning with an answer at tout
c
  240 del = tout - t
      absdel = abs(del)
c
c.......................................................................
c
c   if already past output point, interpolate and return
c
  250 if(abs(x-t) .lt. absdel) go to 260
      call dintp(x,yy,tout,y,ypout,neq,kold,phi,ivc,iv,kgi,gi,
     1                                        alpha,g,w,xold,p)
      idid = 3
      if (x .ne. tout) go to 255
      idid = 2
      intout = .false.
  255 t = tout
      told = t
      return
c
c   if cannot go past tstop and sufficiently close,
c   extrapolate and return
c
  260 if (info(4) .ne. 1) go to 280
      if (abs(tstop-x) .ge. fouru*abs(x)) go to 280
      dt = tout - x
      do 270 l = 1,neq
  270   y(l) = yy(l) + dt*yp(l)
      call df(tout,y,ypout,rpar,ipar)
      idid = 3
      t = tout
      told = t
      return
c
  280 if (info(3) .eq. 0  .or.  .not.intout) go to 300
c
c   intermediate-output mode
c
      idid = 1
      do 290 l = 1,neq
        y(l)=yy(l)
  290   ypout(l) = yp(l)
      t = x
      told = t
      intout = .false.
      return
c
c.......................................................................
c
c     monitor number of steps attempted
c
  300 if (ksteps .le. maxnum) go to 330
c
c                       a significant amount of work has been expended
      idid=-1
      ksteps=0
      if (.not. stiff) go to 310
c
c                       problem appears to be stiff
      idid=-4
      stiff= .false.
      kle4=0
c
  310 do 320 l = 1,neq
        y(l) = yy(l)
  320   ypout(l) = yp(l)
      t = x
      told = t
      info(1) = -1
      intout = .false.
      return
c
c.......................................................................
c
c   limit step size, set weight vector and take a step
c
  330 ha = abs(h)
      if (info(4) .ne. 1) go to 340
      ha = min(ha,abs(tstop-x))
  340 h = sign(ha,h)
      eps = 1.0d0
      ltol = 1
      do 350 l = 1,neq
        if (info(2) .eq. 1) ltol = l
        wt(l) = rtol(ltol)*abs(yy(l)) + atol(ltol)
        if (wt(l) .le. 0.0d0) go to 360
  350   continue
      go to 380
c
c                       relative error criterion inappropriate
  360 idid = -3
      do 370 l = 1,neq
        y(l) = yy(l)
  370   ypout(l) = yp(l)
      t = x
      told = t
      info(1) = -1
      intout = .false.
      return
c
  380 call dsteps(df,neq,yy,x,h,eps,wt,start,hold,kord,kold,crash,phi,p,
     1           yp,psi,alpha,beta,sig,v,w,g,phase1,ns,nornd,ksteps,
     2           twou,fouru,xold,kprev,ivc,iv,kgi,gi,rpar,ipar)
c
c.......................................................................
c
      if(.not.crash) go to 420
c
c                       tolerances too small
      idid = -2
      rtol(1) = eps*rtol(1)
      atol(1) = eps*atol(1)
      if (info(2) .eq. 0) go to 400
      do 390 l = 2,neq
        rtol(l) = eps*rtol(l)
  390   atol(l) = eps*atol(l)
  400 do 410 l = 1,neq
        y(l) = yy(l)
  410   ypout(l) = yp(l)
      t = x
      told = t
      info(1) = -1
      intout = .false.
      return
c
c   (stiffness test) count number of consecutive steps taken with the
c   order of the method being less or equal to four
c
  420 kle4 = kle4 + 1
      if(kold .gt. 4) kle4 = 0
      if(kle4 .ge. 50) stiff = .true.
      intout = .true.
      go to 250
      end
