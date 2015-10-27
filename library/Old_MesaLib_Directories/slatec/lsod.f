*deck lsod
      subroutine lsod (f, neq, t, y, tout, rtol, atol, idid, ypout, yh,
     +   yh1, ewt, savf, acor, wm, iwm, jac, intout, tstop, tolfac,
     +   delsgn, rpar, ipar)
c***begin prologue  lsod
c***subsidiary
c***purpose  subsidiary to debdf
c***library   slatec
c***type      single precision (lsod-s, dlsod-d)
c***author  (unknown)
c***description
c
c   debdf  merely allocates storage for  lsod  to relieve the user of
c   the inconvenience of a long call list.  consequently  lsod  is used
c   as described in the comments for  debdf .
c
c***see also  debdf
c***routines called  hstart, intyd, r1mach, stod, vnwrms, xermsg
c***common blocks    debdf1
c***revision history  (yymmdd)
c   800901  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   900510  convert xerrwv calls to xermsg calls.  (rwc)
c***end prologue  lsod
c
c
      logical intout
c
      dimension y(*),ypout(*),yh(neq,6),yh1(*),ewt(*),savf(*),
     1          acor(*),wm(*),iwm(*),rtol(*),atol(*),rpar(*),ipar(*)
      character*8 xern1
      character*16 xern3, xern4
c
      common /debdf1/ told, rowns(210),
     1   el0, h, hmin, hmxi, hu, x, u,
     2   iquit, init, lyh, lewt, lacor, lsavf, lwm, ksteps,
     3   ibegin, itol, iinteg, itstop, ijac, iband, iowns(6),
     4   ier, jstart, kflag, ldum, meth, miter, maxord, n, nq, nst,
     5   nfe, nje, nqu
c
      external f, jac
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
c***first executable statement  lsod
      if (ibegin .eq. 0) then
c
c        on the first call , perform initialization --
c        define the machine unit roundoff quantity  u  by calling the
c        function routine r1mach. the user must make sure that the
c        values set in r1mach are relevant to the computer being used.
c
         u = r1mach(4)
c                          -- set associated machine dependent parameter
         wm(1) = sqrt(u)
c                          -- set termination flag
         iquit = 0
c                          -- set initialization indicator
         init = 0
c                          -- set counter for attempted steps
         ksteps = 0
c                          -- set indicator for intermediate-output
         intout = .false.
c                          -- set start indicator for stod code
         jstart = 0
c                          -- set bdf method indicator
         meth = 2
c                          -- set maximum order for bdf method
         maxord = 5
c                          -- set iteration matrix indicator
c
         if (ijac .eq. 0 .and. iband .eq. 0) miter = 2
         if (ijac .eq. 1 .and. iband .eq. 0) miter = 1
         if (ijac .eq. 0 .and. iband .eq. 1) miter = 5
         if (ijac .eq. 1 .and. iband .eq. 1) miter = 4
c
c                          -- set other necessary items in common block
         n = neq
         nst = 0
         nje = 0
         hmxi = 0.
         nq = 1
         h = 1.
c                          -- reset ibegin for subsequent calls
         ibegin=1
      endif
c
c.......................................................................
c
c      check validity of input parameters on each entry
c
      if (neq .lt. 1) then
         write (xern1, '(i8)') neq
         call xermsg ('slatec', 'lsod',
     *      'in debdf, the number of equations must be a positive ' //
     *      'integer.$$you have called the code with neq = ' // xern1,
     *      6, 1)
         idid=-33
      endif
c
      nrtolp = 0
      natolp = 0
      do 60 k = 1,neq
         if (nrtolp .le. 0) then
            if (rtol(k) .lt. 0.) then
               write (xern1, '(i8)') k
               write (xern3, '(1pe15.6)') rtol(k)
               call xermsg ('slatec', 'lsod',
     *            'in debdf, the relative error tolerances must ' //
     *            'be non-negative.$$you have called the code with ' //
     *            'rtol(' // xern1 // ') = ' // xern3 // '$$in the ' //
     *            'case of vector error tolerances, no further ' //
     *            'checking of rtol components is done.', 7, 1)
               idid = -33
               if (natolp .gt. 0) go to 70
               nrtolp = 1
            elseif (natolp .gt. 0) then
               go to 50
            endif
         endif
c
         if (atol(k) .lt. 0.) then
            write (xern1, '(i8)') k
            write (xern3, '(1pe15.6)') atol(k)
            call xermsg ('slatec', 'lsod',
     *         'in debdf, the absolute error ' //
     *         'tolerances must be non-negative.$$you have called ' //
     *         'the code with atol(' // xern1 // ') = ' // xern3 //
     *         '$$in the case of vector error tolerances, no further '
     *         // 'checking of atol components is done.', 8, 1)
            idid=-33
            if (nrtolp .gt. 0) go to 70
            natolp=1
         endif
   50    if (itol .eq. 0) go to 70
   60 continue
c
   70 if (itstop .eq. 1) then
         if (sign(1.,tout-t) .ne. sign(1.,tstop-t) .or.
     1      abs(tout-t) .gt. abs(tstop-t)) then
            write (xern3, '(1pe15.6)') tout
            write (xern4, '(1pe15.6)') tstop
            call xermsg ('slatec', 'lsod',
     *         'in debdf, you have called the ' //
     *         'code with tout = ' // xern3 // '$$but you have ' //
     *         'also told the code not to integrate past the point ' //
     *         'tstop = ' // xern4 // ' by setting info(4) = 1.  ' //
     *         'these instructions conflict.', 14, 1)
            idid=-33
         endif
      endif
c
c        check some continuation possibilities
c
      if (init .ne. 0) then
         if (t .eq. tout) then
            write (xern3, '(1pe15.6)') t
            call xermsg ('slatec', 'lsod',
     *         'in debdf, you have called the code with t = tout = ' //
     *         xern3 // '  this is not allowed on continuation calls.',
     *         9, 1)
            idid=-33
         endif
c
         if (t .ne. told) then
            write (xern3, '(1pe15.6)') told
            write (xern4, '(1pe15.6)') t
            call xermsg ('slatec', 'lsod',
     *         'in debdf, you have changed the value of t from ' //
     *         xern3 // ' to ' // xern4 //
     *         '  this is not allowed on continuation calls.', 10, 1)
            idid=-33
         endif
c
         if (init .ne. 1) then
            if (delsgn*(tout-t) .lt. 0.) then
               write (xern3, '(1pe15.6)') tout
               call xermsg ('slatec', 'lsod',
     *            'in debdf, by calling the code with tout = ' //
     *            xern3 // ' you are attempting to change the ' //
     *            'direction of integration.$$' //
     *            'this is not allowed without restarting.', 11, 1)
               idid=-33
            endif
         endif
      endif
c
      if (idid .eq. (-33)) then
         if (iquit .ne. (-33)) then
c                       invalid input detected
            iquit=-33
            ibegin=-1
         else
            call xermsg ('slatec', 'lsod',
     *         'in debdf, invalid input was ' //
     *         'detected on successive entries.  it is impossible ' //
     *         'to proceed because you have not corrected the ' //
     *         'problem, so execution is being terminated.', 12, 2)
         endif
         return
      endif
c
c.......................................................................
c
c     rtol = atol = 0. is allowed as valid input and interpreted as
c     asking for the most accurate solution possible. in this case,
c     the relative error tolerance rtol is reset to the smallest value
c     100*u which is likely to be reasonable for this method and machine
c
      do 170 k=1,neq
        if (rtol(k)+atol(k) .gt. 0.) go to 160
        rtol(k)=100.*u
        idid=-2
  160   if (itol .eq. 0) go to 180
  170   continue
c
  180 if (idid .ne. (-2)) go to 190
c                       rtol=atol=0 on input, so rtol is changed to a
c                                                small positive value
      ibegin=-1
      return
c
c     branch on status of initialization indicator
c            init=0 means initial derivatives and nominal step size
c                   and direction not yet set
c            init=1 means nominal step size and direction not yet set
c            init=2 means no further initialization required
c
  190 if (init .eq. 0) go to 200
      if (init .eq. 1) go to 220
      go to 240
c
c.......................................................................
c
c     more initialization --
c                         -- evaluate initial derivatives
c
  200 init=1
      call f(t,y,yh(1,2),rpar,ipar)
      nfe=1
      if (t .ne. tout) go to 220
      idid=2
      do 210 l = 1,neq
  210    ypout(l) = yh(l,2)
      told=t
      return
c
c                         -- compute initial step size
c                         -- save sign of integration direction
c                         -- set independent and dependent variables
c                                              x and yh(*) for stod
c
  220 ltol = 1
      do 225 l=1,neq
        if (itol .eq. 1) ltol = l
        tol = rtol(ltol)*abs(y(l)) + atol(ltol)
        if (tol .eq. 0.) go to 380
  225   ewt(l) = tol
c
      big = sqrt(r1mach(2))
      call hstart (f,neq,t,tout,y,yh(1,2),ewt,1,u,big,
     1             yh(1,3),yh(1,4),yh(1,5),yh(1,6),rpar,ipar,h)
c
      delsgn = sign(1.0,tout-t)
      x = t
      do 230 l = 1,neq
        yh(l,1) = y(l)
  230   yh(l,2) = h*yh(l,2)
      init = 2
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
  250 if (abs(x-t) .lt. absdel) go to 270
      call intyd(tout,0,yh,neq,y,intflg)
      call intyd(tout,1,yh,neq,ypout,intflg)
      idid = 3
      if (x .ne. tout) go to 260
      idid = 2
      intout = .false.
  260 t = tout
      told = t
      return
c
c   if cannot go past tstop and sufficiently close,
c   extrapolate and return
c
  270 if (itstop .ne. 1) go to 290
      if (abs(tstop-x) .ge. 100.*u*abs(x)) go to 290
      dt = tout - x
      do 280 l = 1,neq
  280   y(l) = yh(l,1) + (dt/h)*yh(l,2)
      call f(tout,y,ypout,rpar,ipar)
      nfe = nfe + 1
      idid = 3
      t = tout
      told = t
      return
c
  290 if (iinteg .eq. 0  .or.  .not.intout) go to 300
c
c   intermediate-output mode
c
      idid = 1
      go to 500
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
      ibegin = -1
      go to 500
c
c.......................................................................
c
c   limit step size and set weight vector
c
  330 hmin = 100.*u*abs(x)
      ha = max(abs(h),hmin)
      if (itstop .ne. 1) go to 340
      ha = min(ha,abs(tstop-x))
  340 h = sign(ha,h)
      ltol = 1
      do 350 l = 1,neq
        if (itol .eq. 1) ltol = l
        ewt(l) = rtol(ltol)*abs(yh(l,1)) + atol(ltol)
        if (ewt(l) .le. 0.0) go to 380
  350   continue
      tolfac = u*vnwrms(neq,yh,ewt)
      if (tolfac .le. 1.) go to 400
c
c                       tolerances too small
      idid = -2
      tolfac = 2.*tolfac
      rtol(1) = tolfac*rtol(1)
      atol(1) = tolfac*atol(1)
      if (itol .eq. 0) go to 370
      do 360 l = 2,neq
        rtol(l) = tolfac*rtol(l)
  360   atol(l) = tolfac*atol(l)
  370 ibegin = -1
      go to 500
c
c                       relative error criterion inappropriate
  380 idid = -3
      ibegin = -1
      go to 500
c
c.......................................................................
c
c     take a step
c
  400 call stod(neq,y,yh,neq,yh1,ewt,savf,acor,wm,iwm,f,jac,rpar,ipar)
c
      jstart = -2
      intout = .true.
      if (kflag .eq. 0) go to 250
c
c.......................................................................
c
      if (kflag .eq. -1) go to 450
c
c                       repeated corrector convergence failures
      idid = -6
      ibegin = -1
      go to 500
c
c                       repeated error test failures
  450 idid = -7
      ibegin = -1
c
c.......................................................................
c
c                       store values before returning to debdf
  500 do 555 l = 1,neq
        y(l) = yh(l,1)
  555   ypout(l) = yh(l,2)/h
      t = x
      told = t
      intout = .false.
      return
      end
