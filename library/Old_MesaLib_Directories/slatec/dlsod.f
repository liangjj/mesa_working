*deck dlsod
      subroutine dlsod (df, neq, t, y, tout, rtol, atol, idid, ypout,
     +   yh, yh1, ewt, savf, acor, wm, iwm, djac, intout, tstop, tolfac,
     +   delsgn, rpar, ipar)
c***begin prologue  dlsod
c***subsidiary
c***purpose  subsidiary to ddebdf
c***library   slatec
c***type      double precision (lsod-s, dlsod-d)
c***author  (unknown)
c***description
c
c   ddebdf  merely allocates storage for  dlsod  to relieve the user of
c   the inconvenience of a long call list.  consequently  dlsod  is used
c   as described in the comments for  ddebdf .
c
c***see also  ddebdf
c***routines called  d1mach, dhstrt, dintyd, dstod, dvnrms, xermsg
c***common blocks    ddebd1
c***revision history  (yymmdd)
c   820301  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   900510  convert xerrwv calls to xermsg calls.  (rwc)
c***end prologue  dlsod
c
      integer iband, ibegin, idid, ier, iinteg, ijac, init, intflg,
     1      iowns, ipar, iquit, itol, itstop, iwm, jstart, k, kflag,
     2      ksteps, l, lacor, ldum, lewt, lsavf, ltol, lwm, lyh, maxnum,
     3      maxord, meth, miter, n, natolp, neq, nfe, nje, nq, nqu,
     4      nrtolp, nst
      double precision absdel, acor, atol, big, d1mach, del,
     1      delsgn, dt, dvnrms, el0, ewt,
     2      h, ha, hmin, hmxi, hu, rowns, rpar, rtol, savf, t, tol,
     3      told, tolfac, tout, tstop, u, wm, x, y, yh, yh1, ypout
      logical intout
      character*8 xern1
      character*16 xern3, xern4
c
      dimension y(*),ypout(*),yh(neq,6),yh1(*),ewt(*),savf(*),
     1          acor(*),wm(*),iwm(*),rtol(*),atol(*),rpar(*),ipar(*)
c
c
      common /ddebd1/ told,rowns(210),el0,h,hmin,hmxi,hu,x,u,iquit,init,
     1                lyh,lewt,lacor,lsavf,lwm,ksteps,ibegin,itol,
     2                iinteg,itstop,ijac,iband,iowns(6),ier,jstart,
     3                kflag,ldum,meth,miter,maxord,n,nq,nst,nfe,nje,nqu
c
      external df, djac
c
c     ..................................................................
c
c       the expense of solving the problem is monitored by counting the
c       number of  steps attempted. when this exceeds  maxnum, the
c       counter is reset to zero and the user is informed about possible
c       excessive work.
      save maxnum
c
      data maxnum /500/
c
c     ..................................................................
c
c***first executable statement  dlsod
      if (ibegin .eq. 0) then
c
c        on the first call , perform initialization --
c        define the machine unit roundoff quantity  u  by calling the
c        function routine d1mach. the user must make sure that the
c        values set in d1mach are relevant to the computer being used.
c
         u = d1mach(4)
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
c                          -- set start indicator for dstod code
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
         hmxi = 0.0d0
         nq = 1
         h = 1.0d0
c                          -- reset ibegin for subsequent calls
         ibegin = 1
      endif
c
c     ..................................................................
c
c      check validity of input parameters on each entry
c
      if (neq .lt. 1) then
         write (xern1, '(i8)') neq
         call xermsg ('slatec', 'dlsod',
     *      'in ddebdf, the number of equations must be a ' //
     *      'positive integer.$$you have called the code with neq = ' //
     *      xern1, 6, 1)
         idid=-33
      endif
c
      nrtolp = 0
      natolp = 0
      do 60 k = 1, neq
         if (nrtolp .le. 0) then
            if (rtol(k) .lt. 0.) then
               write (xern1, '(i8)') k
               write (xern3, '(1pe15.6)') rtol(k)
               call xermsg ('slatec', 'dlsod',
     *            'in ddebdf, the relative error tolerances must ' //
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
            call xermsg ('slatec', 'dlsod',
     *         'in ddebdf, the absolute error ' //
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
         if (sign(1.0d0,tout-t) .ne. sign(1.0d0,tstop-t) .or.
     1      abs(tout-t) .gt. abs(tstop-t)) then
            write (xern3, '(1pe15.6)') tout
            write (xern4, '(1pe15.6)') tstop
            call xermsg ('slatec', 'dlsod',
     *         'in ddebdf, you have called the ' //
     *         'code with tout = ' // xern3 // '$$but you have ' //
     *         'also told the code not to integrate past the point ' //
     *         'tstop = ' // xern4 // ' by setting info(4) = 1.$$' //
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
            call xermsg ('slatec', 'dlsod',
     *         'in ddebdf, you have called the code with t = tout = ' //
     *         xern3 // '$$this is not allowed on continuation calls.',
     *         9, 1)
            idid=-33
         endif
c
         if (t .ne. told) then
            write (xern3, '(1pe15.6)') told
            write (xern4, '(1pe15.6)') t
            call xermsg ('slatec', 'dlsod',
     *         'in ddebdf, you have changed the value of t from ' //
     *         xern3 // ' to ' // xern4 //
     *         '  this is not allowed on continuation calls.', 10, 1)
            idid=-33
         endif
c
         if (init .ne. 1) then
            if (delsgn*(tout-t) .lt. 0.0d0) then
               write (xern3, '(1pe15.6)') tout
               call xermsg ('slatec', 'dlsod',
     *            'in ddebdf, by calling the code with tout = ' //
     *            xern3 // ' you are attempting to change the ' //
     *            'direction of integration.$$this is not allowed ' //
     *            'without restarting.', 11, 1)
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
            call xermsg ('slatec', 'dlsod',
     *         'in ddebdf, invalid input was detected on ' //
     *         'successive entries.  it is impossible to proceed ' //
     *         'because you have not corrected the problem, ' //
     *         'so execution is being terminated.', 12, 2)
         endif
         return
      endif
c
c        ...............................................................
c
c             rtol = atol = 0. is allowed as valid input and interpreted
c             as asking for the most accurate solution possible. in this
c             case, the relative error tolerance rtol is reset to the
c             smallest value 100*u which is likely to be reasonable for
c             this method and machine
c
      do 180 k = 1, neq
         if (rtol(k) + atol(k) .gt. 0.0d0) go to 170
            rtol(k) = 100.0d0*u
            idid = -2
  170    continue
c     ...exit
         if (itol .eq. 0) go to 190
  180 continue
  190 continue
c
      if (idid .ne. (-2)) go to 200
c        rtol=atol=0 on input, so rtol is changed to a
c                                 small positive value
         ibegin = -1
      go to 460
  200 continue
c        begin block permitting ...exits to 450
c           begin block permitting ...exits to 430
c              begin block permitting ...exits to 260
c                 begin block permitting ...exits to 230
c
c                    branch on status of initialization indicator
c                           init=0 means initial derivatives and
c                           nominal step size
c                                  and direction not yet set
c                           init=1 means nominal step size and
c                           direction not yet set init=2 means no
c                           further initialization required
c
                     if (init .eq. 0) go to 210
c                 ......exit
                        if (init .eq. 1) go to 230
c              .........exit
                        go to 260
  210                continue
c
c                    ................................................
c
c                         more initialization --
c                                             -- evaluate initial
c                                             derivatives
c
                     init = 1
                     call df(t,y,yh(1,2),rpar,ipar)
                     nfe = 1
c                 ...exit
                     if (t .ne. tout) go to 230
                     idid = 2
                     do 220 l = 1, neq
                        ypout(l) = yh(l,2)
  220                continue
                     told = t
c        ............exit
                     go to 450
  230             continue
c
c                 -- compute initial step size
c                 -- save sign of integration direction
c                 -- set independent and dependent variables
c                                      x and yh(*) for dstod
c
                  ltol = 1
                  do 240 l = 1, neq
                     if (itol .eq. 1) ltol = l
                     tol = rtol(ltol)*abs(y(l)) + atol(ltol)
                     if (tol .eq. 0.0d0) go to 390
                     ewt(l) = tol
  240             continue
c
                  big = sqrt(d1mach(2))
                  call dhstrt(df,neq,t,tout,y,yh(1,2),ewt,1,u,big,
     1                        yh(1,3),yh(1,4),yh(1,5),yh(1,6),rpar,
     2                        ipar,h)
c
                  delsgn = sign(1.0d0,tout-t)
                  x = t
                  do 250 l = 1, neq
                     yh(l,1) = y(l)
                     yh(l,2) = h*yh(l,2)
  250             continue
                  init = 2
  260          continue
c
c              ......................................................
c
c                 on each call set information which determines the
c                 allowed interval of integration before returning
c                 with an answer at tout
c
               del = tout - t
               absdel = abs(del)
c
c              ......................................................
c
c                 if already past output point, interpolate and
c                 return
c
  270          continue
c                 begin block permitting ...exits to 400
c                    begin block permitting ...exits to 380
                        if (abs(x-t) .lt. absdel) go to 290
                           call dintyd(tout,0,yh,neq,y,intflg)
                           call dintyd(tout,1,yh,neq,ypout,intflg)
                           idid = 3
                           if (x .ne. tout) go to 280
                              idid = 2
                              intout = .false.
  280                      continue
                           t = tout
                           told = t
c        ..................exit
                           go to 450
  290                   continue
c
c                       if cannot go past tstop and sufficiently
c                       close, extrapolate and return
c
                        if (itstop .ne. 1) go to 310
                        if (abs(tstop-x) .ge. 100.0d0*u*abs(x))
     1                     go to 310
                           dt = tout - x
                           do 300 l = 1, neq
                              y(l) = yh(l,1) + (dt/h)*yh(l,2)
  300                      continue
                           call df(tout,y,ypout,rpar,ipar)
                           nfe = nfe + 1
                           idid = 3
                           t = tout
                           told = t
c        ..................exit
                           go to 450
  310                   continue
c
                        if (iinteg .eq. 0 .or. .not.intout) go to 320
c
c                          intermediate-output mode
c
                           idid = 1
                        go to 370
  320                   continue
c
c                       .............................................
c
c                            monitor number of steps attempted
c
                        if (ksteps .le. maxnum) go to 330
c
c                          a significant amount of work has been
c                          expended
                           idid = -1
                           ksteps = 0
                           ibegin = -1
                        go to 370
  330                   continue
c
c                          ..........................................
c
c                             limit step size and set weight vector
c
                           hmin = 100.0d0*u*abs(x)
                           ha = max(abs(h),hmin)
                           if (itstop .eq. 1)
     1                        ha = min(ha,abs(tstop-x))
                           h = sign(ha,h)
                           ltol = 1
                           do 340 l = 1, neq
                              if (itol .eq. 1) ltol = l
                              ewt(l) = rtol(ltol)*abs(yh(l,1))
     1                                 + atol(ltol)
c                    .........exit
                              if (ewt(l) .le. 0.0d0) go to 380
  340                      continue
                           tolfac = u*dvnrms(neq,yh,ewt)
c                 .........exit
                           if (tolfac .le. 1.0d0) go to 400
c
c                          tolerances too small
                           idid = -2
                           tolfac = 2.0d0*tolfac
                           rtol(1) = tolfac*rtol(1)
                           atol(1) = tolfac*atol(1)
                           if (itol .eq. 0) go to 360
                              do 350 l = 2, neq
                                 rtol(l) = tolfac*rtol(l)
                                 atol(l) = tolfac*atol(l)
  350                         continue
  360                      continue
                           ibegin = -1
  370                   continue
c           ............exit
                        go to 430
  380                continue
c
c                    relative error criterion inappropriate
  390                continue
                     idid = -3
                     ibegin = -1
c           .........exit
                     go to 430
  400             continue
c
c                 ...................................................
c
c                      take a step
c
                  call dstod(neq,y,yh,neq,yh1,ewt,savf,acor,wm,iwm,
     1                       df,djac,rpar,ipar)
c
                  jstart = -2
                  intout = .true.
               if (kflag .eq. 0) go to 270
c
c              ......................................................
c
               if (kflag .eq. -1) go to 410
c
c                 repeated corrector convergence failures
                  idid = -6
                  ibegin = -1
               go to 420
  410          continue
c
c                 repeated error test failures
                  idid = -7
                  ibegin = -1
  420          continue
  430       continue
c
c           .........................................................
c
c                                  store values before returning to
c                                  ddebdf
            do 440 l = 1, neq
               y(l) = yh(l,1)
               ypout(l) = yh(l,2)/h
  440       continue
            t = x
            told = t
            intout = .false.
  450    continue
  460 continue
      return
      end
