*deck cdstp
      subroutine cdstp (eps, f, fa, hmax, impl, ierror, jacobn, matdim,
     8   maxord, mint, miter, ml, mu, n, nde, ywt, uround, users, avgh,
     8   avgord, h, hused, jtask, mntold, mtrold, nfe, nje, nqused,
     8   nstep, t, y, yh, a, convrg, dfdy, el, fac, hold, ipvt, jstate,
     8   jstepl, nq, nwait, rc, rmax, save1, save2, tq, trend, iswflg,
     8   mtrsv, mxrdsv)
c***begin prologue  cdstp
c***subsidiary
c***purpose  cdstp performs one step of the integration of an initial
c            value problem for a system of ordinary differential
c            equations.
c***library   slatec (sdrive)
c***type      complex (sdstp-s, ddstp-d, cdstp-c)
c***author  kahaner, d. k., (nist)
c             national institute of standards and technology
c             gaithersburg, md  20899
c           sutherland, c. d., (lanl)
c             mail stop d466
c             los alamos national laboratory
c             los alamos, nm  87545
c***description
c
c  communication with cdstp is done with the following variables:
c
c    yh      an n by maxord+1 array containing the dependent variables
c              and their scaled derivatives.  maxord, the maximum order
c              used, is currently 12 for the adams methods and 5 for the
c              gear methods.  yh(i,j+1) contains the j-th derivative of
c              y(i), scaled by h**j/factorial(j).  only y(i),
c              1 .le. i .le. n, need be set by the calling program on
c              the first entry.  the yh array should not be altered by
c              the calling program.  when referencing yh as a
c              2-dimensional array, use a column length of n, as this is
c              the value used in cdstp.
c    dfdy    a block of locations used for partial derivatives if miter
c              is not 0.  if miter is 1 or 2 its length must be at least
c              n*n.  if miter is 4 or 5 its length must be at least
c              (2*ml+mu+1)*n.
c    ywt     an array of n locations used in convergence and error tests
c    save1
c    save2   arrays of length n used for temporary storage.
c    ipvt    an integer array of length n used by the linear system
c              solvers for the storage of row interchange information.
c    a       a block of locations used to store the matrix a, when using
c              the implicit method.  if impl is 1, a is a matdim by n
c              array.  if miter is 1 or 2 matdim is n, and if miter is 4
c              or 5 matdim is 2*ml+mu+1.  if impl is 2 its length is n.
c              if impl is 3, a is a matdim by nde array.
c    jtask   an integer used on input.
c              it has the following values and meanings:
c                 .eq. 0  perform the first step.  this value enables
c                         the subroutine to initialize itself.
c                .gt. 0  take a new step continuing from the last.
c                         assumes the last step was successful and
c                         user has not changed any parameters.
c                 .lt. 0  take a new step with a new value of h and/or
c                         mint and/or miter.
c    jstate  a completion code with the following meanings:
c                1  the step was successful.
c                2  a solution could not be obtained with h .ne. 0.
c                3  a solution was not obtained in mxtry attempts.
c                4  for impl .ne. 0, the matrix a is singular.
c              on a return with jstate .gt. 1, the values of t and
c              the yh array are as of the beginning of the last
c              step, and h is the last step size attempted.
c
c***routines called  cdcor, cdcst, cdntl, cdpsc, cdpst, cdscl, scnrm2
c***revision history  (yymmdd)
c   790601  date written
c   900329  initial submission to slatec.
c***end prologue  cdstp
      external f, jacobn, fa, users
      integer i, ierror, impl, ipvt(*), iswflg, iter, j, jstate, jstepl,
     8        jtask, matdim, maxord, mint, miter, ml, mntold, mtrold,
     8        mtrsv, mu, mxfail, mxiter, mxrdsv, mxtry, n, nde, ndjstp,
     8        nfail, nfe, nje, nq, nqused, nstep, nsv, ntry, nwait
      complex a(matdim,*), dfdy(matdim,*), fac(*), save1(*), save2(*),
     8        y(*), yh(n,*), ywt(*)
      real avgh, avgord, bias1, bias2, bias3, bnd, ctest, d, denom, d1,
     8     el(13,12), eps, erdn, erup, etest, h, hmax, hn, hold, hs,
     8     hused, numer, rc, rctest, rh, rh1, rh2, rh3, rmax, rmfail,
     8     rmnorm, scnrm2, t, told, tq(3,12), trend, trshld, uround,
     8     y0nrm
      logical convrg, evalfa, evaljc, ier, switch
      parameter(bias1 = 1.3e0, bias2 = 1.2e0, bias3 = 1.4e0, mxfail = 3,
     8          mxiter = 3, mxtry = 50, rctest = .3e0, rmfail = 2.e0,
     8          rmnorm = 10.e0, trshld = 1.e0)
      parameter (ndjstp = 10)
      data ier /.false./
c***first executable statement  cdstp
      nsv = n
      bnd = 0.e0
      switch = .false.
      ntry = 0
      told = t
      nfail = 0
      if (jtask .le. 0) then
        call cdntl (eps, f, fa, hmax, hold, impl, jtask, matdim,
     8              maxord, mint, miter, ml, mu, n, nde, save1, t,
     8              uround, users, y, ywt,  h, mntold, mtrold, nfe, rc,
     8              yh,  a, convrg, el, fac, ier, ipvt, nq, nwait, rh,
     8              rmax, save2, tq, trend, iswflg, jstate)
        if (n .eq. 0) go to 440
        if (h .eq. 0.e0) go to 400
        if (ier) go to 420
      end if
 100  ntry = ntry + 1
      if (ntry .gt. mxtry) go to 410
      t = t + h
      call cdpsc (1, n, nq,  yh)
      evaljc = (((abs(rc - 1.e0) .gt. rctest) .or.
     8  (nstep .ge. jstepl + ndjstp)) .and. (miter .ne. 0))
      evalfa = .not. evaljc
c
 110  iter = 0
      do 115 i = 1,n
 115    y(i) = yh(i,1)
      call f (n, t, y, save2)
      if (n .eq. 0) then
        jstate = 6
        go to 430
      end if
      nfe = nfe + 1
      if (evaljc .or. ier) then
        call cdpst (el, f, fa, h, impl, jacobn, matdim, miter, ml,
     8              mu, n, nde, nq, save2, t, users, y, yh, ywt, uround,
     8              nfe, nje,  a, dfdy, fac, ier, ipvt, save1, iswflg,
     8              bnd, jstate)
        if (n .eq. 0) go to 430
        if (ier) go to 160
        convrg = .false.
        rc = 1.e0
        jstepl = nstep
      end if
      do 125 i = 1,n
 125    save1(i) = 0.e0
c                      up to mxiter corrector iterations are taken.
c                      convergence is tested by requiring the r.m.s.
c                      norm of changes to be less than eps.  the sum of
c                      the corrections is accumulated in the vector
c                      save1(i).  it is approximately equal to the l-th
c                      derivative of y multiplied by
c                      h**l/(factorial(l-1)*el(l,nq)), and is thus
c                      proportional to the actual errors to the lowest
c                      power of h present (h**l).  the yh array is not
c                      altered in the correction loop.  the norm of the
c                      iterate difference is stored in d.  if
c                      iter .gt. 0, an estimate of the convergence rate
c                      constant is stored in trend, and this is used in
c                      the convergence test.
c
 130  call cdcor (dfdy, el, fa, h, ierror, impl, ipvt, matdim, miter,
     8            ml, mu, n, nde, nq, t, users, y, yh, ywt,  evalfa,
     8            save1, save2,  a, d, jstate)
        if (n .eq. 0) go to 430
      if (iswflg .eq. 3 .and. mint .eq. 1) then
        if (iter .eq. 0) then
          numer = scnrm2(n, save1, 1)
          do 132 i = 1,n
 132        dfdy(1,i) = save1(i)
          y0nrm = scnrm2(n, yh, 1)
        else
          denom = numer
          do 134 i = 1,n
 134        dfdy(1,i) = save1(i) - dfdy(1,i)
          numer = scnrm2(n, dfdy, matdim)
          if (el(1,nq)*numer .le. 100.e0*uround*y0nrm) then
            if (rmax .eq. rmfail) then
              switch = .true.
              go to 170
            end if
          end if
          do 136 i = 1,n
 136        dfdy(1,i) = save1(i)
          if (denom .ne. 0.e0)
     8    bnd = max(bnd, numer/(denom*abs(h)*el(1,nq)))
        end if
      end if
      if (iter .gt. 0) trend = max(.9e0*trend, d/d1)
      d1 = d
      ctest = min(2.e0*trend, 1.e0)*d
      if (ctest .le. eps) go to 170
      iter = iter + 1
      if (iter .lt. mxiter) then
        do 140 i = 1,n
 140      y(i) = yh(i,1) + el(1,nq)*save1(i)
        call f (n, t, y, save2)
        if (n .eq. 0) then
          jstate = 6
          go to 430
        end if
        nfe = nfe + 1
        go to 130
      end if
c                     the corrector iteration failed to converge in
c                     mxiter tries.  if partials are involved but are
c                     not up to date, they are reevaluated for the next
c                     try.  otherwise the yh array is retracted to its
c                     values before prediction, and h is reduced, if
c                     possible.  if not, a no-convergence exit is taken.
      if (convrg) then
        evaljc = .true.
        evalfa = .false.
        go to 110
      end if
 160  t = told
      call cdpsc (-1, n, nq,  yh)
      nwait = nq + 2
      if (jtask .ne. 0 .and. jtask .ne. 2) rmax = rmfail
      if (iter .eq. 0) then
        rh = .3e0
      else
        rh = .9e0*(eps/ctest)**(.2e0)
      end if
      if (rh*h .eq. 0.e0) go to 400
      call cdscl (hmax, n, nq, rmax,  h, rc, rh, yh)
      go to 100
c                          the corrector has converged.  convrg is set
c                          to .true. if partial derivatives were used,
c                          to indicate that they may need updating on
c                          subsequent steps.  the error test is made.
 170  convrg = (miter .ne. 0)
      if (ierror .eq. 1 .or. ierror .eq. 5) then
        do 180 i = 1,nde
 180      save2(i) = save1(i)/ywt(i)
      else
        do 185 i = 1,nde
 185      save2(i) = save1(i)/max(abs(y(i)), abs(ywt(i)))
      end if
      etest = scnrm2(nde, save2, 1)/(tq(2,nq)*sqrt(real(nde)))
c
c                           the error test failed.  nfail keeps track of
c                           multiple failures.  restore t and the yh
c                           array to their previous values, and prepare
c                           to try the step again.  compute the optimum
c                           step size for this or one lower order.
      if (etest .gt. eps) then
        t = told
        call cdpsc (-1, n, nq,  yh)
        nfail = nfail + 1
        if (nfail .lt. mxfail .or. nq .eq. 1) then
          if (jtask .ne. 0 .and. jtask .ne. 2) rmax = rmfail
          rh2 = 1.e0/(bias2*(etest/eps)**(1.e0/(nq+1)))
          if (nq .gt. 1) then
            if (ierror .eq. 1 .or. ierror .eq. 5) then
              do 190 i = 1,nde
 190            save2(i) = yh(i,nq+1)/ywt(i)
            else
              do 195 i = 1,nde
 195            save2(i) = yh(i,nq+1)/max(abs(y(i)), abs(ywt(i)))
            end if
            erdn = scnrm2(nde, save2, 1)/(tq(1,nq)*sqrt(real(nde)))
            rh1 = 1.e0/max(1.e0, bias1*(erdn/eps)**(1.e0/nq))
            if (rh2 .lt. rh1) then
              nq = nq - 1
              rc = rc*el(1,nq)/el(1,nq+1)
              rh = rh1
            else
              rh = rh2
            end if
          else
            rh = rh2
          end if
          nwait = nq + 2
          if (rh*h .eq. 0.e0) go to 400
          call cdscl (hmax, n, nq, rmax,  h, rc, rh, yh)
          go to 100
        end if
c                control reaches this section if the error test has
c                failed mxfail or more times.  it is assumed that the
c                derivatives that have accumulated in the yh array have
c                errors of the wrong order.  hence the first derivative
c                is recomputed, the order is set to 1, and the step is
c                retried.
        nfail = 0
        jtask = 2
        do 215 i = 1,n
 215      y(i) = yh(i,1)
        call cdntl (eps, f, fa, hmax, hold, impl, jtask, matdim,
     8              maxord, mint, miter, ml, mu, n, nde, save1, t,
     8              uround, users, y, ywt,  h, mntold, mtrold, nfe, rc,
     8              yh,  a, convrg, el, fac, ier, ipvt, nq, nwait, rh,
     8              rmax, save2, tq, trend, iswflg, jstate)
        rmax = rmnorm
        if (n .eq. 0) go to 440
        if (h .eq. 0.e0) go to 400
        if (ier) go to 420
        go to 100
      end if
c                          after a successful step, update the yh array.
      nstep = nstep + 1
      hused = h
      nqused = nq
      avgh = ((nstep-1)*avgh + h)/nstep
      avgord = ((nstep-1)*avgord + nq)/nstep
      do 230 j = 1,nq+1
        do 230 i = 1,n
 230      yh(i,j) = yh(i,j) + el(j,nq)*save1(i)
      do 235 i = 1,n
 235    y(i) = yh(i,1)
c                                          if iswflg is 3, consider
c                                          changing integration methods.
      if (iswflg .eq. 3) then
        if (bnd .ne. 0.e0) then
          if (mint .eq. 1 .and. nq .le. 5) then
            hn = abs(h)/max(uround, (etest/eps)**(1.e0/(nq+1)))
            hn = min(hn, 1.e0/(2.e0*el(1,nq)*bnd))
            hs = abs(h)/max(uround,
     8      (etest/(eps*el(nq+1,1)))**(1.e0/(nq+1)))
            if (hs .gt. 1.2e0*hn) then
              mint = 2
              mntold = mint
              miter = mtrsv
              mtrold = miter
              maxord = min(mxrdsv, 5)
              rc = 0.e0
              rmax = rmnorm
              trend = 1.e0
              call cdcst (maxord, mint, iswflg, el, tq)
              nwait = nq + 2
            end if
          else if (mint .eq. 2) then
            hs = abs(h)/max(uround, (etest/eps)**(1.e0/(nq+1)))
            hn = abs(h)/max(uround,
     8      (etest*el(nq+1,1)/eps)**(1.e0/(nq+1)))
            hn = min(hn, 1.e0/(2.e0*el(1,nq)*bnd))
            if (hn .ge. hs) then
              mint = 1
              mntold = mint
              miter = 0
              mtrold = miter
              maxord = min(mxrdsv, 12)
              rmax = rmnorm
              trend = 1.e0
              convrg = .false.
              call cdcst (maxord, mint, iswflg, el, tq)
              nwait = nq + 2
            end if
          end if
        end if
      end if
      if (switch) then
        mint = 2
        mntold = mint
        miter = mtrsv
        mtrold = miter
        maxord = min(mxrdsv, 5)
        nq = min(nq, maxord)
        rc = 0.e0
        rmax = rmnorm
        trend = 1.e0
        call cdcst (maxord, mint, iswflg, el, tq)
        nwait = nq + 2
      end if
c                           consider changing h if nwait = 1.  otherwise
c                           decrease nwait by 1.  if nwait is then 1 and
c                           nq.lt.maxord, then save1 is saved for use in
c                           a possible order increase on the next step.
c
      if (jtask .eq. 0 .or. jtask .eq. 2) then
        rh = 1.e0/max(uround, bias2*(etest/eps)**(1.e0/(nq+1)))
        if (rh.gt.trshld) call cdscl (hmax, n, nq, rmax, h, rc, rh, yh)
      else if (nwait .gt. 1) then
        nwait = nwait - 1
        if (nwait .eq. 1 .and. nq .lt. maxord) then
          do 250 i = 1,nde
 250        yh(i,maxord+1) = save1(i)
        end if
c             if a change in h is considered, an increase or decrease in
c             order by one is considered also.  a change in h is made
c             only if it is by a factor of at least trshld.  factors
c             rh1, rh2, and rh3 are computed, by which h could be
c             multiplied at order nq - 1, order nq, or order nq + 1,
c             respectively.  the largest of these is determined and the
c             new order chosen accordingly.  if the order is to be
c             increased, we compute one additional scaled derivative.
c             if there is a change of order, reset nq and the
c             coefficients.  in any case h is reset according to rh and
c             the yh array is rescaled.
      else
        if (nq .eq. 1) then
          rh1 = 0.e0
        else
          if (ierror .eq. 1 .or. ierror .eq. 5) then
            do 270 i = 1,nde
 270          save2(i) = yh(i,nq+1)/ywt(i)
          else
            do 275 i = 1,nde
 275          save2(i) = yh(i,nq+1)/max(abs(y(i)), abs(ywt(i)))
          end if
          erdn = scnrm2(nde, save2, 1)/(tq(1,nq)*sqrt(real(nde)))
          rh1 = 1.e0/max(uround, bias1*(erdn/eps)**(1.e0/nq))
        end if
        rh2 = 1.e0/max(uround, bias2*(etest/eps)**(1.e0/(nq+1)))
        if (nq .eq. maxord) then
          rh3 = 0.e0
        else
          if (ierror .eq. 1 .or. ierror .eq. 5) then
            do 290 i = 1,nde
 290          save2(i) = (save1(i) - yh(i,maxord+1))/ywt(i)
          else
            do 295 i = 1,nde
              save2(i) = (save1(i) - yh(i,maxord+1))/
     8        max(abs(y(i)), abs(ywt(i)))
 295          continue
          end if
          erup = scnrm2(nde, save2, 1)/(tq(3,nq)*sqrt(real(nde)))
          rh3 = 1.e0/max(uround, bias3*(erup/eps)**(1.e0/(nq+2)))
        end if
        if (rh1 .gt. rh2 .and. rh1 .ge. rh3) then
          rh = rh1
          if (rh .le. trshld) go to 380
          nq = nq - 1
          rc = rc*el(1,nq)/el(1,nq+1)
        else if (rh2 .ge. rh1 .and. rh2 .ge. rh3) then
          rh = rh2
          if (rh .le. trshld) go to 380
        else
          rh = rh3
          if (rh .le. trshld) go to 380
          do 360 i = 1,n
 360        yh(i,nq+2) = save1(i)*el(nq+1,nq)/(nq+1)
          nq = nq + 1
          rc = rc*el(1,nq)/el(1,nq-1)
        end if
        if (iswflg .eq. 3 .and. mint .eq. 1) then
          if (bnd.ne.0.e0) rh = min(rh, 1.e0/(2.e0*el(1,nq)*bnd*abs(h)))
        end if
        call cdscl (hmax, n, nq, rmax,  h, rc, rh, yh)
        rmax = rmnorm
 380    nwait = nq + 2
      end if
c               all returns are made through this section.  h is saved
c               in hold to allow the caller to change h on the next step
      jstate = 1
      hold = h
      return
c
 400  jstate = 2
      hold = h
      do 405 i = 1,n
 405    y(i) = yh(i,1)
      return
c
 410  jstate = 3
      hold = h
      return
c
 420  jstate = 4
      hold = h
      return
c
 430  t = told
      call cdpsc (-1, nsv, nq,  yh)
      do 435 i = 1,nsv
 435    y(i) = yh(i,1)
 440  hold = h
      return
      end
