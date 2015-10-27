*deck ddntl
      subroutine ddntl (eps, f, fa, hmax, hold, impl, jtask, matdim,
     8   maxord, mint, miter, ml, mu, n, nde, save1, t, uround, users,
     8   y, ywt, h, mntold, mtrold, nfe, rc, yh, a, convrg, el, fac,
     8   ier, ipvt, nq, nwait, rh, rmax, save2, tq, trend, iswflg,
     8   jstate)
c***begin prologue  ddntl
c***subsidiary
c***purpose  subroutine ddntl is called to set parameters on the first
c            call to ddstp, on an internal restart, or when the user has
c            altered mint, miter, and/or h.
c***library   slatec (sdrive)
c***type      double precision (sdntl-s, ddntl-d, cdntl-c)
c***author  kahaner, d. k., (nist)
c             national institute of standards and technology
c             gaithersburg, md  20899
c           sutherland, c. d., (lanl)
c             mail stop d466
c             los alamos national laboratory
c             los alamos, nm  87545
c***description
c
c  on the first call, the order is set to 1 and the initial derivatives
c  are calculated.  rmax is the maximum ratio by which h can be
c  increased in one step.  it is initially rminit to compensate
c  for the small initial h, but then is normally equal to rmnorm.
c  if a failure occurs (in corrector convergence or error test), rmax
c  is set at rmfail for the next increase.
c  if the caller has changed mint, or if jtask = 0, ddcst is called
c  to set the coefficients of the method.  if the caller has changed h,
c  yh must be rescaled.  if h or mint has been changed, nwait is
c  reset to nq + 2 to prevent further increases in h for that many
c  steps.  also, rc is reset.  rc is the ratio of new to old values of
c  the coefficient l(0)*h.  if the caller has changed miter, rc is
c  set to 0 to force the partials to be updated, if partials are used.
c
c***routines called  ddcst, ddscl, dgbfa, dgbsl, dgefa, dgesl, dnrm2
c***revision history  (yymmdd)
c   790601  date written
c   900329  initial submission to slatec.
c***end prologue  ddntl
      integer i, iflag, impl, info, iswflg, jstate, jtask, matdim,
     8        maxord, mint, miter, ml, mntold, mtrold, mu, n, nde, nfe,
     8        nq, nwait
      double precision a(matdim,*), el(13,12), eps, fac(*), h, hmax,
     8     hold, oldl0, rc, rh, rmax, rminit, save1(*), save2(*), dnrm2,
     8     sum, t, tq(3,12), trend, uround, y(*), yh(n,*), ywt(*)
      integer ipvt(*)
      logical convrg, ier
      parameter(rminit = 10000.d0)
c***first executable statement  ddntl
      ier = .false.
      if (jtask .ge. 0) then
        if (jtask .eq. 0) then
          call ddcst (maxord, mint, iswflg,  el, tq)
          rmax = rminit
        end if
        rc = 0.d0
        convrg = .false.
        trend = 1.d0
        nq = 1
        nwait = 3
        call f (n, t, y, save2)
        if (n .eq. 0) then
          jstate = 6
          return
        end if
        nfe = nfe + 1
        if (impl .ne. 0) then
          if (miter .eq. 3) then
            iflag = 0
            call users (y, yh, ywt, save1, save2, t, h, el, impl, n,
     8                  nde, iflag)
            if (iflag .eq. -1) then
              ier = .true.
              return
            end if
            if (n .eq. 0) then
              jstate = 10
              return
            end if
          else if (impl .eq. 1) then
            if (miter .eq. 1 .or. miter .eq. 2) then
              call fa (n, t, y, a, matdim, ml, mu, nde)
              if (n .eq. 0) then
                jstate = 9
                return
              end if
              call dgefa (a, matdim, n, ipvt, info)
              if (info .ne. 0) then
                ier = .true.
                return
              end if
              call dgesl (a, matdim, n, ipvt, save2, 0)
            else if (miter .eq. 4 .or. miter .eq. 5) then
              call fa (n, t, y, a(ml+1,1), matdim, ml, mu, nde)
              if (n .eq. 0) then
                jstate = 9
                return
              end if
              call dgbfa (a, matdim, n, ml, mu, ipvt, info)
              if (info .ne. 0) then
                ier = .true.
                return
              end if
              call dgbsl (a, matdim, n, ml, mu, ipvt, save2, 0)
            end if
          else if (impl .eq. 2) then
            call fa (n, t, y, a, matdim, ml, mu, nde)
            if (n .eq. 0) then
              jstate = 9
              return
            end if
            do 150 i = 1,nde
              if (a(i,1) .eq. 0.d0) then
                ier = .true.
                return
              else
                save2(i) = save2(i)/a(i,1)
              end if
 150          continue
            do 155 i = nde+1,n
 155          a(i,1) = 0.d0
          else if (impl .eq. 3) then
            if (miter .eq. 1 .or. miter .eq. 2) then
              call fa (n, t, y, a, matdim, ml, mu, nde)
              if (n .eq. 0) then
                jstate = 9
                return
              end if
              call dgefa (a, matdim, nde, ipvt, info)
              if (info .ne. 0) then
                ier = .true.
                return
              end if
              call dgesl (a, matdim, nde, ipvt, save2, 0)
            else if (miter .eq. 4 .or. miter .eq. 5) then
              call fa (n, t, y, a(ml+1,1), matdim, ml, mu, nde)
              if (n .eq. 0) then
                jstate = 9
                return
              end if
              call dgbfa (a, matdim, nde, ml, mu, ipvt, info)
              if (info .ne. 0) then
                ier = .true.
                return
              end if
              call dgbsl (a, matdim, nde, ml, mu, ipvt, save2, 0)
            end if
          end if
        end if
        do 170 i = 1,nde
 170      save1(i) = save2(i)/max(1.d0, ywt(i))
        sum = dnrm2(nde, save1, 1)/sqrt(dble(nde))
        if (sum .gt. eps/abs(h)) h = sign(eps/sum, h)
        do 180 i = 1,n
 180      yh(i,2) = h*save2(i)
        if (miter .eq. 2 .or. miter .eq. 5 .or. iswflg .eq. 3) then
          do 20 i = 1,n
 20         fac(i) = sqrt(uround)
        end if
      else
        if (miter .ne. mtrold) then
          mtrold = miter
          rc = 0.d0
          convrg = .false.
        end if
        if (mint .ne. mntold) then
          mntold = mint
          oldl0 = el(1,nq)
          call ddcst (maxord, mint, iswflg,  el, tq)
          rc = rc*el(1,nq)/oldl0
          nwait = nq + 2
        end if
        if (h .ne. hold) then
          nwait = nq + 2
          rh = h/hold
          call ddscl (hmax, n, nq, rmax,  hold, rc, rh, yh)
        end if
      end if
      return
      end
