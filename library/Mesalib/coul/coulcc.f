*deck coulcc
      subroutine coulcc(xx,eta1,zlmin,nl, fc,gc,fcp,gcp, sig,
     x                  mode1,kfn,ifail)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c  complex coulomb wavefunction program using steed's method           c
c                                                                      c
c  a. r. barnett           manchester  march   1981                    c
c  modified i.j. thompson  daresbury, sept. 1983 for complex functions c
c                                                                      c
c  original program  rcwfn       in    cpc  8 (1974) 377-395           c
c                 +  rcwff       in    cpc 11 (1976) 141-142           c
c                 +  coulfg      in    cpc 27 (1982) 147-166           c
c  description of real algorithm in    cpc 21 (1981) 297-314           c
c  description of complex algorithm    cpc xx (1984) yyy-zzz           c
c  this version written up       in    cpc xx (1984) yyy-zzz           c
c                                                                      c
c  coulcc returns f,g,g',g',sig for complex xx, eta1, and zlmin,       c
c   for nl integer-spaced lambda values zlmin to zlmin+nl-1 inclusive, c
c   thus giving  complex-energy solutions to the coulomb schrodinger   c
c   equation,to the klein-gordon equation and to suitable forms of     c
c   the dirac equation ,also spherical & cylindrical bessel equations  c
c                                                                      c
c  if /mode1/= 1  get f,g,f',g'   for integer-spaced lambda values     c
c            = 2      f,g      unused arrays must be dimensioned in    c
c            = 3      f,  f'          call to at least length (1)      c
c            = 4      f                                                c
c            = 11 get f,h+,f',h+' ) if kfn=0, h+ = g + i.f        )    c
c            = 12     f,h+        )       >0, h+ = j + i.y = h(1) ) in c
c            = 21 get f,h-,f',h-' ) if kfn=0, h- = g - i.f        ) gc c
c            = 22     f,h-        )       >0, h- = j - i.y = h(2) )    c
c                                                                      c
c     if mode1<0 then the values returned are scaled by an exponential c
c                factor (dependent only on xx) to bring nearer unity   c
c                the functions for large /xx/, small eta & /zl/ < /xx/ c
c        define scale = (  0        if mode1 > 0                       c
c                       (  aimag(xx) if mode1 < 0  &  kfn < 3
c                       (  real(xx) if mode1 < 0  &  kfn = 3           c
c        then fc = exp(-abs(scale)) * ( f, j, j, or i)                 c
c         and gc = exp(-abs(scale)) * ( g, y, or y )                   c
c               or exp(scale)       * ( h+, h(1), or k)                c
c               or exp(-scale)      * ( h- or h(2) )                   c
c                                                                      c
c  if  kfn  =  0,-1  complex coulomb functions are returned   f & g    c
c           =  1   spherical bessel      "      "     "       j & y    c
c           =  2 cylindrical bessel      "      "     "       j & y    c
c           =  3 modified cyl. bessel    "      "     "       i & k    c
c                                                                      c
c          and where coulomb phase shifts put in sig if kfn=0 (not -1) c
c                                                                      c
c  the use of mode and kfn is independent                              c
c    (except that for kfn=3,  h(1) & h(2) are not given)               c
c                                                                      c
c  with negative orders lambda, coulcc can still be used but with      c
c    reduced accuracy as cf1 becomes unstable. the user is thus        c
c    strongly advised to use reflection formulae based on              c
c    h+-(zl,,) = h+-(-zl-1,,) * exp +-i(sig(zl)-sig(-zl-1)-(zl+1/2)pi) c
c                                                                      c
c  precision:  results to within 2-3 decimals of 'machine accuracy',   c
c               but if cf1a fails because x too small or eta too large c
c               the f solution  is less accurate if it decreases with  c
c               decreasing lambda (e.g. for lambda.le.-1 & eta.ne.0)   c
c              rerr in common/steed/ traces the main roundoff errors.  c
c                                                                      c
c   coulcc is coded for real on ibm or equivalent  accur >= 10**-14  c
c          with a section of doubled real*16 for less roundoff errors. c
c          (if no doubled precision available, increase jmax to eg 100)c
c   use implicit complex*32 & real*16 on vs compiler accur >= 10**-32  c
c   for single precision cdc (48 bits) reassign real=real etc.       c
c                                                                      c
c   ifail  on input   = 0 : no printing of error messages              c
c                    ne 0 : print error messages on file 6             c
c   ifail  in output = -2 : argument out of range                      c
c                    = -1 : one of the continued fractions failed,     c
c                           or arithmetic check before final recursion c
c                    =  0 : all calculations satisfactory              c
c                    ge 0 : results available for orders up to & at    c
c                             position nl-ifail in the output arrays.  c
c                    = -3 : values at zlmin not found as over/underflowc
c                    = -4 : roundoff errors make results meaningless   c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c     machine dependent constants :                                    c
c                                                                      c
c     accur    target bound on relative error (except near 0 crossings)c
c               (accur should be at least 100 * acc8)                  c
c     acc8     smallest number with 1+acc8 .ne.1 in real  arithmetic c
c     acc16    smallest number with 1+acc16.ne.1 in real*16 arithmetic c
c     fpmax    magnitude of largest floating point number * acc8       c
c     fpmin    magnitude of smallest floating point number / acc8      c
c     fplmax   log(fpmax)                                              c
c     fplmin   log(fpmin)                                              c
c                                                                      c
c     routines called :       logam/clogam/cdigam,                     c
c                             f20, cf1a, rcf, cf1c, cf2, f11, cf1r     c
c     intrinsic functions :   min, max, sqrt, real, imag, abs, log, exp,
c      (generic names)        nint, mod, atan, atan2, cos, sin, cmplx,
c                             sign, conjg, int, tanh                   c
c     note: statement fntn.   nintc = integer nearest to a complex no. c
c                                                                      c
c     parameters determining region of calculations :                  c
c                                                                      c
c        r20      estimate of (2f0 iterations)/(cf2 iterations)        c
c        asym     minimum x/(eta**2+l) for cf1a to converge easily     c
c        xnear    minimum abs(x) for cf2 to converge accurately        c
c        limit    maximum no. iterations for cf1, cf2, and 1f1 series  c
c        jmax     size of work arrays for pade accelerations           c
c        ndrop    number of successive decrements to define instabilityc
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit complex (a-h,o-z)
      parameter(jmax=50)
      dimension fc(nl),gc(nl),fcp(nl),gcp(nl),sig(nl),xrcf(jmax,4)
      logical pr,etane0,ifcp,rlel,donem,unstab,zlneg,axial,nocf2
      real err,rerr,absc,accur,acct,acc8,acch,acc16,accb, xnear,cf1r,
     x       zero,one,two,half,hpi,tlog,fpmax,fpmin,fplmin,fplmax,
     x       paccq,eps,off,scale,sf,sfsh,ta,rk,omega,r20,asym,absx
c
      common       /steed/ rerr,nfp,n11,npq(2),n20,kas(2)
      common /io/ inp, iout
c***  common blocks are for information & storage only.
c     (they are not essential to working of the code)
      common /rcfcm1/ pk,ek,clgaa,clgab,clgbb,dsig,tpk1,w,rl,fcl1,q,gam,
     x                hcl,hpl,fcm,hcl1,alpha,beta,pl
      equivalence            (pk,xrcf(1,1))
c
      data zero,one,two,limit /0.0e+0, 1.0e+0, 2.0e+0, 20000 /,
     x     half, ci / 0.5e+0, (0e+0, 1e+0) /,
     x     fpmax,fpmin,fplmax,fplmin / 1e+308,1e-308,850e+0,-850e+0/,
     x     r20,asym,xnear,ndrop / 3., 3., .5, 5 /,
     x     accur, acc8, acc16 / 1e-10, 4e-15, 2e-29 /
      nintc(w) = nint(real(w))
      absc(w) = abs(real(w)) + abs(aimag(w))
c
      mode = mod(abs(mode1),10)
      ifcp = mod(mode,2).eq.1
      pr = ifail.ne.0
      ifail = -2
      n11   = 0
      nfp   = 0
      kas(1)   = 0
      kas(2)   = 0
      npq(1)   = 0
      npq(2)   = 0
      n20 = 0
      hpi = two*atan(one)
      tlog = log(two)
      accur = max(accur, 10*acc8)
      acct = accur * .5
c                       initialise the log-gamma function :
      call logam(acc8)
      acch  = sqrt(accur)
      accb  = sqrt(acch)
      rerr = acct
c
      cik = one
         if(kfn.ge.3) cik = ci * sign(one,acc8-aimag(xx))
      x     = xx * cik
      eta   = eta1
      if(kfn .gt. 0) eta = zero
         etane0  = absc(eta).gt.acc8
         etai = eta*ci
      dell  = zero
      if(kfn .ge. 2)  dell = half
      zm1   = zlmin - dell
      scale = zero
      if(mode1.lt.0) scale = aimag(x)
c
      m1 = 1
      l1  = m1 + nl - 1
      rlel = abs(aimag(eta)) + abs(aimag(zm1)) .lt. acc8
      absx = abs(x)
      axial = rlel .and. abs(aimag(x)) .lt. acc8 * absx
      if(mode.le.2 .and. absx.lt.fpmin) go to 310
      xi  = one/x
      xlog = log(x)
c            log with cut along the negative real axis| see also omega
      id = 1
      donem = .false.
         unstab = .false.
      lf = m1
      ifail = -1
   10    zlm = zm1 + lf - m1
         zll = zm1 + l1 - m1
c
c ***       zll  is final lambda value, or 0.5 smaller for j,y bessels
c
              z11 = zll
              if(id.lt.0) z11 = zlm
              p11 = ci*sign(one,acc8-aimag(eta))
      last = l1
c
c ***       find phase shifts and gamow factor at lambda = zll
c
      pk = zll + one
      aa = pk - etai
      ab = pk + etai
      bb = two*pk
         zlneg = real(bb).le.zero .and. abs(bb-nintc(bb)).lt.accb
                     clgaa = clogam(aa)
                     clgab = clgaa
         if(etane0.and..not.rlel)  clgab = clogam(ab)
         if(etane0.and.     rlel)  clgab = conjg(clgaa)
         sigma = (clgaa - clgab) * ci*half
         if(kfn.eq.0) sig(l1) = sigma
         if(.not.zlneg) cll = zll*tlog- hpi*eta - clogam(bb)
     x                                          + (clgaa+clgab)*half
              theta  = x - eta*(xlog+tlog) - zll*hpi + sigma
c
        ta =(aimag(aa)**2+aimag(ab)**2+abs(real(aa))+abs(real(ab)))*half
      if(id.gt.0 .and. absx .lt. ta*asym .and. .not.zlneg) go to 20
c
c ***         use cf1 instead of cf1a, if predicted to converge faster,
c                 (otherwise using cf1a as it treats negative lambda &
c                  recurrence-unstable cases properly)
c
           rk = sign(one, real(x) + acc8)
           p =  theta
           if(rk.lt.0) p = -x + eta*(log(-x)+tlog)-zll*hpi-sigma
      f = rk * cf1a(x*rk,eta*rk,zll,p,acct,jmax,nfp,fest,err,fpmax,xrcf,
     x                                      xrcf(1,3), xrcf(1,4))
      fesl = log(fest) + abs(aimag(x))
         nfp = - nfp
      if(nfp.lt.0   .or.(unstab.and.err.lt.accb)) go to 40
      if(.not.zlneg .or. unstab.and.err.gt.accb)  go to 20
         if(pr) write(iout,1060) '-l',err
         if(err.gt.accb) go to 280
         go to 40
c
c ***    evaluate cf1  =  f   =  f'(zll,eta,x)/f(zll,eta,x)
c
   20 if(axial) then
c                                                        real version
      f = cf1r(x,eta,zll,acc8,sf ,rk,  etane0,limit,err,nfp,
     x         acch,fpmin,fpmax,pr,'coulcc')
          fcl = sf
          tpk1= rk
         else
c                                                        complex version
      f = cf1c(x,eta,zll,acc8,fcl,tpk1,etane0,limit,err,nfp,
     x         acch,fpmin,fpmax,pr,'coulcc')
         endif
      if(err.gt.one) go to 390
c
c ***  make a simple check for cf1 being badly unstable:
c
      if(id.lt.0) go to 30
      unstab = real((one-eta*xi)*ci*aimag(theta)/f).gt.zero
     x .and..not.axial .and. abs(aimag(theta)).gt.-log(acc8)*.5
     x .and. absc(eta)+absc(zll).lt.absc(x)
      if(unstab) go to 60
c
c *** compare accumulated phase fcl with asymptotic phase for g(k+1) :
c     to determine estimate of f(zll) (with correct sign) to start recur
c
   30 w   =  x*x  *(half/tpk1 + one/tpk1**2) + eta*(eta-two*x)/tpk1
      fesl   = (zll+one) * xlog + cll - w - log(fcl)
   40 fesl = fesl - abs(scale)
          rk   =        max(real(fesl), fplmin*half)
          fesl = cmplx(min(rk,   fplmax*half ) , aimag(fesl))
      fest= exp(fesl)
c
           rerr = max(rerr, err, acc8 * abs(real(theta)) )
c
      fcl = fest
      fpl = fcl*f
      if(ifcp) fcp(l1) = fpl
               fc (l1) = fcl
c
c *** downward recurrence to lambda = zlm. array gc,if present,stores rl
c
      i  = max(-id, 0)
      zl  = zll + i
         mono = 0
        off = abs(fcl)
         ta = absc(sigma)
      do 70  l  = l1-id,lf,-id
         if(etane0) then
               if(rlel) then
                    dsig = atan2(real(eta),real(zl))
                    rl = sqrt(real(zl)**2 + real(eta)**2)
                  else
                    aa = zl - etai
                    bb = zl + etai
                    if(absc(aa).lt.acch.or.absc(bb).lt.acch) goto 50
                    dsig = (log(aa) - log(bb)) * ci*half
                    rl = aa * exp(ci*dsig)
                 endif
             if(absc(sigma).lt.ta*half) then
c               re-calculate sigma because of accumulating roundoffs:
                sl =(clogam(zl+i-etai)-clogam(zl+i+etai))*ci*half
                rl = (zl - etai) * exp(ci*id*(sigma - sl))
                sigma = sl
                ta = zero
              else
                sigma = sigma - dsig * id
              endif
                ta = max(ta, absc(sigma))
             sl    =  eta  + zl*zl*xi
                pl = zero
                if(absc(zl).gt.acch) pl = (sl*sl - rl*rl)/zl
             fcl1  = (fcl *sl + id*zl*fpl)/rl
              sf = abs(fcl1)
                       if(sf.gt.fpmax) go to 350
             fpl   = (fpl *sl + id*pl*fcl)/rl
             if(mode .le. 1) gcp(l+id)= pl * id
        else
c                               eta = 0, including bessels.  nb rl==sl
           rl = zl* xi
           fcl1 = fcl * rl + fpl*id
              sf = abs(fcl1)
                      if(sf.gt.fpmax) go to 350
           fpl  =(fcl1* rl - fcl) * id
        endif
c             if(absc(fcl1).lt.absc(fcl)) then
              if(sf.lt.off) then
                 mono = mono + 1
                else
                 mono = 0
                endif
         fcl   =  fcl1
           off = sf
         fc(l) =  fcl
         if(ifcp) fcp(l)  = fpl
           if(kfn.eq.0) sig(l) = sigma
           if(mode .le. 2) gc(l+id) = rl
      zl = zl - id
      if(mono.lt.ndrop) go to 70
      if(axial .or. real(zlm)*id.gt.-ndrop.and..not.etane0) go to 70
         unstab = .true.
c
c ***    take action if cannot or should not recur below this zl:
   50    zlm = zl
         lf = l
            if(id.lt.0) go to 380
         if(.not.unstab) lf = l + 1
         if(l+mono.lt.l1-2 .or. id.lt.0 .or. .not.unstab) go to 80
c             otherwise, all l values (for stability) should be done
c                        in the reverse direction:
         go to 60
   70 continue
      go to 80
   60       id = -1
            lf = l1
            l1 = m1
            rerr = acct
            go to 10
   80 if(fcl .eq. zero) fcl = + acc8
      f  = fpl/fcl
c
c *** check, if second time around, that the 'f' values agree
c
      if(id.gt.0) first = f
      if(donem) rerr = max(rerr, absc(f-first)/absc(f))
      if(donem) go to 90
c
       nocf2 = .false.
      thetam  = x - eta*(xlog+tlog) - zlm*hpi + sigma
c
c *** on left x-plane, determine omega by requiring cut on -x axis
c     on right x-plane, choose omega (using estimate based on thetam)
c       so h(omega) is smaller and recurs upwards accurately.
c     (x-plane boundary is shifted to give cf2(lh) a chance to converge)
c
                           omega = sign(one,aimag(x)+acc8)
      if(real(x).ge.xnear) omega = sign(one,aimag(thetam)+acc8)
c
         sfsh = exp(omega*scale - abs(scale))
         off=exp(min(two * max(abs(aimag(x)),abs(aimag(thetam)),
     x                         abs(aimag(zlm))*3 ) , fplmax) )
          eps = max(acc8 , acct * half / off)
c
c ***    try first estimated omega, then its opposite,
c        to find the h(omega) linearly independent of f
c        i.e. maximise  cf1-cf2 = 1/(f h(omega)) , to minimise h(omega)
c
   90 do 100 l=1,2
         lh = 1
         if(omega.lt.zero) lh = 2
      pm = ci*omega
      etap = eta * pm
         if(donem) go to 130
         pq1 = zero
         paccq = one
         kase = 0
c
c ***            check for small x, i.e. whether to avoid cf2 :
c
      if(mode.ge.3 .and. absx.lt.one ) go to 190
      if(mode.lt.3 .and. (nocf2 .or. absx.lt.xnear .and.
     x   absc(eta)*absx .lt. 5 .and. absc(zlm).lt.4)) then
        kase = 5
        go to 120
        endif
c
c ***  evaluate   cf2 : pq1 = p + i.omega.q  at lambda = zlm
c
         pq1 = cf2(x,eta,zlm,pm,eps,limit,err,npq(lh),acc8,acch,
     x             pr,accur,dell,'coulcc')
c
       err = err * max(one,absc(pq1)/max(absc(f-pq1),acc8) )
       if(err.lt.acch)       go to 110
c
c *** check if impossible to get f-pq accurately because of cancellation
             nocf2 = real(x).lt.xnear .and. abs(aimag(x)).lt.-log(acc8)
c                original guess for omega (based on thetam) was wrong
c                use kase 5 or 6 if necessary if re(x) < xnear
  100            omega = - omega
                if(unstab) go to 360
                if(real(x).lt.-xnear .and. pr) write(iout,1060) '-x',err
  110     rerr = max(rerr,err)
c
c ***  establish case of calculation required for irregular solution
c
  120 if(kase.ge.5) go to 130
      if(real(x) .gt. xnear) then
c          estimate errors if kase 2 or 3 were to be used:
         paccq = eps * off * absc(pq1) / max(abs(aimag(pq1)),acc8)
        endif
      if(paccq .lt. accur) then
          kase = 2
          if(axial) kase = 3
      else
          kase = 1
          if(npq(1) * r20 .lt. jmax)     kase = 4
c             i.e. change to kase=4 if the 2f0 predicted to converge
      endif
  130 go to (190,140,150,170,190,190),  abs(kase)
  140    if(.not.donem)
c
c ***  evaluate   cf2 : pq2 = p - i.omega.q  at lambda = zlm   (kase 2)
c
     x  pq2 = cf2(x,eta,zlm,-pm,eps,limit,err,npq(3-lh),acc8,acch,
     x             pr,accur,dell,'coulcc')
c
        p     = (pq2 + pq1) * half
        q     = (pq2 - pq1) * half*pm
      go to 160
  150   p     = real(pq1)
        q     = aimag(pq1)
c
c ***   with kase = 3 on the real axes, p and q are real & pq2 = pq1*
c
        pq2 = conjg(pq1)
c
c *** solve for fcm = f at lambda = zlm,then find norm factor w=fcm/fcl
c
  160 w   = (pq1 - f) * (pq2 - f)
         sf = exp(-abs(scale))
      fcm = sqrt(q / w) * sf
c                  any sqrt given here is corrected by
c                  using sign for fcm nearest to phase of fcl
      if(real(fcm/fcl).lt.zero) fcm  = - fcm
      gam = (f - p)/q
         ta = absc(gam + pm)
         paccq= eps * max(ta,one/ta)
      hcl = fcm * (gam + pm) * (sfsh/(sf*sf))
c
      if(paccq.gt.max(accur,100*acc8) .and. kase.gt.0) then
c                                    consider a kase = 1 calculation
          f11v= f11(x,eta,z11,p11,acct,limit,0,err,n11,fpmax,acc8,acc16)
          if(err.lt.paccq) go to 200
          endif
      rerr=max(rerr,paccq)
      go to 230
c
c *** arrive here if kase = 4
c     to evaluate the exponentially decreasing h(lh) directly.
c
  170  if(donem) go to 180
      aa = etap - zlm
      bb = etap + zlm + one
      f20v = f20(aa,bb,-half*pm*xi, acct,jmax,err,fpmax,n20,xrcf)
        if(n20.le.0) go to 190
        rerr = max(rerr,err)
         hcl = fpmin
         if(abs(real(pm*thetam)+omega*scale).gt.fplmax) go to 330
  180 hcl = f20v * exp(pm * thetam + omega*scale)
      fcm = sfsh / ((f - pq1) * hcl )
      go to 230
c
c *** arrive here if kase=1   (or if 2f0 tried mistakenly & failed)
c
c           for small values of x, calculate f(x,sl) directly from 1f1
c               using real*16 arithmetic if possible.
c           where z11 = zll if id>0, or = zlm if id<0
c
  190 f11v = f11(x,eta,z11,p11,acct,limit,0,err,n11,fpmax,acc8,acc16)
c
  200       if(n11.lt.0) then
c                               f11 failed from bb = negative integer
               write(iout,1060) '-l',one
               go to 390
               endif
            if(err.gt.paccq .and. paccq.lt.accb) then
c                               consider a kase 2 or 3 calculation :
                kase = -2
                if(axial) kase = -3
                go to 130
                endif
         rerr = max(rerr, err)
         if(err.gt.fpmax) go to 370
         if(id.lt.0) cll = z11*tlog- hpi*eta - clogam(bb)
     x                       + clogam(z11 + one + p11*eta) - p11*sigma
      ek   = (z11+one)*xlog - p11*x + cll  - abs(scale)
      if(id.gt.0) ek = ek - fesl + log(fcl)
         if(real(ek).gt.fplmax) go to 350
         if(real(ek).lt.fplmin) go to 340
      fcm = f11v * exp(ek)
c
      if(kase.ge.5) then
        if(absc(zlm+zlm-nintc(zlm+zlm)).lt.accb) kase = 6
c
c ***  for abs(x) < xnear, then cf2 may not converge accurately, so
c ***      use an expansion for irregular soln from origin :
c
         sl = zlm
            zlneg = real(zlm) .lt. -one + accb
         if(kase.eq.5 .or. zlneg) sl = - zlm - one
         pk = sl + one
            aa = pk - etap
            ab = pk + etap
            bb = two*pk
                     clgaa = clogam(aa)
                     clgab = clgaa
         if(etane0)  clgab = clogam(ab)
                     clgbb = clogam(bb)
          cll = sl*tlog- hpi*eta - clgbb + (clgaa + clgab) * half
          dsig = (clgaa - clgab) * pm*half
             if(kase.eq.6) p11 = - pm
          ek  = pk * xlog - p11*x + cll  - abs(scale)
                     sf = exp(-abs(scale))
                     chi = zero
       if(.not.( kase.eq.5 .or. zlneg ) ) go to 210
c
c *** use  g(l)  =  (cos(chi) * f(l) - f(-l-1)) /  sin(chi)
c
c      where chi = sig(l) - sig(-l-1) - (2l+1)*pi/2
c
         chi = sigma - dsig - (zlm-sl) * hpi
         f11v=f11(x,eta,sl,p11,acct,limit,0,err,npq(1),fpmax,acc8,acc16)
                    rerr = max(rerr,err)
            if(kase.eq.6) go to 210
         fesl = f11v * exp( ek )
         fcl1 = exp(pm*chi) * fcm
         hcl = fcl1 - fesl
               rerr=max(rerr,acct*max(absc(fcl1),absc(fesl))/absc(hcl))
         hcl = hcl / sin(chi) * (sfsh/(sf*sf))
       go to 220
c
c *** use the logarithmic expansion for the irregular solution (kase 6)
c        for the case that bb is integral so sin(chi) would be zero.
c
  210    rl = bb - one
         n  = nintc(rl)
         zlog = xlog + tlog - pm*hpi
         chi = chi + pm * thetam + omega * scale + ab * zlog
            aa  = one - aa
         if(abs(nintc(aa)-aa).lt.accur .and. real(aa).lt.half) then
            hcl = zero
         else
               if(id.gt.0 .and. .not.zlneg) f11v = fcm * exp(-ek)
            hcl = exp(chi - clgbb - clogam(aa)) * (-1)**(n+1)
     x              * ( f11v * zlog +
     x      f11(x,eta,sl,-pm,acct,limit,2,err,npq(2),fpmax,acc8,acc16))
                rerr = max(rerr,err)
            endif
         if(n.gt.0) then
             ek = chi + clogam(rl) - clgab - rl*zlog
             df =f11(x,eta,-sl-one,-pm,zero,n,0,err,l,fpmax,acc8,acc16)
             hcl = hcl + exp(ek) * df
            endif
c
  220    pq1 = f - sfsh/(fcm * hcl)
      else
           if(mode.le.2) hcl = sfsh/((f - pq1) * fcm)
           kase = 1
      endif
c
c ***  now have absolute normalisations for coulomb functions
c          fcm & hcl  at lambda = zlm
c      so determine linear transformations for functions required :
c
  230 ih = abs(mode1) / 10
        if(kfn.eq.3) ih = (3-aimag(cik))/2  + half
      p11 = one
      if(ih.eq.1) p11 = ci
      if(ih.eq.2) p11 = -ci
                  df = - pm
      if(ih.ge.1) df = - pm + p11
          if(absc(df).lt.acch) df = zero
c
c *** normalisations for spherical or cylindrical bessel functions
c
                          alpha = zero
          if(kfn  .eq. 1) alpha = xi
          if(kfn  .ge. 2) alpha = xi*half
                          beta  = one
          if(kfn  .eq. 1) beta  = xi
          if(kfn  .ge. 2) beta  = sqrt(xi/hpi)
          if(kfn  .ge. 2 .and. real(beta).lt.zero) beta  = - beta
c
      aa = one
      if(kfn.gt.0) aa = -p11 * beta
      if(kfn.ge.3) then
c                        calculate rescaling factors for i & k output
         p = exp((zlm+dell) * hpi * cik)
         aa= beta * hpi * p
         beta = beta / p
         q = cik * id
        endif
c                        calculate rescaling factors for gc output
      if(ih.eq.0) then
         ta = abs(scale) + aimag(pm)*scale
         rk = zero
         if(ta.lt.fplmax) rk = exp(-ta)
       else
         ta = abs(scale) + aimag(p11)*scale
c
         if(absc(df).gt.acch .and. ta.gt.fplmax) go to 320
         if(absc(df).gt.acch) df = df * exp(ta)
         rk = exp(two * (lh-ih) * scale)
      endif
c
         kas((3-id)/2) = kase
      w = fcm / fcl
         if(log(absc(w))+log(absc(fc(lf))) .lt. fplmin) go to 340
         if(mode.ge.3) go to 240
            if(absc(f-pq1) .lt. acch*absc(f) .and. pr)
     x                             write(iout,1020) lh,zlm+dell
      hpl = hcl * pq1
         if(absc(hpl).lt.fpmin.or.absc(hcl).lt.fpmin) go to 330
c
c *** idward recurrence from hcl,hpl(lf) (stored gc(l) is rl if reqd)
c *** renormalise fc,fcp at each lambda
c ***    zl   = zlm - min(id,0) here
c
  240 do 270 l = lf,l1,id
                     fcl = w* fc(l)
                      if(absc(fcl).lt.fpmin) go to 340
            if(ifcp) fpl = w*fcp(l)
                     fc(l)  = beta * fcl
            if(ifcp) fcp(l) = beta * (fpl - alpha * fcl) * cik
       if(mode .ge. 3) go to 260
       if(l.eq.lf)  go to 250
                      zl = zl + id
                      zid= zl * id
                      rl = gc(l)
         if(etane0)   then
                      sl = eta + zl*zl*xi
            if(mode.eq.1) then
              pl = gcp(l)
            else
              pl = zero
              if(absc(zl).gt.acch) pl = (sl*sl - rl*rl)/zid
            endif
           hcl1     = (sl*hcl - zid*hpl) / rl
           hpl      = (sl*hpl - pl *hcl) / rl
         else
           hcl1 = rl * hcl - hpl * id
           hpl  = (hcl - rl * hcl1) * id
         endif
         hcl      = hcl1
         if(absc(hcl).gt.fpmax) go to 320
  250    gc(l) = aa * (rk * hcl + df * fcl)
      if(mode.eq.1) gcp(l) = (aa *(rk*hpl +df*fpl) - alpha * gc(l)) *cik
         if(kfn.ge.3) aa = aa * q
  260    if(kfn.ge.3) beta = - beta * q
  270  last = min(last,(l1 - l)*id)
c
c *** come here after all soft errors to determine how many l values ok
c
  280  if(id.gt.0 .or.  last.eq.0) ifail = last
       if(id.lt.0 .and. last.ne.0) ifail = -3
c
c *** come here after all errors for this l range (zlm,zll)
c
  290 if(id.gt.0 .and. lf.ne.m1) go to 300
         if(ifail.lt.0) return
         if(rerr.gt.accb) write(iout,1070) rerr
         if(rerr.gt.0.1) ifail = -4
         return
c
c *** so on first block, 'f' started decreasing monotonically,
c                        or hit bound states for low zl.
c     thus redo m1 to lf-1 in reverse direction
c      i.e. do cf1a at zlmin & cf2 at zlm (midway between zlmin & zlmax)
c
  300 id = -1
      if(.not.unstab) lf = lf - 1
      donem = unstab
      lf = min(lf,l1)
      l1 = m1
      go to 10
c
c ***    error messages
c
  310 if(pr) write (iout,1000) xx
 1000 format(/' coulcc: cannot calculate irregular solutions for x =',
     x 1p,2d10.2,', as abs(x) is too small'/)
      return
  320 if(pr) write(iout,1010) zl+dell,'ir',hcl,'more',fpmax
 1010 format(' coulcc: at zl =',2f8.3,' ',a2,'regular solution (',1p,
     x 2e10.1,') will be ',a4,' than',e10.1)
      go to 280
  330 if(pr) write(iout,1010) zl+dell,'ir',hcl,'less',fpmin
      go to 280
  340 if(pr) write(iout,1010) zl+dell,'  ',fcl,'less',fpmin
      go to 280
  350 if(pr) write(iout,1010) zl+dell,'  ',fcl,'more',fpmax
      go to 280
 1020 format('0coulcc warning: linear independence between ''f'' and ''h
     x(',i1,')'' is lost at zl =',2f7.2,' (eg. coulomb eigenstate, or cf
     x1 unstable)'/)
  360 if(pr) write(iout,1030) zll+dell
 1030 format(' coulcc: (eta&l)/x too large for cf1a, and cf1 unstable at
     x l =',2f8.2)
      go to 280
  370 if(pr) write(iout,1040) z11,i
 1040 format(' coulcc: overflow in 1f1 series at zl =',2f8.3,' at term',
     x i5)
      go to 390
  380 if(pr) write(iout,1050) zlmin,zlm,zlm+one,zlmin+nl-one
 1050 format(' coulcc: both bound-state poles and f-instabilities occur'
     x ,', or multiple instabilities present.'
     x,/,' try calling twice,  first for zl from',2f8.3,' to',2f8.3,
     x ' (incl.)',/,20x,     'second for zl from',2f8.3,' to',2f8.3)
c     go to 390
  390 ifail = -1
      go to 290
 1060 format('0coulcc warning: as ''',a2,''' reflection rules not used,
     #errors can be up to',1p,d12.2/)
 1070 format('0coulcc warning: overall roundoff error approx.',1p,e11.1)
      end
