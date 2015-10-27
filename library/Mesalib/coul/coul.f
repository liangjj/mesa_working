*deck cf1a
      function cf1a(rho,eta,xl,psi,eps,nmax,nused,fcl,re,fpmax,xx,g,c)
c
c     evaluate the asymptotic expansion for the
c            logarithmic derivative of the regular solution
c
c ***        cf1a  =  f   =  f'(xl,eta,rho)/f(xl,eta,rho)
c
c      that is valid for real(rho)>0, and best for rho >> eta**2, xl,
c      and is derived from the 2f0 expansions for h+ and h-
c      e.g. by froeberg (rev. mod. physics vol 27, p399 , 1955)
c      some lines of this subprogram are for convenience copied from
c           takemasa, tamura & wolter cpc 17 (1979) 351.
c
c     evaluate to accuracy eps with at most nmax terms.
c
c     if the terms start diverging,
c     the corresponding continued fraction is found by rcf
c     & evaluated progressively by steed's method to obtain convergence.
c
c      useful number also input:  fpmax = near-largest f.p. number
c
      implicit complex(a-h,o-z)
      dimension xx(2,nmax),g(nmax),c(nmax)
      real re,eps,t1,t2,t3,zero,one,two,at,atl,absc,fpmax
      data zero,one,two,ci / 0e+0, 1e+0, 2e+0, (0e+0,1e+0) /
      absc(w) = abs(real(w)) + abs(aimag(w))
c
      hpi = two*atan(one)
      t1 = sin(real(psi))
      t2 = cos(real(psi))
      atl= tanh(aimag(psi))
c             give cos(psi)/cosh(im(psi)), which always has correct sign
          cosl = cmplx( t2 , -t1 * atl )
      tanl = cmplx(t1,t2*atl) / cosl
      re = zero
      xll1= xl*(xl+one)
      etasq = eta*eta
      sl1=one
      sl=sl1
      sc1=zero
      sc=sc1
      tl1=sc
      tl=tl1
      tc1=one-eta/rho
      tc=tc1
      fcl  = tl + sl*tanl
      g(1) = (tc + sc*tanl) / fcl
      glast = g(1)
      atl = absc(glast)
         f    = glast
         d = one
         df   = glast
      j = 0
      do 10 n=2,nmax
      t1=n-1
      t2=two*t1-one
      t3=t1*(t1-one)
      denom=two*rho*t1
      c1=(eta*t2)/denom
      c2=(etasq+xll1-t3)/denom
      sl2=c1*sl1-c2*tl1
      tl2=c1*tl1+c2*sl1
      sc2=c1*sc1-c2*tc1-sl2/rho
      tc2=c1*tc1+c2*sc1-tl2/rho
      sl=sl+sl2
      tl=tl+tl2
      sc=sc+sc2
      tc=tc+tc2
      sl1=sl2
      tl1=tl2
      sc1=sc2
      tc1=tc2
      fcl  =  tl + sl*tanl
         if(absc(fcl).gt.fpmax .or. absc(fcl).lt.1./fpmax) go to 40
      gsum = (tc + sc*tanl) / fcl
      g(n) = gsum - glast
      glast = gsum
         at = absc(g(n))
         if(at.lt.absc(gsum)*eps) go to 20
      if(j.gt.0 .or. at.gt.atl .or. n.ge.nmax-2) j = j + 1
         if(j.eq.0) go to 10
            call rcf(g,c,j,n,xx,eps)
              if(n.lt.0) go to 40
            do 60 k=max(j,2),n
               d = one/(d*c(k) + one)
               df = df*(d - one)
               f = f + df
         if(absc(df) .lt. absc(f)*eps) go to 30
         if(df.eq.zero.and.f.eq.zero.and.n.ge.4) go to 30
   60         continue
         j = n
   10    atl = at
      k = -nmax
      go to 30
   20 fcl = fcl * cosl
         cf1a = gsum
         re = at / absc(gsum)
         nused = n
         return
   30 cf1a = f
      fcl = fcl * cosl
         re = absc(df) / absc(f)
         nused = k
      return
   40 cf1a = g(1)
      fcl = 1.0
      re = 1.0
      nused = 0
      return
      end
*deck cf1c
      function cf1c(x,eta,zl,eps,fcl,tpk1,etane0,limit,err,nfp,
     x              acch,fpmin,fpmax,pr,caller)
      implicit complex(a-h,o-z)
      common /io/ inp, iout
      logical pr,etane0
      real one,two,eps,err,acch,fpmin,fpmax,absc,small,rk,px
      character*6 caller
      data one,two / 1e+0, 2e+0 /
      absc(w) = abs(real(w)) + abs(aimag(w))
c
c
c ***    evaluate cf1  =  f   =  f'(zl,eta,x)/f(zl,eta,x)
c
c        using complex arithmetic
c
      fcl = one
      xi = one/x
      pk  = zl + one
      px  = pk  + limit
   10 ek  = eta / pk
        rk2 =          one + ek*ek
      f   = (ek + pk*xi)*fcl + (fcl - one)*xi
      pk1 =  pk + one
         tpk1 = pk + pk1
      tk  = tpk1*(xi + ek/pk1)
      if(etane0) then
c ***   test ensures b1 .ne. zero for negative eta etc.; fixup is exact.
             if(absc(tk) .gt. acch)  go to 20
             fcl  = rk2/(one + (eta/pk1)**2)
             sl   = tpk1*xi * (tpk1+two)*xi
             pk   =  two + pk
             go to 10
         endif
   20 d   =  one/tk
      df  = -fcl*rk2*d
            if(real(pk).gt.real(zl)+two) fcl = - rk2 * sl
            fcl = fcl * d * tpk1 * xi
      f   =  f  + df
c
c ***   begin cf1 loop on pk = k = lambda + 1
c
      rk    = one
      small    = sqrt(fpmin)
   30 pk    = pk1
        pk1 = pk1 + one
         tpk1 = pk + pk1
         if(etane0) then
           ek  = eta / pk
           rk2 =          one + ek*ek
          endif
        tk  = tpk1*(xi + ek/pk1)
        d   =  tk - d*rk2
              if(absc(d) .gt. acch)             go to 40
              if(pr) write (iout,1000) caller,d,df,acch,pk,ek,eta,x
              rk= rk +   one
              if( rk .gt. two )                  go to 50
   40 d     = one/d
            fcl = fcl * d * tpk1*xi
            if(absc(fcl).lt.small) fcl = fcl / small
            if(absc(fcl).gt.fpmax) fcl = fcl / fpmax
        df  = df*(d*tk - one)
        f   = f  + df
              if( real(pk) .gt. px ) go to 50
      if(absc(df) .ge. absc(f)*eps)             go to 30
                nfp = pk - zl - 1
                  err = eps * sqrt(real(nfp))
      cf1c = f
      return
 1000 format(/' ',a6,': cf1 accuracy loss: d,df,acch,k,eta/k,eta,x = ',
     x    /1x,1p,13d9.2/)
   50 if(pr) write (iout,1010) caller,limit,abs(x)
 1010 format(' ',a6,': cf1 has failed to converge after ',i10  ,' iterat
     xions as abs(x) =',f15.0)
      err = two
      return
      end
*deck cf1r
      function cf1r(x,eta,zl,eps,fcl,tpk1,etane0,limit,err,nfp,
     x              acch,fpmin,fpmax,pr,caller)
      common /io/ inp, iout
      logical pr,etane0
      character*6 caller
      data one,two / 1e+0, 2e+0 /
c
c
c ***    evaluate cf1  =  f   =  f'(zl,eta,x)/f(zl,eta,x)
c
c        using real arithmetic
c
      fcl = one
      xi = one/x
      pk  = zl + one
      px  = pk  + limit
   10 ek  = eta / pk
        rk2 =          one + ek*ek
      f   = (ek + pk*xi)*fcl + (fcl - one)*xi
      pk1 =  pk + one
         tpk1 = pk + pk1
      tk  = tpk1*(xi + ek/pk1)
      if(etane0) then
c ***   test ensures b1 .ne. zero for negative eta etc.; fixup is exact.
             if(abs(tk) .gt. acch)  go to 20
             fcl  = rk2/(one + (eta/pk1)**2)
             sl   = tpk1*xi * (tpk1+two)*xi
             pk   =  two + pk
             go to 10
         endif
   20 d   =  one/tk
      df  = -fcl*rk2*d
            if(pk.gt.zl+two) fcl = - rk2 * sl
            fcl = fcl * d * tpk1 * xi
      f   =  f  + df
c
c ***   begin cf1 loop on pk = k = lambda + 1
c
      rk    = one
      small    = sqrt(fpmin)
   30 pk    = pk1
        pk1 = pk1 + one
         tpk1 = pk + pk1
         if(etane0) then
           ek  = eta / pk
           rk2 =          one + ek*ek
          endif
        tk  = tpk1*(xi + ek/pk1)
        d   =  tk - d*rk2
              if(abs(d) .gt. acch)             go to 40
              if(pr) write (iout,1000) caller,d,df,acch,pk,ek,eta,x
              rk= rk +   one
              if( rk .gt. two )                  go to 50
   40 d     = one/d
            fcl = fcl * d * tpk1*xi
            if(abs(fcl).lt.small) fcl = fcl / small
            if(abs(fcl).gt.fpmax) fcl = fcl / fpmax
        df  = df*(d*tk - one)
        f   = f  + df
              if( pk .gt. px ) go to 50
      if(abs(df) .ge. abs(f)*eps)             go to 30
                nfp = pk - zl - 1
                  err = eps * sqrt(real(nfp))
      cf1r = f
      return
 1000 format(/' ',a6,': cf1 accuracy loss: d,df,acch,k,eta/k,eta,x = ',
     x    /1x,1p,7d9.2/)
   50 if(pr) write (iout,1010) caller,limit,abs(x)
 1010 format(' ',a6,': cf1 has failed to converge after ',i10  ,' iterat
     xions as abs(x) =',f15.0)
      err = two
      return
      end
*deck cf2
      function cf2(x,eta,zl,pm,eps,limit,err,npq,acc8,acch,
     x             pr,accur,dell,caller)
      implicit complex(a-h,o-z)
      common /io/ inp, iout
      logical pr
      real eps,err,acc8,acch,accur,ta,rk,
     x       absc,zero,half,one,two
      character*6 caller
      data zero,half,one,two / 0e+0, .5e+0, 1e+0, 2e+0 /
      absc(w) = abs(real(w)) + abs(aimag(w))
c
c                                    (omega)        (omega)
c *** evaluate  cf2  = p + pm.q  =  h   (eta,x)' / h   (eta,x)
c                                    zl             zl
c     where pm = omega.i
c
      ta = two*limit
      e2mm1 = eta*eta + zl*zl + zl
      etap = eta * pm
      xi = one/x
      wi = two*etap
      rk = zero
      pq = (one - eta*xi) * pm
      aa = -e2mm1 + etap
      bb = two*(x - eta + pm)
         rl = xi * pm
      if(absc(bb).lt.acch) then
         rl = rl * aa / (aa + rk + wi)
         pq = pq + rl * (bb + two*pm)
            aa = aa + two*(rk+one+wi)
            bb = bb + (two+two)*pm
            rk = rk + (two+two)
         endif
      dd = one/bb
      dl = aa*dd* rl
   10 pq    = pq + dl
         rk = rk + two
         aa = aa + rk + wi
         bb = bb + two*pm
         dd = one/(aa*dd + bb)
         dl = dl*(bb*dd - one)
            err = absc(dl)/absc(pq)
         if(err.ge.max(eps,acc8*rk*half) .and. rk.le.ta) go to 10
c
         npq   = rk/two
         pq    = pq + dl
           if(pr.and.npq.ge.limit-1 .and. err.gt.accur)
     x            write(iout,1000) caller,int(aimag(pm)),npq,err,zl+dell
 1000 format(' ',a6,': cf2(',i2,') not converged fully in ',i7,
     x' iterations, so error in irregular solution =',1p,d11.2,' at zl
     x=', 0p,2f8.3)
      cf2 = pq
      return
      end
*deck clogam
      function clogam(z)
      common /io/ inp, iout
c
c     this routine computes the logarithm of the gamma function gamma(z)
c     for any complex argument 'z' to any accuracy preset by call logam
c
      complex z,u,v,h,r,clogam,cdigam,ser
      dimension b(15),bn(15),bd(15)
c
      data lerr /6/, nx0 /6/, nb /15/,
     x  zero,one,two,four,half,quart /0e+0,1e+0,2e+0,4e+0,.5e+0,.25e+0/
      data bn(1),bd(1)    / +1e+0,   6e+0 /,
     x     bn(2),bd(2)    / -1e+0,  30e+0 /,
     x     bn(3),bd(3)    / +1e+0,  42e+0 /,
     x     bn(4),bd(4)    / -1e+0,  30e+0 /,
     x     bn(5),bd(5)    / +5e+0,  66e+0 /,
     x     bn(6),bd(6)    /          -691e+0,  2730e+0/,
     x     bn(7),bd(7)    /          +  7e+0,     6e+0/,
     x     bn(8),bd(8)    /         -3617e+0,   510e+0/,
     x     bn(9),bd(9)    /         43867e+0,   798e+0/,
     x     bn(10),bd(10)  /       -174611e+0,   330e+0/,
     x     bn(11),bd(11)  /        854513e+0,   138e+0/,
     x     bn(12),bd(12)  /    -236364091e+0,  2730e+0/,
     x     bn(13),bd(13)  /     + 8553103e+0,     6e+0/,
     x     bn(14),bd(14)  /  -23749461029e+0,   870e+0/,
     x     bn(15),bd(15)  / 8615841276005e+0, 14322e+0/
      data fplmin / -850e+0 /
c
      x=real(z)
      t=aimag(z)
      mx = int(real(accur*100 - x))
      if(abs(abs(x)-mx) + abs(t).lt.accur*100) go to 60
      f=abs(t)
      v=cmplx(x,f)
      if(x .lt. zero) v=one-v
      h=zero
      c=real(v)
      n=nx0-int(c)
      if(n .lt. 0) go to 30
      h=v
      d=aimag(v)
      a=atan2(d,c)
      if(n .eq. 0) go to 20
      do 10 i = 1,n
      c=c+one
      v=cmplx(c,d)
      h=h*v
   10 a=a+atan2(d,c)
   20 h=cmplx(half*log(real(h)**2+aimag(h)**2),a)
      v=v+one
   30 r=one/v**2
      ser = b(nt)
      do 40 j=2,nt
        k = nt+1 - j
   40 ser = b(k) + r*ser
      clogam = hl2p+(v-half)*log(v)-v + ser/v - h
      if(x .ge. zero) go to 50
c
      a= int(x)-one
      c=pi*(x-a)
      d=pi*f
c     e=exp(-two*d)
        e = zero
        f = -two*d
        if(f.gt.fplmin) e = exp(f)
      f=sin(c)
      e= d + half*log(e*f**2+quart*(one-e)**2)
      f=atan2(cos(c)*tanh(d),f)-a*pi
      clogam=alpi-cmplx(e,f)-clogam
c
   50 if(sign(one,t) .lt. -half) clogam=conjg(clogam)
      return
c
   60 write(iout,1000) 'clogam',x
 1000 format(1x,a6,' ... argument is non positive integer = ',f20.2)
      clogam = zero
      return
c
      entry cdigam(z)
c
c     this routine computes the logarithmic derivative of the gamma
c     function  psi(z) = digamma(z) = d (ln gamma(z))/dz  for any
c     complex argument z, to any accuracy preset by call logam(acc)
c
      u=z
      x=real(u)
      a=abs(x)
      if(abs(aimag(u)) + abs(a + int(x)) .lt. accur) go to 110
      if(x .lt. zero) u=-u
      v=u
      h=zero
      n=nx0-int(a)
      if(n .lt. 0) go to 90
      h=one/v
      if(n .eq. 0) go to 80
      do 70 i = 1,n
      v=v+one
   70 h=h+one/v
   80 v=v+one
   90 r=one/v**2
      ser = b(nt) * (2*nt-1)
      do 100 j=2,nt
        k = nt+1 - j
  100 ser = b(k)*(2*k-1) + r*ser
      cdigam = log(v) - half/v - r*ser - h
      if(x .ge. zero) return
      h=pi*u
      cdigam = cdigam + one/u + pi*cos(h)/sin(h)
      return
c
  110 write(iout,1000) 'cdigam',x
      cdigam=zero
      return
c
      entry logam(acc)
c
c      initialisation call for calculations to accuracy 'acc'
c
      nx0 = 6
      x0  = nx0 + one
      pi = four*atan(one)
      alpi = log(pi)
      hl2p = log(two*pi) * half
      accur = acc
      do 120 k=1,nb
       f21 = k*2 - one
       b(k) = bn(k) / (bd(k) * k*two * f21)
       err = abs(b(k)) * k*two / x0**f21
  120 if(err.lt.acc) go to 130
       nx0 = int((err/acc)**(one/f21) * x0)
       k = nb
  130 nt = k
c     print *,' logam requires k = ',k ,' with cutoff at x =',nx0+1
      return
      end
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
     x     fpmax,fpmin,fplmax,fplmin / 1e+2000,1e-2000,850e+0,-850e+0/,
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
*deck f11
      function f11(x,eta,zl,p,eps,limit,kind,err,nits,fpmax,acc8,acc16)
      common /io/ inp, iout
      complex x,eta,zl,p,aa,bb,z,f11,cdigam,ci
       complex dd,g,f,ai,bi,t
      logical zllin
      real*16 ar,br,gr,gi,dr,di,tr,ti,ur,ui,fi,fi1,den
      data zero,one,two / 0e+0, 1e+0, 2e+0 /, ci / (0e+0, 1e+0) /
      absc(aa) = abs(real(aa)) + abs(aimag(aa))
      nintc(aa) = nint(real(aa))
c
c *** evaluate the hypergeometric function 1f1
c                                        i
c            f (aa;bb; z) = sum  (aa)   z / ( (bb)  i| )
c           1 1              i       i            i
c
c     to accuracy eps with at most limit terms.
c  if kind = 0 : using extended precision but real arithmetic only,
c            1 : using normal precision in complex arithmetic,
c   or       2 : using normal complex arithmetic, but with cdigam factor
c
c  where
         aa = zl+one - eta*p
         bb = two*(zl+one)
c  and
         z  = two*p*x
c
         zllin = real(bb).le.zero .and. abs(bb-nintc(bb)).lt.acc8**0.25
             if(.not.zllin.or.real(bb)+limit.lt.1.5) go to 10
                nits = -1
                return
   10 if(limit.le.0) then
         f11 = zero
         err = zero
         nits= 1
         return
         endif
      ta = one
      rk = one
      if(kind.le.0.and.absc(z)*absc(aa).gt.absc(bb) * 1.0) then
         dr = one
         di = zero
         gr = one
         gi = zero
         ar = real(aa)
         br = real(bb)
         fi = zero
      do 20 i=2,limit
         fi1 = fi + one
         tr = br * fi1
         ti = aimag(bb) * fi1
         den= one / (tr*tr + ti*ti)
         ur = (ar*tr + aimag(aa)*ti) * den
         ui = (aimag(aa)*tr - ar*ti) * den
         tr = ur*gr - ui*gi
         ti = ur*gi + ui*gr
         gr = real(z) * tr - aimag(z)*ti
         gi = real(z) * ti + aimag(z)*tr
         dr = dr + gr
         di = di + gi
            err = abs(gr) + abs(gi)
               if(err.gt.fpmax) go to 60
            rk  = abs(dr) + abs(di)
            ta = max(ta,rk)
         if(err.lt.rk*eps .or. i.ge.4.and.err.lt.acc16) go to 30
         fi = fi1
         ar = ar + one
   20    br = br + one
c
   30    f11 = dr + ci * di
         err = acc16 * ta / rk
c
      else
c* ---------------------------------- alternative code
c*    if real*16 arithmetic is not available, (or already using it|),
c*    then use kind > 0
         g = one
          f = one
          if(kind.ge.2) f = cdigam(aa) - cdigam(bb) - cdigam(g)
         dd = f
         do 40 i=2,limit
            ai = aa + (i-2)
            bi = bb + (i-2)
            r  = i-one
         g = g * z * ai / (bi * r)
         if(kind.ge.2)
c                              multiply by (psi(a+r)-psi(b+r)-psi(1+r))
     x        f = f + one/ai - one/bi - one/r
         t  = g * f
         dd = dd + t
            err = absc(t)
               if(err.gt.fpmax) go to 60
            rk = absc(dd)
         ta = max(ta,rk)
         if(err.lt.rk*eps.or.err.lt.acc8.and.i.ge.4) go to 50
   40    continue
 
   50    err = acc8 * ta / rk
         f11 = dd
c* ------------------------------------------- end of alternative code
      endif
   60    nits = i
      return
      end
*deck f20
      function f20(aa,bb,z,eps,jmax,re,fpmax,n,x)
c
c     evaluate the hypergeometric function 2f0
c                                             i
c            f (aa,bb;;z) = sum  (aa)  (bb)  z / i
c           2 0              i       i     i
c
c     to accuracy eps with at most jmax terms.
c
c     if the terms start diverging,
c     the corresponding continued fraction is found by rcf
c     & evaluated progressively by steed's method to obtain convergence.
c
c      useful number also input:  fpmax = near-largest f.p. number
c
      implicit complex(a-h,o-z)
      dimension x(jmax,4)
      logical finite
      real ep,eps,at,atl,absc,re,fpmax
      data one,zero / (1e+0,0e+0), (0e+0,0e+0) /
      absc(w) = abs(real(w)) + abs(aimag(w))
      nintc(w) = nint(real(w))
c
      re = 0.0
      x(1,1) = one
      sum = x(1,1)
      atl = absc(x(1,1))
         f    = sum
         d = one
         df   = sum
      j = 0
      ep = eps * jmax *10.
      ma = - nintc(aa)
      mb = - nintc(bb)
      finite = abs(abs(real(aa))-ma).lt.ep .and. abs(aimag(aa)).lt.ep
     x    .or. abs(abs(real(bb))-mb).lt.ep .and. abs(aimag(bb)).lt.ep
      imax = jmax
      if(finite.and.ma.ge.0) imax = min(ma+1,imax)
      if(finite.and.mb.ge.0) imax = min(mb+1,imax)
      do 10 i=2,imax
      x(i,1) = x(i-1,1) * z * (aa+i-2) * (bb+i-2) / (i-1)
         if(absc(x(i,1)).gt.fpmax) go to 40
      at = absc(x(i,1))
         if(j.eq.0) then
                 sum = sum + x(i,1)
                 if(at .lt. absc(sum)*eps) go to 20
               endif
      if(finite) go to 10
      if(j.gt.0 .or. at.gt.atl .or. i.ge.jmax-2) j = j + 1
         if(j.eq.0) go to 10
         call rcf(x(1,1),x(1,2),j,i,x(1,3),eps)
              if(i.lt.0) go to 40
            do 50 k=max(j,2),i
            d = one/(d*x(k,2) + one)
            df = df*(d - one)
            f = f + df
            if(absc(df) .lt. absc(f)*eps) go to 30
            if(df.eq.zero.and.f.eq.zero.and.i.ge.4) go to 30
   50       continue
         j = i
   10 atl = at
      if(.not.finite) i = -jmax
   20 n = i
       f20 = sum
       if(.not.finite) re  = at / absc(sum)
       return
   30 f20 = f
      re = absc(df) / absc(f)
      n = k
      return
   40 i = 0
      go to 20
      end
c----------------------------------------------------------------------c
c                 coulomb function code of barnett from                c
c                             cpc library                              c
c----------------------------------------------------------------------c
*deck grncal
      subroutine grncal (rn,ln,en,charge,gr1,gr1p,gr2,gr2p)
c
      parameter (nc8=12, np8=61, ns8=4, nr8=10, no8=10, nx8=78, mx8=10)
c
      dimension fc(np8), gc(np8), fcp(np8), gcp(np8), sig(np8)
c
      complex cark, xx, eta, zlmin, fc, gc, fcp, gcp, sig, civv, sc, sc1
      complex fc1, gc1, civ, cn1, cn11, ws, cexp2, fc1p, gc1p
      data stw /1.4142135e+0/
      data zero, one, pi2 /0.e0,1.e0,1.5707963e0/
c
c this subroutine constructs the regular and irregular modulating
c functions(gr1 and gr2) and their respective derivatives(gr1p,gr2p)
c as a function of r, the radial variable, and the channel
c orbital angular momentum, l, from the products of the general
c coulomb package , coulcc , of thompson and barnett (cpc 36,
c 363(1985).
c
c these modulating functions solve the following second-order,
c ordinary differential equation:
c
c ( d2/d2x + ( 1. - (2.*eta)/x - (l*(l+1))/x*x)) gi(x,eta,l) = 0
c
c where x and eta can be complex.
c
c description of parameters :
c
c* calling variables
c
c rn = radius in a0
c ln = orbital angular momentum
c en = energy in rydbergs(k2)
c charge = residual tartget charge = z
c gr1,gr1p - regular function and derivative
c gr2,gr2p - irregular function and derivative
c
c** conventions
c
c* 1) z = 0 , k2 .ge. 0
c
c     x = k*r     eta = 0
c
c     g1 = ricatti-bessel function of order l = x*jl(x)
c     g2 = ricatti-neumann function of order l = x*nl(x)
c
c     jl(x) = spherical bessel function of order l
c     nl(x) = spherical neumann function of order l
c
c coulcc: jl(x) , nl(x) , jl(x)' , nl(x)'
c              where ' indicates a derivative wrt x (not r)
c
c asymptotic forms
c
c      g1 = sin(x - l*pi/2)
c      g2 = cos(x - l*pi/2)
c
c* 2) z = 0 , k2 .lt. 0
c
c    x = i*k*r    eta = 0    i = sqrt(-1)
c
c    g1 = sqrt(2)*((-i)**l)*x*jl(x)
c    g2 = (-(i)**l)*x*hl+(x)/sqrt(2)
c
c    jl(x) = spherical bessel function of order l (x imaginary)
c    hl+(x) = jl(x) + i*nl(x)
c    nl(x) = spherical neumann function of order l
c
c       jl and hl+ are either purely real or purely imaginary quantities
c      the scaling factors proportional to (i)**l are used to
c      obtain a real quantitiy
c
c coulcc: jl(x) , hl+(x) , jl(x)' , hl+(x)'
c
c asymptotic form
c
c     g1 = exp(+k*r)/sqrt(2)
c     g2 = exp(-k*r)/sqrt(2)
c
c* 3) z .ne. 0 , k2 .ge. 0
c
c     x = k*r    eta = -z/k
c
c     g1 = regular coulomb function of order l = fl(x)
c     g2 = irregular coulomb function of order l = gl(x)
c
c coulc c: fl(x) , gl(x) , fl(x)' , gl(x)' , sig(l)
c
c asymptotic form
c    g1 = sin(thetal)
c    g2 = cos(thetal)
c
c    thetal = x - eta*alog(2*x) - pi*l/2 + sig(l)
c    sig(l) = coulomb phase shift
c
c* 4) z .ne. 0 , k2 .lt. 0
c
c     x = i*k*r    eta = i*z/k
c
c     g1 = regular whitaker function = (-i*sqrt(2)/ws)*fl(x)
c     g2 = irregualr whitaker function = ws*hl+(x)/sqrt(2)
c
c       ws = exp(i(pi/2(l + i*eta) - sig(l)))
c       fl(x) = regaular coulomb function of imaginary argument
c       gl(x) = irregular coulomb function
c       hl+(x) = fl(x) + i*gl(x)
c
c coulcc: fl(x) , hl+(x) , fl(x)' , hl+(x)' , sig(l)
c
c asymptotic forms
c
c     g1 = exp(+kr)*((2*kr)**etaa)/sqrt((2)
c     g2 = -exp(-kr)*((2*kr)**-etaa)/sqrt((2)
c
c          etaa = -z/k
c
c***********************************************************************
c
      civ=cmplx(zero,one)
      civv=cmplx(zero,-one)
      if (charge.ne.0.) go to 10
      eta=cmplx(zero,zero)
      l=ln
      l1=l+1
      xl=float(l)
      zlmin=cmplx(xl,0.)
      e=en
      xk=sqrt(abs(e))
      if (e.lt.zero) then
      sc=cmplx(one,zero)
      if (l.ne.0) sc=civv**l
      sc1=cmplx(-one,zero)
      if (l.ne.0) sc1=-(civ**l)
      endif
      r=rn
      ark=xk*r
c
      if (e.ge.zero) then
      xx=cmplx(ark,zero)
      call coulcc (xx,eta,zlmin,1,fc,gc,fcp,gcp,sig,1,1,ifail)
      gr1=ark*real(fc(1))
      gr2=-ark*real(gc(1))
      gr1p=xk*(ark*real(fcp(1))+real(fc(1)))
      gr2p=-xk*(ark*real(gcp(1))+real(gc(1)))
      else
      xx=cmplx(zero,ark)
      call coulcc (xx,eta,zlmin,1,fc,gc,fcp,gcp,sig,11,1,ifail)
      fc1=sc*fc(1)
      gc1=sc1*gc(1)
      gr1=ark*real(fc1)*stw
      gr2=ark*real(gc1)/stw
      fc1p=fc(1)+xx*fcp(1)
      fc1p=sc*fc1p
      gc1p=gc(1)+xx*gcp(1)
      gc1p=sc1*gc1p
      gr1p=xk*real(fc1p)*stw
      gr2p=xk*real(gc1p)/stw
      endif
      go to 20
   10 l=ln
      zlmin=l
      e=en
      xk=sqrt(abs(e))
      etaa=charge/xk
      if (e.ge.zero) then
      eta=cmplx(-etaa,zero)
      else
      eta=cmplx(zero,etaa)
      cn11=civ*eta
      cn11=pi2*(zlmin+cn11)
      endif
      r=rn
      ark=xk*r
      if (e.ge.zero) then
      xx=cmplx(ark,zero)
      call coulcc (xx,eta,zlmin,1,fc,gc,fcp,gcp,sig,1,0,ifail)
      gr1=real(fc(1))+aimag(fc(1))
      gr2=real(gc(1))+aimag(gc(1))
      gr1p=xk*(real(fcp(1))+aimag(fcp(1)))
      gr2p=xk*(real(gcp(1))+aimag(gcp(1)))
      else
      xx=cmplx(zero,ark)
      call coulcc (xx,eta,zlmin,1,fc,gc,fcp,gcp,sig,11,0,ifail)
      cn1=cn11-sig(1)
      cn1=civ*cn1
      ws=cexp(cn1)
      gc1=ws*gc(1)
      gc1p=ws*civ*gcp(1)
      gr2=real(gc1)/stw
      gr2p=xk*real(gc1p)/stw
      cexp2=-civ/ws
      fc1=cexp2*fc(1)
      fc1p=civ*cexp2*fcp(1)
      gr1=real(fc1)*stw
      gr1p=xk*real(fc1p)*stw
      endif
c
   20 return
      end
*deck rcf
      subroutine rcf(a,b,ibeg,inum,xx,eps)
c
c*******************************************************************
c
c  rcf converts polynomial a to the corresponding continued
c         fraction, in 'normal'  form with coefficients b
c         by the 'p algorithmn' of patry & gupta
c
c   a(z) = a1/z + a2/z**3 + a3/z**5 + ... + an/z**(2n-1)
c
c   b(z) = b1/z+ b2/z+ b3/z+ .../(z+ bn/z)
c
c  data:
c   a     vector a(k), k=1,inum         input
c   b     vector b(k), k=ibeg,inum      output
c   ibeg  order of first coef. calc.    input
c   inum  order of a, even or odd       input
c   xx    auxiliary vector of length .ge. length of vector b
c         caller provides space for a,b,xx
c     note that neither of the first two terms a(1) a(2) should be zero
c             & the user can start the calculation with any value of
c                ibeg provided the c.f. coefs have been already
c                calculated up to inum = ibeg-1
c             & the method breaks down as soon as the absolute value
c                of a c.f. coef. is less than eps.    at the time of the
c                break up xx(1) has been replaced by 1e-50, and inum has
c                been replaced by minus times the number of this coef.
c   algorithm: j.patry & s.gupta,
c              eir-bericht nr. 247,
c              eidg. institut fur reaktorforschung wuerenlingen
c              wueringlingen, schweiz.
c              november 1973
c   see also:  haenggi,roesel & trautmann,
c              jnl. computational physics, vol 137, pp242-258 (1980)
c   note:      restart procedure modified by i.j.thompson
c
c*******************************************************************
c
      implicit complex(a-h,o-z)
      dimension a(100),b(100),xx(2,100)
      logical even
      real eps
      common /io/ inp, iout
      common /rcfcm2/ x1,m2m1,mp12,even,m
c     ibn = ibeg + inum - 1
      ibn = inum
c                             b(ibn) is last value set on this call
      if(ibeg.gt.4 .and. m .ne. ibeg-1) go to 90
c                             b(m) is last value set in previous call
      if(ibeg.gt.4) go to 50
      if(ibeg.eq.4) go to 20
      b(1) = a(1)
      if(ibn.ge.2) b(2) = - a(2)/a(1)
      if(ibn.lt.3) go to 10
      x0 = a(3) / a(2)
      xx(2,1) = b(2)
      xx(1,1) = - x0
      xx(1,2) = 0.
      b(3) = -x0 - b(2)
      x0 = -b(3) * a(2)
      m = 3
      mp12 = 2
      even = .true.
      if(ibn.gt.3) go to 20
   10 return
   20 if(abs(b(3)) .lt. eps*abs(x0)) goto 80
      m = 4
   30 x1 = a(m)
      m2m1 = mp12
      mp12 = m2m1 + 1
      if(even) mp12 = m2m1
      do 40 k=2,mp12
   40 x1 = x1 + a(m-k+1) * xx(1,k-1)
      b(m) = - x1/x0
      if(m.ge.ibn) return
   50 if(abs(b(m)).lt.eps*abs(x0)) go to 80
      k = m2m1
   60 xx(2,k) = xx(1,k) + b(m) * xx(2,k-1)
      k = k-1
      if(k.gt.1) go to 60
      xx(2,1) = xx(1,1) + b(m)
      do 70 k=1,m2m1
      x0 = xx(2,k)
      xx(2,k) = xx(1,k)
   70 xx(1,k) = x0
      x0 = x1
      xx(1,m2m1+1) = 0.
      m = m+1
      even = .not.even
      go to 30
   80 inum = -m
c     xx(1,1) = 1.e-50
c     print 1000,m
c1000 format('0rcf: zero cf coefficient at position ',i4/)
      return
   90 print 1000,m,ibeg-1
 1000 format('0rcf: last call set m =',i4,', but restart requires',i4)
      stop
      end
