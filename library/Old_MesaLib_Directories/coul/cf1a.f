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
