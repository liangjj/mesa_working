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
