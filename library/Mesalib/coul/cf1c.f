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
