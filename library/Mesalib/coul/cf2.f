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
