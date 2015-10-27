*deck @(#)bess.f	5.1  11/6/94
      function bess (z,l)
      implicit real*8(a-h,o-z)
c
c     evaluates modified spherical bessel function.
c
      common/dfac/dfac(30)
      common/fact/fac(17),fprod(9,9)

      parameter (zero=0.0d+00,one=1.0d+00,two=2.0d+00,five=5.0d+00)
      parameter (f16pt1=16.1d+00)
      parameter (fivem14=5.0d-14)
c
      if(z.lt.zero) then
         bess=zero
      else if(z.eq.0) then
         if(l.eq.0) then
            bess=one
         else
            bess=zero
         endif
      else if(z.le.five) then
         zp=z*z/two
         term=(z**l)/dfac(l+l+3)
         bess=term
         j=0
   10    j=j+1
            term=term*zp/float(j*(l+l+j+j+1))
            bess=bess+term
            if (abs(term/bess).gt.fivem14) go to 10
            bess=bess*exp(-z)
      else if(z.le.f16pt1) then
         rp=zero
         rm=zero
         tzp=two*z
         tzm=-tzp
         l1=l+1
         do 60 k1=1,l1
            k=k1-1
            rp=rp+fprod(k1,l1)/tzp**k
            rm=rm+fprod(k1,l1)/tzm**k
   60    continue
         bess=(rm-((-one)**l)*rp*exp(tzm))/tzp
      else
         rm=zero
         tzm=-two*z
         l1=l+1
         do 80 k1=1,l1
            k=k1-1
            rm=rm+fprod(k1,l1)/tzm**k
   80    continue
         bess=rm/(-tzm)
      endif
c
c
      return
      end
