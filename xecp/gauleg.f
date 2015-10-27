      function gauleg(ntheta,nphi,pl,wt,f,xmphi,wphi,pi)
      implicit integer(a-z)
c
      real*8 gauleg
      real*8 pl(ntheta),wt(ntheta)
      real*8 f(ntheta,nphi)
      real*8 xmphi(nphi)
      real*8 wphi,pi
c
      real*8 a,b,zero,two
c
      parameter(zero=0.0d+00,two=2.0d+00)
c
      a=zero 
      do 20 theta=1,ntheta
         b=zero
         do 10 phi=1,nphi
            b=b+f(theta,phi)*xmphi(phi)
   10    continue
         a=a+b*pl(theta)*wt(theta)
   20 continue
c***** look into this below, it was originally there.
c     a=a*sqrt(two/pi)*wphi
c     should it be this in order to use this in both routines
c     can include normalization and this term in wphi
      a=a*wphi
c
      gauleg=a
c
c
      return
      end
