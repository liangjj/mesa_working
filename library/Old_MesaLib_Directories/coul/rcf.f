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
