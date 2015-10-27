*deck dpoch
      double precision function dpoch (a, x)
c***begin prologue  dpoch
c***purpose  evaluate a generalization of pochhammer's symbol.
c***library   slatec (fnlib)
c***category  c1, c7a
c***type      double precision (poch-s, dpoch-d)
c***keywords  fnlib, pochhammer, special functions
c***author  fullerton, w., (lanl)
c***description
c
c evaluate a double precision generalization of pochhammer's symbol
c (a)-sub-x = gamma(a+x)/gamma(a) for double precision a and x.
c for x a non-negative integer, poch(a,x) is just pochhammer's symbol.
c this is a preliminary version that does not handle wrong arguments
c properly and may not properly handle the case when the result is
c computed to less than half of double precision.
c
c***references  (none)
c***routines called  d9lgmc, dfac, dgamma, dgamr, dlgams, dlnrel, xermsg
c***revision history  (yymmdd)
c   770701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890911  removed unnecessary intrinsics.  (wrb)
c   890911  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900727  added external statement.  (wrb)
c***end prologue  dpoch
      double precision a, x, absa, absax, alnga, alngax, ax, b, pi,
     1  sgnga, sgngax, dfac, dlnrel, d9lgmc, dgamma, dgamr, dcot
      external dgamma
      save pi
      data pi / 3.1415926535 8979323846 2643383279 503 d0 /
c***first executable statement  dpoch
      ax = a + x
      if (ax.gt.0.0d0) go to 30
      if (aint(ax).ne.ax) go to 30
c
      if (a .gt. 0.0d0 .or. aint(a) .ne. a) call xermsg ('slatec',
     +   'dpoch', 'a+x is non-positive integer but a is not', 2, 2)
c
c we know here that both a+x and a are non-positive integers.
c
      dpoch = 1.0d0
      if (x.eq.0.d0) return
c
      n = x
      if (min(a+x,a).lt.(-20.0d0)) go to 20
c
      ia = a
      dpoch = (-1.0d0)**n * dfac(-ia)/dfac(-ia-n)
      return
c
 20   dpoch = (-1.0d0)**n * exp ((a-0.5d0)*dlnrel(x/(a-1.0d0))
     1  + x*log(-a+1.0d0-x) - x + d9lgmc(-a+1.0d0) - d9lgmc(-a-x+1.d0))
      return
c
c a+x is not zero or a negative integer.
c
 30   dpoch = 0.0d0
      if (a.le.0.0d0 .and. aint(a).eq.a) return
c
      n = abs(x)
      if (dble(n).ne.x .or. n.gt.20) go to 50
c
c x is a small non-positive integer, presummably a common case.
c
      dpoch = 1.0d0
      if (n.eq.0) return
      do 40 i=1,n
        dpoch = dpoch * (a+i-1)
 40   continue
      return
c
 50   absax = abs(a+x)
      absa = abs(a)
      if (max(absax,absa).gt.20.0d0) go to 60
      dpoch = dgamma(a+x) * dgamr(a)
      return
c
 60   if (abs(x).gt.0.5d0*absa) go to 70
c
c abs(x) is small and both abs(a+x) and abs(a) are large.  thus,
c a+x and a must have the same sign.  for negative a, we use
c gamma(a+x)/gamma(a) = gamma(-a+1)/gamma(-a-x+1) *
c sin(pi*a)/sin(pi*(a+x))
c
      b = a
      if (b.lt.0.0d0) b = -a - x + 1.0d0
      dpoch = exp ((b-0.5d0)*dlnrel(x/b) + x*log(b+x) - x
     1  + d9lgmc(b+x) - d9lgmc(b) )
      if (a.lt.0.0d0 .and. dpoch.ne.0.0d0) dpoch =
     1  dpoch/(cos(pi*x) + dcot(pi*a)*sin(pi*x) )
      return
c
 70   call dlgams (a+x, alngax, sgngax)
      call dlgams (a, alnga, sgnga)
      dpoch = sgngax * sgnga * exp(alngax-alnga)
c
      return
      end
