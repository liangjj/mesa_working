*deck poch
      function poch (a, x)
c***begin prologue  poch
c***purpose  evaluate a generalization of pochhammer's symbol.
c***library   slatec (fnlib)
c***category  c1, c7a
c***type      single precision (poch-s, dpoch-d)
c***keywords  fnlib, pochhammer, special functions
c***author  fullerton, w., (lanl)
c***description
c
c evaluate a generalization of pochhammer's symbol
c (a)-sub-x = gamma(a+x)/gamma(a).  for x a non-negative integer,
c poch(a,x) is just pochhammer's symbol.  a and x are single precision.
c this is a preliminary version.  error handling when poch(a,x) is
c less than half precision is probably incorrect.  grossly incorrect
c arguments are not handled properly.
c
c***references  (none)
c***routines called  algams, alnrel, fac, gamma, gamr, r9lgmc, xermsg
c***revision history  (yymmdd)
c   770701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900727  added external statement.  (wrb)
c***end prologue  poch
      external gamma
      save pi
      data pi / 3.1415926535 89793238 e0 /
c***first executable statement  poch
      ax = a + x
      if (ax.gt.0.0) go to 30
      if (aint(ax).ne.ax) go to 30
c
      if (a .gt. 0.0 .or. aint(a) .ne. a) call xermsg ('slatec', 'poch',
     +   'a+x is non-positive integer but a is not', 2, 2)
c
c we know here that both a+x and a are non-positive integers.
c
      poch = 1.0
      if (x.eq.0.0) return
c
      n = x
      if (min(a+x,a).lt.(-20.0)) go to 20
c
      poch = (-1.0)**n * fac(-int(a))/fac(-int(a)-n)
      return
c
 20   poch = (-1.0)**n * exp ((a-0.5)*alnrel(x/(a-1.0))
     1  + x*log(-a+1.0-x) - x + r9lgmc(-a+1.) - r9lgmc(-a-x+1.) )
      return
c
c here we know a+x is not zero or a negative integer.
c
 30   poch = 0.0
      if (a.le.0.0 .and. aint(a).eq.a) return
c
      n = abs(x)
      if (real(n).ne.x .or. n.gt.20) go to 50
c
c x is a small non-positive integer, presummably a common case.
c
      poch = 1.0
      if (n.eq.0) return
      do 40 i=1,n
        poch = poch * (a+i-1)
 40   continue
      return
c
 50   absax = abs(a+x)
      absa = abs(a)
      if (max(absax,absa).gt.20.0) go to 60
      poch = gamma(a+x)*gamr(a)
      return
c
 60   if (abs(x).gt.0.5*absa) go to 70
c
c here abs(x) is small and both abs(a+x) and abs(a) are large.  thus,
c a+x and a must have the same sign.  for negative a, we use
c gamma(a+x)/gamma(a) = gamma(-a+1)/gamma(-a-x+1) *
c sin(pi*a)/sin(pi*(a+x))
c
      b = a
      if (b.lt.0.0) b = -a - x + 1.0
      poch = exp ((b-0.5)*alnrel(x/b) + x*log(b+x) - x +
     1  r9lgmc(b+x) - r9lgmc(b) )
      if (a.lt.0.0 .and. poch.ne.0.0) poch = poch/(cos(pi*x) +
     1  cot(pi*a)*sin(pi*x))
      return
c
 70   call algams (a+x, alngax, sgngax)
      call algams (a, alnga, sgnga)
      poch = sgngax * sgnga * exp(alngax-alnga)
c
      return
      end
