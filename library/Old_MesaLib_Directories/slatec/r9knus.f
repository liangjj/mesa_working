*deck r9knus
      subroutine r9knus (xnu, x, bknu, bknu1, iswtch)
c***begin prologue  r9knus
c***subsidiary
c***purpose  compute bessel functions exp(x)*k-sub-xnu(x) and exp(x)*
c            k-sub-xnu+1(x) for 0.0 .le. xnu .lt. 1.0.
c***library   slatec (fnlib)
c***category  c10b3
c***type      single precision (r9knus-s, d9knus-d)
c***keywords  bessel function, fnlib, special functions
c***author  fullerton, w., (lanl)
c***description
c
c compute bessel functions exp(x) * k-sub-xnu (x)  and
c exp(x) * k-sub-xnu+1 (x) for 0.0 .le. xnu .lt. 1.0 .
c
c series for c0k        on the interval  0.          to  2.50000d-01
c                                        with weighted error   1.60e-17
c                                         log weighted error  16.79
c                               significant figures required  15.99
c                                    decimal places required  17.40
c
c series for znu1       on the interval -7.00000d-01 to  0.
c                                        with weighted error   1.43e-17
c                                         log weighted error  16.85
c                               significant figures required  16.08
c                                    decimal places required  17.38
c
c***references  (none)
c***routines called  csevl, gamma, inits, r1mach, xermsg
c***revision history  (yymmdd)
c   770601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900720  routine changed from user-callable to subsidiary.  (wrb)
c   900727  added external statement.  (wrb)
c   920618  removed space from variable names.  (rwc, wrb)
c***end prologue  r9knus
      dimension alpha(15), beta(15), a(15), c0kcs(16), znu1cs(12)
      logical first
      external gamma
      save c0kcs, znu1cs, euler, sqpi2, aln2, ntc0k, ntznu1,
     1 xnusml, xsml, alnsml, alnbig, alneps, first
      data c0kcs( 1) /    .0601830572 42626108e0 /
      data c0kcs( 2) /   -.1536487143 3017286e0 /
      data c0kcs( 3) /   -.0117511760 08210492e0 /
      data c0kcs( 4) /   -.0008524878 88919795e0 /
      data c0kcs( 5) /   -.0000613298 38767496e0 /
      data c0kcs( 6) /   -.0000044052 28124551e0 /
      data c0kcs( 7) /   -.0000003163 12467283e0 /
      data c0kcs( 8) /   -.0000000227 10719382e0 /
      data c0kcs( 9) /   -.0000000016 30564460e0 /
      data c0kcs(10) /   -.0000000001 17069392e0 /
      data c0kcs(11) /   -.0000000000 08405206e0 /
      data c0kcs(12) /   -.0000000000 00603466e0 /
      data c0kcs(13) /   -.0000000000 00043326e0 /
      data c0kcs(14) /   -.0000000000 00003110e0 /
      data c0kcs(15) /   -.0000000000 00000223e0 /
      data c0kcs(16) /   -.0000000000 00000016e0 /
      data znu1cs( 1) /    .2033067569 9419173e0 /
      data znu1cs( 2) /    .1400779334 1321977e0 /
      data znu1cs( 3) /    .0079167969 61001613e0 /
      data znu1cs( 4) /    .0003398011 82532104e0 /
      data znu1cs( 5) /    .0000117419 75688989e0 /
      data znu1cs( 6) /    .0000003393 57570612e0 /
      data znu1cs( 7) /    .0000000084 25941769e0 /
      data znu1cs( 8) /    .0000000001 83336677e0 /
      data znu1cs( 9) /    .0000000000 03549698e0 /
      data znu1cs(10) /    .0000000000 00061903e0 /
      data znu1cs(11) /    .0000000000 00000981e0 /
      data znu1cs(12) /    .0000000000 00000014e0 /
      data euler / 0.5772156649 0153286e0 /
      data sqpi2 / 1.253314137 3155003e0 /
      data aln2 / 0.693147180 55994531e0 /
      data first /.true./
c***first executable statement  r9knus
      if (first) then
         ntc0k = inits (c0kcs, 16, 0.1*r1mach(3))
         ntznu1 = inits (znu1cs, 12, 0.1*r1mach(3))
c
         xnusml = sqrt (r1mach(3)/8.0)
         xsml = 0.1*r1mach(3)
         alnsml = log (r1mach(1))
         alnbig = log (r1mach(2))
         alneps = log (0.1*r1mach(3))
      endif
      first = .false.
c
      if (xnu .lt. 0. .or. xnu .ge. 1.0) call xermsg ('slatec',
     +   'r9knus', 'xnu must be ge 0 and lt 1', 1, 2)
      if (x .le. 0.) call xermsg ('slatec', 'r9knus', 'x must be gt 0',
     +   2, 2)
c
      iswtch = 0
      if (x.gt.2.0) go to 50
c
c x is small.  compute k-sub-xnu (x) and the derivative of k-sub-xnu (x)
c then find k-sub-xnu+1 (x).  xnu is reduced to the interval (-.5,+.5)
c then to (0., .5), because k of negative order (-nu) = k of positive
c order (+nu).
c
      v = xnu
      if (xnu.gt.0.5) v = 1.0 - xnu
c
c carefully find (x/2)**xnu and z**xnu where z = x*x/4.
      alnz = 2.0 * (log(x) - aln2)
c
      if (x.gt.xnu) go to 20
      if (-0.5*xnu*alnz-aln2-log(xnu) .gt. alnbig) call xermsg
     +   ('slatec', 'r9knus', 'x so small bessel k-sub-xnu overflows',
     +   3, 2)
c
 20   vlnz = v*alnz
      x2tov = exp (0.5*vlnz)
      ztov = 0.0
      if (vlnz.gt.alnsml) ztov = x2tov**2
c
      a0 = 0.5*gamma(1.0+v)
      b0 = 0.5*gamma(1.0-v)
      c0 = -euler
      if (ztov.gt.0.5 .and. v.gt.xnusml) c0 = -0.75 +
     1  csevl ((8.0*v)*v-1., c0kcs, ntc0k)
c
      if (ztov.le.0.5) alpha(1) = (a0-ztov*b0)/v
      if (ztov.gt.0.5) alpha(1) = c0 - alnz*(0.75 +
     1  csevl (vlnz/0.35+1.0, znu1cs, ntznu1))*b0
      beta(1) = -0.5*(a0+ztov*b0)
c
      z = 0.0
      if (x.gt.xsml) z = 0.25*x*x
      nterms = max (2.0, 11.0+(8.*alnz-25.19-alneps)/(4.28-alnz))
      do 30 i=2,nterms
        xi = i - 1
        a0 = a0/(xi*(xi-v))
        b0 = b0/(xi*(xi+v))
        alpha(i) = (alpha(i-1)+2.0*xi*a0)/(xi*(xi+v))
        beta(i) = (xi-0.5*v)*alpha(i) - ztov*b0
 30   continue
c
      bknu = alpha(nterms)
      bknud = beta(nterms)
      do 40 ii=2,nterms
        i = nterms + 1 - ii
        bknu = alpha(i) + bknu*z
        bknud = beta(i) + bknud*z
 40   continue
c
      expx = exp(x)
      bknu = expx*bknu/x2tov
c
      if (-0.5*(xnu+1.)*alnz-2.0*aln2.gt.alnbig) iswtch = 1
      if (iswtch.eq.1) return
      bknud = expx*bknud*2.0/(x2tov*x)
c
      if (xnu.le.0.5) bknu1 = v*bknu/x - bknud
      if (xnu.le.0.5) return
c
      bknu0 = bknu
      bknu = -v*bknu/x - bknud
      bknu1 = 2.0*xnu*bknu/x + bknu0
      return
c
c x is large.  find k-sub-xnu (x) and k-sub-xnu+1 (x) with y. l. luke-s
c rational expansion.
c
 50   sqrtx = sqrt(x)
      if (x.gt.1.0/xsml) go to 90
      an = -1.56 + 4.0/x
      bn = -0.29 - 0.22/x
      nterms = min (15, max1 (3.0, an+bn*alneps))
c
      do 80 inu=1,2
        xmu = 0.
        if (inu.eq.1 .and. xnu.gt.xnusml) xmu = (4.0*xnu)*xnu
        if (inu.eq.2) xmu = 4.0*(abs(xnu)+1.)**2
c
        a(1) = 1.0 - xmu
        a(2) = 9.0 - xmu
        a(3) = 25.0 - xmu
        if (a(2).eq.0.) result = sqpi2*(16.*x+xmu+7.)/(16.*x*sqrtx)
        if (a(2).eq.0.) go to 70
c
        alpha(1) = 1.0
        alpha(2) = (16.*x+a(2))/a(2)
        alpha(3) = ((768.*x+48.*a(3))*x + a(2)*a(3))/(a(2)*a(3))
c
        beta(1) = 1.0
        beta(2) = (16.*x+(xmu+7.))/a(2)
        beta(3) = ((768.*x+48.*(xmu+23.))*x + ((xmu+62.)*xmu+129.))
     1    / (a(2)*a(3))
c
        if (nterms.lt.4) go to 65
        do 60 i=4,nterms
          n = i - 1
          x2n = 2*n - 1
c
          a(i) = (x2n+2.)**2 - xmu
          qq = 16.*x2n/a(i)
          p1 = -x2n*(12*n*n-20*n-a(1))/((x2n-2.)*a(i)) - qq*x
          p2 = (12*n*n-28*n+8-a(1))/a(i) - qq*x
          p3 = -x2n*a(i-3)/((x2n-2.)*a(i))
c
          alpha(i) = -p1*alpha(i-1) - p2*alpha(i-2) - p3*alpha(i-3)
          beta(i) = -p1*beta(i-1) - p2*beta(i-2) - p3*beta(i-3)
 60     continue
c
 65     result = sqpi2*beta(nterms)/(sqrtx*alpha(nterms))
c
 70     if (inu.eq.1) bknu = result
        if (inu.eq.2) bknu1 = result
 80   continue
      return
c
 90   bknu = sqpi2/sqrtx
      bknu1 = bknu
      return
c
      end
