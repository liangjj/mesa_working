*deck r9gmic
      function r9gmic (a, x, alx)
c***begin prologue  r9gmic
c***subsidiary
c***purpose  compute the complementary incomplete gamma function for a
c            near a negative integer and for small x.
c***library   slatec (fnlib)
c***category  c7e
c***type      single precision (r9gmic-s, d9gmic-d)
c***keywords  complementary incomplete gamma function, fnlib, small x,
c             special functions
c***author  fullerton, w., (lanl)
c***description
c
c compute the complementary incomplete gamma function for a near
c a negative integer and for small x.
c
c***references  (none)
c***routines called  alngam, r1mach, xermsg
c***revision history  (yymmdd)
c   770701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900720  routine changed from user-callable to subsidiary.  (wrb)
c***end prologue  r9gmic
      save euler, eps, bot
      data euler / .5772156649 015329 e0 /
      data eps, bot / 2*0.0 /
c***first executable statement  r9gmic
      if (eps.eq.0.0) eps = 0.5*r1mach(3)
      if (bot.eq.0.0) bot = log(r1mach(1))
c
      if (a .gt. 0.0) call xermsg ('slatec', 'r9gmic',
     +   'a must be near a negative integer', 2, 2)
      if (x .le. 0.0) call xermsg ('slatec', 'r9gmic',
     +   'x must be gt zero', 3, 2)
c
      ma = a - 0.5
      fm = -ma
      m = -ma
c
      te = 1.0
      t = 1.0
      s = t
      do 20 k=1,200
        fkp1 = k + 1
        te = -x*te/(fm+fkp1)
        t = te/fkp1
        s = s + t
        if (abs(t).lt.eps*s) go to 30
 20   continue
      call xermsg ('slatec', 'r9gmic',
     +   'no convergence in 200 terms of continued fraction', 4, 2)
c
 30   r9gmic = -alx - euler + x*s/(fm+1.0)
      if (m.eq.0) return
c
      if (m.eq.1) r9gmic = -r9gmic - 1.0 + 1.0/x
      if (m.eq.1) return
c
      te = fm
      t = 1.0
      s = t
      mm1 = m - 1
      do 40 k=1,mm1
        fk = k
        te = -x*te/fk
        t = te/(fm-fk)
        s = s + t
        if (abs(t).lt.eps*abs(s)) go to 50
 40   continue
c
 50   do 60 k=1,m
        r9gmic = r9gmic + 1.0/k
 60   continue
c
      sgng = 1.0
      if (mod(m,2).eq.1) sgng = -1.0
      alng = log(r9gmic) - alngam(fm+1.0)
c
      r9gmic = 0.0
      if (alng.gt.bot) r9gmic = sgng*exp(alng)
      if (s.ne.0.0) r9gmic = r9gmic + sign (exp(-fm*alx+log(abs(s)/fm))
     1  , s)
c
      if (r9gmic .eq. 0.0 .and. s .eq. 0.0) call xermsg ('slatec',
     +   'r9gmic', 'result underflows', 1, 1)
      return
c
      end
