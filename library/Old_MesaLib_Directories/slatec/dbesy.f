*deck dbesy
      subroutine dbesy (x, fnu, n, y)
c***begin prologue  dbesy
c***purpose  implement forward recursion on the three term recursion
c            relation for a sequence of non-negative order bessel
c            functions y/sub(fnu+i-1)/(x), i=1,...,n for real, positive
c            x and non-negative orders fnu.
c***library   slatec
c***category  c10a3
c***type      double precision (besy-s, dbesy-d)
c***keywords  special functions, y bessel function
c***author  amos, d. e., (snla)
c***description
c
c     abstract  **** a double precision routine ****
c         dbesy implements forward recursion on the three term
c         recursion relation for a sequence of non-negative order bessel
c         functions y/sub(fnu+i-1)/(x), i=1,n for real x .gt. 0.0d0 and
c         non-negative orders fnu.  if fnu .lt. nulim, orders fnu and
c         fnu+1 are obtained from dbsynu which computes by a power
c         series for x .le. 2, the k bessel function of an imaginary
c         argument for 2 .lt. x .le. 20 and the asymptotic expansion for
c         x .gt. 20.
c
c         if fnu .ge. nulim, the uniform asymptotic expansion is coded
c         in dasyjy for orders fnu and fnu+1 to start the recursion.
c         nulim is 70 or 100 depending on whether n=1 or n .ge. 2.  an
c         overflow test is made on the leading term of the asymptotic
c         expansion before any extensive computation is done.
c
c         the maximum number of significant digits obtainable
c         is the smaller of 14 and the number of digits carried in
c         double precision arithmetic.
c
c     description of arguments
c
c         input
c           x      - x .gt. 0.0d0
c           fnu    - order of the initial y function, fnu .ge. 0.0d0
c           n      - number of members in the sequence, n .ge. 1
c
c         output
c           y      - a vector whose first n components contain values
c                    for the sequence y(i)=y/sub(fnu+i-1)/(x), i=1,n.
c
c     error conditions
c         improper input arguments - a fatal error
c         overflow - a fatal error
c
c***references  f. w. j. olver, tables of bessel functions of moderate
c                 or large orders, npl mathematical tables 6, her
c                 majesty's stationery office, london, 1962.
c               n. m. temme, on the numerical evaluation of the modified
c                 bessel function of the third kind, journal of
c                 computational physics 19, (1975), pp. 324-337.
c               n. m. temme, on the numerical evaluation of the ordinary
c                 bessel function of the second kind, journal of
c                 computational physics 21, (1976), pp. 343-350.
c***routines called  d1mach, dasyjy, dbesy0, dbesy1, dbsynu, dyairy,
c                    i1mach, xermsg
c***revision history  (yymmdd)
c   800501  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890911  removed unnecessary intrinsics.  (wrb)
c   890911  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   920501  reformatted the references section.  (wrb)
c***end prologue  dbesy
c
      external dyairy
      integer i, iflw, j, n, nb, nd, nn, nud, nulim
      integer i1mach
      double precision azn,cn,dnu,elim,flgjy,fn,fnu,ran,s,s1,s2,tm,trx,
     1           w,wk,w2n,x,xlim,xxn,y
      double precision dbesy0, dbesy1, d1mach
      dimension w(2), nulim(2), y(*), wk(7)
      save nulim
      data nulim(1),nulim(2) / 70 , 100 /
c***first executable statement  dbesy
      nn = -i1mach(15)
      elim = 2.303d0*(nn*d1mach(5)-3.0d0)
      xlim = d1mach(1)*1.0d+3
      if (fnu.lt.0.0d0) go to 140
      if (x.le.0.0d0) go to 150
      if (x.lt.xlim) go to 170
      if (n.lt.1) go to 160
c
c     nd is a dummy variable for n
c
      nd = n
      nud = int(fnu)
      dnu = fnu - nud
      nn = min(2,nd)
      fn = fnu + n - 1
      if (fn.lt.2.0d0) go to 100
c
c     overflow test  (leading exponential of asymptotic expansion)
c     for the last order, fnu+n-1.ge.nulim
c
      xxn = x/fn
      w2n = 1.0d0-xxn*xxn
      if(w2n.le.0.0d0) go to 10
      ran = sqrt(w2n)
      azn = log((1.0d0+ran)/xxn) - ran
      cn = fn*azn
      if(cn.gt.elim) go to 170
   10 continue
      if (nud.lt.nulim(nn)) go to 20
c
c     asymptotic expansion for orders fnu and fnu+1.ge.nulim
c
      flgjy = -1.0d0
      call dasyjy(dyairy,x,fnu,flgjy,nn,y,wk,iflw)
      if(iflw.ne.0) go to 170
      if (nn.eq.1) return
      trx = 2.0d0/x
      tm = (fnu+fnu+2.0d0)/x
      go to 80
c
   20 continue
      if (dnu.ne.0.0d0) go to 30
      s1 = dbesy0(x)
      if (nud.eq.0 .and. nd.eq.1) go to 70
      s2 = dbesy1(x)
      go to 40
   30 continue
      nb = 2
      if (nud.eq.0 .and. nd.eq.1) nb = 1
      call dbsynu(x, dnu, nb, w)
      s1 = w(1)
      if (nb.eq.1) go to 70
      s2 = w(2)
   40 continue
      trx = 2.0d0/x
      tm = (dnu+dnu+2.0d0)/x
c     forward recur from dnu to fnu+1 to get y(1) and y(2)
      if (nd.eq.1) nud = nud - 1
      if (nud.gt.0) go to 50
      if (nd.gt.1) go to 70
      s1 = s2
      go to 70
   50 continue
      do 60 i=1,nud
        s = s2
        s2 = tm*s2 - s1
        s1 = s
        tm = tm + trx
   60 continue
      if (nd.eq.1) s1 = s2
   70 continue
      y(1) = s1
      if (nd.eq.1) return
      y(2) = s2
   80 continue
      if (nd.eq.2) return
c     forward recur from fnu+2 to fnu+n-1
      do 90 i=3,nd
        y(i) = tm*y(i-1) - y(i-2)
        tm = tm + trx
   90 continue
      return
c
  100 continue
c     overflow test
      if (fn.le.1.0d0) go to 110
      if (-fn*(log(x)-0.693d0).gt.elim) go to 170
  110 continue
      if (dnu.eq.0.0d0) go to 120
      call dbsynu(x, fnu, nd, y)
      return
  120 continue
      j = nud
      if (j.eq.1) go to 130
      j = j + 1
      y(j) = dbesy0(x)
      if (nd.eq.1) return
      j = j + 1
  130 continue
      y(j) = dbesy1(x)
      if (nd.eq.1) return
      trx = 2.0d0/x
      tm = trx
      go to 80
c
c
c
  140 continue
      call xermsg ('slatec', 'dbesy', 'order, fnu, less than zero', 2,
     +   1)
      return
  150 continue
      call xermsg ('slatec', 'dbesy', 'x less than or equal to zero',
     +   2, 1)
      return
  160 continue
      call xermsg ('slatec', 'dbesy', 'n less than one', 2, 1)
      return
  170 continue
      call xermsg ('slatec', 'dbesy',
     +   'overflow, fnu or n too large or x too small', 6, 1)
      return
      end
