*deck bvalu
      function bvalu (t, a, n, k, ideriv, x, inbv, work)
c***begin prologue  bvalu
c***purpose  evaluate the b-representation of a b-spline at x for the
c            function value or any of its derivatives.
c***library   slatec
c***category  e3, k6
c***type      single precision (bvalu-s, dbvalu-d)
c***keywords  differentiation of b-spline, evaluation of b-spline
c***author  amos, d. e., (snla)
c***description
c
c     written by carl de boor and modified by d. e. amos
c
c     abstract
c         bvalu is the bvalue function of the reference.
c
c         bvalu evaluates the b-representation (t,a,n,k) of a b-spline
c         at x for the function value on ideriv = 0 or any of its
c         derivatives on ideriv = 1,2,...,k-1.  right limiting values
c         (right derivatives) are returned except at the right end
c         point x=t(n+1) where left limiting values are computed.  the
c         spline is defined on t(k) .le. x .le. t(n+1).  bvalu returns
c         a fatal error message when x is outside of this interval.
c
c         to compute left derivatives or left limiting values at a
c         knot t(i), replace n by i-1 and set x=t(i), i=k+1,n+1.
c
c         bvalu calls intrv
c
c     description of arguments
c         input
c          t       - knot vector of length n+k
c          a       - b-spline coefficient vector of length n
c          n       - number of b-spline coefficients
c                    n = sum of knot multiplicities-k
c          k       - order of the b-spline, k .ge. 1
c          ideriv  - order of the derivative, 0 .le. ideriv .le. k-1
c                    ideriv=0 returns the b-spline value
c          x       - argument, t(k) .le. x .le. t(n+1)
c          inbv    - an initialization parameter which must be set
c                    to 1 the first time bvalu is called.
c
c         output
c          inbv    - inbv contains information for efficient process-
c                    ing after the initial call and inbv must not
c                    be changed by the user.  distinct splines require
c                    distinct inbv parameters.
c          work    - work vector of length 3*k.
c          bvalu   - value of the ideriv-th derivative at x
c
c     error conditions
c         an improper input is a fatal error
c
c***references  carl de boor, package for calculating with b-splines,
c                 siam journal on numerical analysis 14, 3 (june 1977),
c                 pp. 441-472.
c***routines called  intrv, xermsg
c***revision history  (yymmdd)
c   800901  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  bvalu
c
      integer i,ideriv,iderp1,ihi,ihmkmj,ilo,imk,imkpj, inbv, ipj,
     1 ip1, ip1mj, j, jj, j1, j2, k, kmider, kmj, km1, kpk, mflag, n
      real a, fkmj, t, work, x
c     dimension t(n+k), work(3*k)
      dimension t(*), a(*), work(*)
c***first executable statement  bvalu
      bvalu = 0.0e0
      if(k.lt.1) go to 102
      if(n.lt.k) go to 101
      if(ideriv.lt.0 .or. ideriv.ge.k) go to 110
      kmider = k - ideriv
c
c *** find *i* in (k,n) such that t(i) .le. x .lt. t(i+1)
c     (or, .le. t(i+1) if t(i) .lt. t(i+1) = t(n+1)).
      km1 = k - 1
      call intrv(t, n+1, x, inbv, i, mflag)
      if (x.lt.t(k)) go to 120
      if (mflag.eq.0) go to 20
      if (x.gt.t(i)) go to 130
   10 if (i.eq.k) go to 140
      i = i - 1
      if (x.eq.t(i)) go to 10
c
c *** difference the coefficients *ideriv* times
c     work(i) = aj(i), work(k+i) = dp(i), work(k+k+i) = dm(i), i=1.k
c
   20 imk = i - k
      do 30 j=1,k
        imkpj = imk + j
        work(j) = a(imkpj)
   30 continue
      if (ideriv.eq.0) go to 60
      do 50 j=1,ideriv
        kmj = k - j
        fkmj = kmj
        do 40 jj=1,kmj
          ihi = i + jj
          ihmkmj = ihi - kmj
          work(jj) = (work(jj+1)-work(jj))/(t(ihi)-t(ihmkmj))*fkmj
   40   continue
   50 continue
c
c *** compute value at *x* in (t(i),(t(i+1)) of ideriv-th derivative,
c     given its relevant b-spline coeff. in aj(1),...,aj(k-ideriv).
   60 if (ideriv.eq.km1) go to 100
      ip1 = i + 1
      kpk = k + k
      j1 = k + 1
      j2 = kpk + 1
      do 70 j=1,kmider
        ipj = i + j
        work(j1) = t(ipj) - x
        ip1mj = ip1 - j
        work(j2) = x - t(ip1mj)
        j1 = j1 + 1
        j2 = j2 + 1
   70 continue
      iderp1 = ideriv + 1
      do 90 j=iderp1,km1
        kmj = k - j
        ilo = kmj
        do 80 jj=1,kmj
          work(jj) = (work(jj+1)*work(kpk+ilo)+work(jj)
     1              *work(k+jj))/(work(kpk+ilo)+work(k+jj))
          ilo = ilo - 1
   80   continue
   90 continue
  100 bvalu = work(1)
      return
c
c
  101 continue
      call xermsg ('slatec', 'bvalu', 'n does not satisfy n.ge.k', 2,
     +   1)
      return
  102 continue
      call xermsg ('slatec', 'bvalu', 'k does not satisfy k.ge.1', 2,
     +   1)
      return
  110 continue
      call xermsg ('slatec', 'bvalu',
     +   'ideriv does not satisfy 0.le.ideriv.lt.k', 2, 1)
      return
  120 continue
      call xermsg ('slatec', 'bvalu',
     +   'x is n0t greater than or equal to t(k)', 2, 1)
      return
  130 continue
      call xermsg ('slatec', 'bvalu',
     +   'x is not less than or equal to t(n+1)', 2, 1)
      return
  140 continue
      call xermsg ('slatec', 'bvalu',
     +   'a left limiting value cannot be obtained at t(k)', 2, 1)
      return
      end
