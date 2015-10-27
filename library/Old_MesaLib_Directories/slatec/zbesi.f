*deck zbesi
      subroutine zbesi (zr, zi, fnu, kode, n, cyr, cyi, nz, ierr)
c***begin prologue  zbesi
c***purpose  compute a sequence of the bessel functions i(a,z) for
c            complex argument z and real nonnegative orders a=b,b+1,
c            b+2,... where b>0.  a scaling option is available to
c            help avoid overflow.
c***library   slatec
c***category  c10b4
c***type      complex (cbesi-c, zbesi-c)
c***keywords  bessel functions of complex argument, i bessel functions,
c             modified bessel functions
c***author  amos, d. e., (snl)
c***description
c
c                    ***a double precision routine***
c         on kode=1, zbesi computes an n-member sequence of complex
c         bessel functions cy(l)=i(fnu+l-1,z) for real nonnegative
c         orders fnu+l-1, l=1,...,n and complex z in the cut plane
c         -pi<arg(z)<=pi where z=zr+i*zi.  on kode=2, cbesi returns
c         the scaled functions
c
c            cy(l) = exp(-abs(x))*i(fnu+l-1,z), l=1,...,n and x=re(z)
c
c         which removes the exponential growth in both the left and
c         right half-planes as z goes to infinity.
c
c         input
c           zr     - double precision real part of argument z
c           zi     - double precision imag part of argument z
c           fnu    - double precision initial order, fnu>=0
c           kode   - a parameter to indicate the scaling option
c                    kode=1  returns
c                            cy(l)=i(fnu+l-1,z), l=1,...,n
c                        =2  returns
c                            cy(l)=exp(-abs(x))*i(fnu+l-1,z), l=1,...,n
c                            where x=re(z)
c           n      - number of terms in the sequence, n>=1
c
c         output
c           cyr    - double precision real part of result vector
c           cyi    - double precision imag part of result vector
c           nz     - number of underflows set to zero
c                    nz=0    normal return
c                    nz>0    cy(l)=0, l=n-nz+1,...,n
c           ierr   - error flag
c                    ierr=0  normal return     - computation completed
c                    ierr=1  input error       - no computation
c                    ierr=2  overflow          - no computation
c                            (re(z) too large on kode=1)
c                    ierr=3  precision warning - computation completed
c                            (result has half precision or less
c                            because abs(z) or fnu+n-1 is large)
c                    ierr=4  precision error   - no computation
c                            (result has no precision because
c                            abs(z) or fnu+n-1 is too large)
c                    ierr=5  algorithmic error - no computation
c                            (termination condition not met)
c
c *long description:
c
c         the computation of i(a,z) is carried out by the power series
c         for small abs(z), the asymptotic expansion for large abs(z),
c         the miller algorithm normalized by the wronskian and a
c         neumann series for intermediate magnitudes of z, and the
c         uniform asymptotic expansions for i(a,z) and j(a,z) for
c         large orders a.  backward recurrence is used to generate
c         sequences or reduce orders when necessary.
c
c         the calculations above are done in the right half plane and
c         continued into the left half plane by the formula
c
c            i(a,z*exp(t)) = exp(t*a)*i(a,z), re(z)>0
c                        t = i*pi or -i*pi
c
c         for negative orders, the formula
c
c            i(-a,z) = i(a,z) + (2/pi)*sin(pi*a)*k(a,z)
c
c         can be used.  however, for large orders close to integers the
c         the function changes radically.  when a is a large positive
c         integer, the magnitude of i(-a,z)=i(a,z) is a large
c         negative power of ten. but when a is not an integer,
c         k(a,z) dominates in magnitude with a large positive power of
c         ten and the most that the second term can be reduced is by
c         unit roundoff from the coefficient. thus, wide changes can
c         occur within unit roundoff of a large integer for a. here,
c         large means a>abs(z).
c
c         in most complex variable computation, one must evaluate ele-
c         mentary functions.  when the magnitude of z or fnu+n-1 is
c         large, losses of significance by argument reduction occur.
c         consequently, if either one exceeds u1=sqrt(0.5/ur), then
c         losses exceeding half precision are likely and an error flag
c         ierr=3 is triggered where ur=max(d1mach(4),1.0d-18) is double
c         precision unit roundoff limited to 18 digits precision.  also,
c         if either is larger than u2=0.5/ur, then all significance is
c         lost and ierr=4.  in order to use the int function, arguments
c         must be further restricted not to exceed the largest machine
c         integer, u3=i1mach(9).  thus, the magnitude of z and fnu+n-1
c         is restricted by min(u2,u3).  in ieee arithmetic, u1,u2, and
c         u3 approximate 2.0e+3, 4.2e+6, 2.1e+9 in single precision
c         and 4.7e+7, 2.3e+15 and 2.1e+9 in double precision.  this
c         makes u2 limiting in single precision and u3 limiting in
c         double precision.  this means that one can expect to retain,
c         in the worst cases on ieee machines, no digits in single pre-
c         cision and only 6 digits in double precision.  similar con-
c         siderations hold for other machines.
c
c         the approximate relative error in the magnitude of a complex
c         bessel function can be expressed as p*10**s where p=max(unit
c         roundoff,1.0e-18) is the nominal precision and 10**s repre-
c         sents the increase in error due to argument reduction in the
c         elementary functions.  here, s=max(1,abs(log10(abs(z))),
c         abs(log10(fnu))) approximately (i.e., s=max(1,abs(exponent of
c         abs(z),abs(exponent of fnu)) ).  however, the phase angle may
c         have only absolute accuracy.  this is most likely to occur
c         when one component (in magnitude) is larger than the other by
c         several orders of magnitude.  if one component is 10**k larger
c         than the other, then one can expect only max(abs(log10(p))-k,
c         0) significant digits; or, stated another way, when k exceeds
c         the exponent of p, no significant digits remain in the smaller
c         component.  however, the phase angle retains absolute accuracy
c         because, in complex arithmetic with precision p, the smaller
c         component will not (as a rule) decrease below p times the
c         magnitude of the larger component.  in these extreme cases,
c         the principal phase angle is on the order of +p, -p, pi/2-p,
c         or -pi/2+p.
c
c***references  1. m. abramowitz and i. a. stegun, handbook of mathe-
c                 matical functions, national bureau of standards
c                 applied mathematics series 55, u. s. department
c                 of commerce, tenth printing (1972) or later.
c               2. d. e. amos, computation of bessel functions of
c                 complex argument, report sand83-0086, sandia national
c                 laboratories, albuquerque, nm, may 1983.
c               3. d. e. amos, computation of bessel functions of
c                 complex argument and large order, report sand83-0643,
c                 sandia national laboratories, albuquerque, nm, may
c                 1983.
c               4. d. e. amos, a subroutine package for bessel functions
c                 of a complex argument and nonnegative order, report
c                 sand85-1018, sandia national laboratory, albuquerque,
c                 nm, may 1985.
c               5. d. e. amos, a portable package for bessel functions
c                 of a complex argument and nonnegative order, acm
c                 transactions on mathematical software, 12 (september
c                 1986), pp. 265-273.
c
c***routines called  d1mach, i1mach, zabs, zbinu
c***revision history  (yymmdd)
c   830501  date written
c   890801  revision date from version 3.2
c   910415  prologue converted to version 4.0 format.  (bab)
c   920128  category corrected.  (wrb)
c   920811  prologue revised.  (dwl)
c***end prologue  zbesi
c     complex cone,csgn,cw,cy,czero,z,zn
      double precision aa, alim, arg, conei, coner, csgni, csgnr, cyi,
     * cyr, dig, elim, fnu, fnul, pi, rl, r1m5, str, tol, zi, zni, znr,
     * zr, d1mach, az, bb, fn, zabs, ascle, rtol, atol, sti
      integer i, ierr, inu, k, kode, k1,k2,n,nz,nn, i1mach
      dimension cyr(n), cyi(n)
      external zabs
      data pi /3.14159265358979324d0/
      data coner, conei /1.0d0,0.0d0/
c
c***first executable statement  zbesi
      ierr = 0
      nz=0
      if (fnu.lt.0.0d0) ierr=1
      if (kode.lt.1 .or. kode.gt.2) ierr=1
      if (n.lt.1) ierr=1
      if (ierr.ne.0) return
c-----------------------------------------------------------------------
c     set parameters related to machine constants.
c     tol is the approximate unit roundoff limited to 1.0e-18.
c     elim is the approximate exponential over- and underflow limit.
c     exp(-elim).lt.exp(-alim)=exp(-elim)/tol    and
c     exp(elim).gt.exp(alim)=exp(elim)*tol       are intervals near
c     underflow and overflow limits where scaled arithmetic is done.
c     rl is the lower boundary of the asymptotic expansion for large z.
c     dig = number of base 10 digits in tol = 10**(-dig).
c     fnul is the lower boundary of the asymptotic series for large fnu.
c-----------------------------------------------------------------------
      tol = max(d1mach(4),1.0d-18)
      k1 = i1mach(15)
      k2 = i1mach(16)
      r1m5 = d1mach(5)
      k = min(abs(k1),abs(k2))
      elim = 2.303d0*(k*r1m5-3.0d0)
      k1 = i1mach(14) - 1
      aa = r1m5*k1
      dig = min(aa,18.0d0)
      aa = aa*2.303d0
      alim = elim + max(-aa,-41.45d0)
      rl = 1.2d0*dig + 3.0d0
      fnul = 10.0d0 + 6.0d0*(dig-3.0d0)
c-----------------------------------------------------------------------
c     test for proper range
c-----------------------------------------------------------------------
      az = zabs(zr,zi)
      fn = fnu+(n-1)
      aa = 0.5d0/tol
      bb=i1mach(9)*0.5d0
      aa = min(aa,bb)
      if (az.gt.aa) go to 260
      if (fn.gt.aa) go to 260
      aa = sqrt(aa)
      if (az.gt.aa) ierr=3
      if (fn.gt.aa) ierr=3
      znr = zr
      zni = zi
      csgnr = coner
      csgni = conei
      if (zr.ge.0.0d0) go to 40
      znr = -zr
      zni = -zi
c-----------------------------------------------------------------------
c     calculate csgn=exp(fnu*pi*i) to minimize losses of significance
c     when fnu is large
c-----------------------------------------------------------------------
      inu = fnu
      arg = (fnu-inu)*pi
      if (zi.lt.0.0d0) arg = -arg
      csgnr = cos(arg)
      csgni = sin(arg)
      if (mod(inu,2).eq.0) go to 40
      csgnr = -csgnr
      csgni = -csgni
   40 continue
      call zbinu(znr, zni, fnu, kode, n, cyr, cyi, nz, rl, fnul, tol,
     * elim, alim)
      if (nz.lt.0) go to 120
      if (zr.ge.0.0d0) return
c-----------------------------------------------------------------------
c     analytic continuation to the left half plane
c-----------------------------------------------------------------------
      nn = n - nz
      if (nn.eq.0) return
      rtol = 1.0d0/tol
      ascle = d1mach(1)*rtol*1.0d+3
      do 50 i=1,nn
c       str = cyr(i)*csgnr - cyi(i)*csgni
c       cyi(i) = cyr(i)*csgni + cyi(i)*csgnr
c       cyr(i) = str
        aa = cyr(i)
        bb = cyi(i)
        atol = 1.0d0
        if (max(abs(aa),abs(bb)).gt.ascle) go to 55
          aa = aa*rtol
          bb = bb*rtol
          atol = tol
   55   continue
        str = aa*csgnr - bb*csgni
        sti = aa*csgni + bb*csgnr
        cyr(i) = str*atol
        cyi(i) = sti*atol
        csgnr = -csgnr
        csgni = -csgni
   50 continue
      return
  120 continue
      if(nz.eq.(-2)) go to 130
      nz = 0
      ierr=2
      return
  130 continue
      nz=0
      ierr=5
      return
  260 continue
      nz=0
      ierr=4
      return
      end
