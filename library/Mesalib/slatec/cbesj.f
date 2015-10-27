*deck cbesj
      subroutine cbesj (z, fnu, kode, n, cy, nz, ierr)
c***begin prologue  cbesj
c***purpose  compute a sequence of the bessel functions j(a,z) for
c            complex argument z and real nonnegative orders a=b,b+1,
c            b+2,... where b>0.  a scaling option is available to
c            help avoid overflow.
c***library   slatec
c***category  c10a4
c***type      complex (cbesj-c, zbesj-c)
c***keywords  bessel functions of complex argument,
c             bessel functions of the first kind, j bessel functions
c***author  amos, d. e., (snl)
c***description
c
c         on kode=1, cbesj computes an n member sequence of complex
c         bessel functions cy(l)=j(fnu+l-1,z) for real nonnegative
c         orders fnu+l-1, l=1,...,n and complex z in the cut plane
c         -pi<arg(z)<=pi.  on kode=2, cbesj returns the scaled functions
c
c            cy(l) = exp(-abs(y))*j(fnu+l-1,z),  l=1,...,n and y=im(z)
c
c         which remove the exponential growth in both the upper and
c         lower half planes as z goes to infinity.  definitions and
c         notation are found in the nbs handbook of mathematical
c         functions (ref. 1).
c
c         input
c           z      - argument of type complex
c           fnu    - initial order of type real, fnu>=0
c           kode   - a parameter to indicate the scaling option
c                    kode=1  returns
c                            cy(l)=j(fnu+l-1,z), l=1,...,n
c                        =2  returns
c                            cy(l)=j(fnu+l-1,z)*exp(-abs(y)), l=1,...,n
c                            where y=im(z)
c           n      - number of terms in the sequence, n>=1
c
c         output
c           cy     - result vector of type complex
c           nz     - number of underflows set to zero
c                    nz=0    normal return
c                    nz>0    cy(l)=0, l=n-nz+1,...,n
c           ierr   - error flag
c                    ierr=0  normal return     - computation completed
c                    ierr=1  input error       - no computation
c                    ierr=2  overflow          - no computation
c                            (im(z) too large on kode=1)
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
c         the computation is carried out by the formulae
c
c            j(a,z) = exp( a*pi*i/2)*i(a,-i*z),  im(z)>=0
c
c            j(a,z) = exp(-a*pi*i/2)*i(a, i*z),  im(z)<0
c
c         where the i bessel function is computed as described in the
c         prologue to cbesi.
c
c         for negative orders, the formula
c
c            j(-a,z) = j(a,z)*cos(a*pi) - y(a,z)*sin(a*pi)
c
c         can be used.  however, for large orders close to integers, the
c         the function changes radically.  when a is a large positive
c         integer, the magnitude of j(-a,z)=j(a,z)*cos(a*pi) is a
c         large negative power of ten.  but when a is not an integer,
c         y(a,z) dominates in magnitude with a large positive power of
c         ten and the most that the second term can be reduced is by
c         unit roundoff from the coefficient.  thus, wide changes can
c         occur within unit roundoff of a large integer for a.  here,
c         large means a>abs(z).
c
c         in most complex variable computation, one must evaluate ele-
c         mentary functions.  when the magnitude of z or fnu+n-1 is
c         large, losses of significance by argument reduction occur.
c         consequently, if either one exceeds u1=sqrt(0.5/ur), then
c         losses exceeding half precision are likely and an error flag
c         ierr=3 is triggered where ur=r1mach(4)=unit roundoff.  also,
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
c***routines called  cbinu, i1mach, r1mach
c***revision history  (yymmdd)
c   830501  date written
c   890801  revision date from version 3.2
c   910415  prologue converted to version 4.0 format.  (bab)
c   920128  category corrected.  (wrb)
c   920811  prologue revised.  (dwl)
c***end prologue  cbesj
c
      complex ci, csgn, cy, z, zn
      real aa, alim, arg, dig, elim, fnu, fnul, hpi, rl, r1, r1m5, r2,
     * tol, yy, r1mach, az, fn, bb, ascle, rtol, atol
      integer i, ierr, inu, inuh, ir, kode, k1, k2, n, nl, nz, i1mach, k
      dimension cy(n)
      data hpi /1.57079632679489662e0/
c
c***first executable statement  cbesj
      ierr = 0
      nz=0
      if (fnu.lt.0.0e0) ierr=1
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
      tol = max(r1mach(4),1.0e-18)
      k1 = i1mach(12)
      k2 = i1mach(13)
      r1m5 = r1mach(5)
      k = min(abs(k1),abs(k2))
      elim = 2.303e0*(k*r1m5-3.0e0)
      k1 = i1mach(11) - 1
      aa = r1m5*k1
      dig = min(aa,18.0e0)
      aa = aa*2.303e0
      alim = elim + max(-aa,-41.45e0)
      rl = 1.2e0*dig + 3.0e0
      fnul = 10.0e0 + 6.0e0*(dig-3.0e0)
      ci = cmplx(0.0e0,1.0e0)
      yy = aimag(z)
      az = abs(z)
c-----------------------------------------------------------------------
c     test for range
c-----------------------------------------------------------------------
      aa = 0.5e0/tol
      bb=i1mach(9)*0.5e0
      aa=min(aa,bb)
      fn=fnu+(n-1)
      if(az.gt.aa) go to 140
      if(fn.gt.aa) go to 140
      aa=sqrt(aa)
      if(az.gt.aa) ierr=3
      if(fn.gt.aa) ierr=3
c-----------------------------------------------------------------------
c     calculate csgn=exp(fnu*hpi*i) to minimize losses of significance
c     when fnu is large
c-----------------------------------------------------------------------
      inu = fnu
      inuh = inu/2
      ir = inu - 2*inuh
      arg = (fnu-(inu-ir))*hpi
      r1 = cos(arg)
      r2 = sin(arg)
      csgn = cmplx(r1,r2)
      if (mod(inuh,2).eq.1) csgn = -csgn
c-----------------------------------------------------------------------
c     zn is in the right half plane
c-----------------------------------------------------------------------
      zn = -z*ci
      if (yy.ge.0.0e0) go to 40
      zn = -zn
      csgn = conjg(csgn)
      ci = conjg(ci)
   40 continue
      call cbinu(zn, fnu, kode, n, cy, nz, rl, fnul, tol, elim, alim)
      if (nz.lt.0) go to 120
      nl = n - nz
      if (nl.eq.0) return
      rtol = 1.0e0/tol
      ascle = r1mach(1)*rtol*1.0e+3
      do 50 i=1,nl
c       cy(i)=cy(i)*csgn
        zn=cy(i)
        aa=real(zn)
        bb=aimag(zn)
        atol=1.0e0
        if (max(abs(aa),abs(bb)).gt.ascle) go to 55
          zn = zn*cmplx(rtol,0.0e0)
          atol = tol
   55   continue
        zn = zn*csgn
        cy(i) = zn*cmplx(atol,0.0e0)
        csgn = csgn*ci
   50 continue
      return
  120 continue
      if(nz.eq.(-2)) go to 130
      nz = 0
      ierr = 2
      return
  130 continue
      nz=0
      ierr=5
      return
  140 continue
      nz=0
      ierr=4
      return
      end
