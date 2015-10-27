*deck zbesh
      subroutine zbesh (zr, zi, fnu, kode, m, n, cyr, cyi, nz, ierr)
c***begin prologue  zbesh
c***purpose  compute a sequence of the hankel functions h(m,a,z)
c            for superscript m=1 or 2, real nonnegative orders a=b,
c            b+1,... where b>0, and nonzero complex argument z.  a
c            scaling option is available to help avoid overflow.
c***library   slatec
c***category  c10a4
c***type      complex (cbesh-c, zbesh-c)
c***keywords  bessel functions of complex argument,
c             bessel functions of the third kind, h bessel functions,
c             hankel functions
c***author  amos, d. e., (snl)
c***description
c
c                      ***a double precision routine***
c         on kode=1, zbesh computes an n member sequence of complex
c         hankel (bessel) functions cy(l)=h(m,fnu+l-1,z) for super-
c         script m=1 or 2, real nonnegative orders fnu+l-1, l=1,...,
c         n, and complex nonzero z in the cut plane -pi<arg(z)<=pi
c         where z=zr+i*zi.  on kode=2, cbesh returns the scaled
c         functions
c
c            cy(l) = h(m,fnu+l-1,z)*exp(-(3-2*m)*z*i),  i**2=-1
c
c         which removes the exponential behavior in both the upper
c         and lower half planes.  definitions and notation are found
c         in the nbs handbook of mathematical functions (ref. 1).
c
c         input
c           zr     - double precision real part of nonzero argument z
c           zi     - double precision imag part of nonzero argument z
c           fnu    - double precision initial order, fnu>=0
c           kode   - a parameter to indicate the scaling option
c                    kode=1  returns
c                            cy(l)=h(m,fnu+l-1,z), l=1,...,n
c                        =2  returns
c                            cy(l)=h(m,fnu+l-1,z)*exp(-(3-2m)*z*i),
c                            l=1,...,n
c           m      - superscript of hankel function, m=1 or 2
c           n      - number of terms in the sequence, n>=1
c
c         output
c           cyr    - double precision real part of result vector
c           cyi    - double precision imag part of result vector
c           nz     - number of underflows set to zero
c                    nz=0    normal return
c                    nz>0    cy(l)=0 for nz values of l (if m=1 and
c                            im(z)>0 or if m=2 and im(z)<0, then
c                            cy(l)=0 for l=1,...,nz; in the com-
c                            plementary half planes, the underflows
c                            may not be in an uninterrupted sequence)
c           ierr   - error flag
c                    ierr=0  normal return     - computation completed
c                    ierr=1  input error       - no computation
c                    ierr=2  overflow          - no computation
c                            (abs(z) too small and/or fnu+n-1
c                            too large)
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
c         the computation is carried out by the formula
c
c            h(m,a,z) = (1/t)*exp(-a*t)*k(a,z*exp(-t))
c                   t = (3-2*m)*i*pi/2
c
c         where the k bessel function is computed as described in the
c         prologue to cbesk.
c
c         exponential decay of h(m,a,z) occurs in the upper half z
c         plane for m=1 and the lower half z plane for m=2.  exponential
c         growth occurs in the complementary half planes.  scaling
c         by exp(-(3-2*m)*z*i) removes the exponential behavior in the
c         whole z plane as z goes to infinity.
c
c         for negative orders, the formula
c
c            h(m,-a,z) = h(m,a,z)*exp((3-2*m)*a*pi*i)
c
c         can be used.
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
c***routines called  d1mach, i1mach, zabs, zacon, zbknu, zbunk, zuoik
c***revision history  (yymmdd)
c   830501  date written
c   890801  revision date from version 3.2
c   910415  prologue converted to version 4.0 format.  (bab)
c   920128  category corrected.  (wrb)
c   920811  prologue revised.  (dwl)
c***end prologue  zbesh
c
c     complex cy,z,zn,zt,csgn
      double precision aa, alim, aln, arg, az, cyi, cyr, dig, elim,
     * fmm, fn, fnu, fnul, hpi, rhpi, rl, r1m5, sgn, str, tol, ufl, zi,
     * zni, znr, zr, zti, d1mach, zabs, bb, ascle, rtol, atol, sti,
     * csgnr, csgni
      integer i, ierr, inu, inuh, ir, k, kode, k1, k2, m,
     * mm, mr, n, nn, nuf, nw, nz, i1mach
      dimension cyr(n), cyi(n)
      external zabs
c
      data hpi /1.57079632679489662d0/
c
c***first executable statement  zbesh
      ierr = 0
      nz=0
      if (zr.eq.0.0d0 .and. zi.eq.0.0d0) ierr=1
      if (fnu.lt.0.0d0) ierr=1
      if (m.lt.1 .or. m.gt.2) ierr=1
      if (kode.lt.1 .or. kode.gt.2) ierr=1
      if (n.lt.1) ierr=1
      if (ierr.ne.0) return
      nn = n
c-----------------------------------------------------------------------
c     set parameters related to machine constants.
c     tol is the approximate unit roundoff limited to 1.0e-18.
c     elim is the approximate exponential over- and underflow limit.
c     exp(-elim).lt.exp(-alim)=exp(-elim)/tol    and
c     exp(elim).gt.exp(alim)=exp(elim)*tol       are intervals near
c     underflow and overflow limits where scaled arithmetic is done.
c     rl is the lower boundary of the asymptotic expansion for large z.
c     dig = number of base 10 digits in tol = 10**(-dig).
c     fnul is the lower boundary of the asymptotic series for large fnu
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
      fnul = 10.0d0 + 6.0d0*(dig-3.0d0)
      rl = 1.2d0*dig + 3.0d0
      fn = fnu + (nn-1)
      mm = 3 - m - m
      fmm = mm
      znr = fmm*zi
      zni = -fmm*zr
c-----------------------------------------------------------------------
c     test for proper range
c-----------------------------------------------------------------------
      az = zabs(zr,zi)
      aa = 0.5d0/tol
      bb = i1mach(9)*0.5d0
      aa = min(aa,bb)
      if (az.gt.aa) go to 260
      if (fn.gt.aa) go to 260
      aa = sqrt(aa)
      if (az.gt.aa) ierr=3
      if (fn.gt.aa) ierr=3
c-----------------------------------------------------------------------
c     overflow test on the last member of the sequence
c-----------------------------------------------------------------------
      ufl = d1mach(1)*1.0d+3
      if (az.lt.ufl) go to 230
      if (fnu.gt.fnul) go to 90
      if (fn.le.1.0d0) go to 70
      if (fn.gt.2.0d0) go to 60
      if (az.gt.tol) go to 70
      arg = 0.5d0*az
      aln = -fn*log(arg)
      if (aln.gt.elim) go to 230
      go to 70
   60 continue
      call zuoik(znr, zni, fnu, kode, 2, nn, cyr, cyi, nuf, tol, elim,
     * alim)
      if (nuf.lt.0) go to 230
      nz = nz + nuf
      nn = nn - nuf
c-----------------------------------------------------------------------
c     here nn=n or nn=0 since nuf=0,nn, or -1 on return from cuoik
c     if nuf=nn, then cy(i)=czero for all i
c-----------------------------------------------------------------------
      if (nn.eq.0) go to 140
   70 continue
      if ((znr.lt.0.0d0) .or. (znr.eq.0.0d0 .and. zni.lt.0.0d0 .and.
     * m.eq.2)) go to 80
c-----------------------------------------------------------------------
c     right half plane computation, xn.ge.0. .and. (xn.ne.0. .or.
c     yn.ge.0. .or. m=1)
c-----------------------------------------------------------------------
      call zbknu(znr, zni, fnu, kode, nn, cyr, cyi, nz, tol, elim, alim)
      go to 110
c-----------------------------------------------------------------------
c     left half plane computation
c-----------------------------------------------------------------------
   80 continue
      mr = -mm
      call zacon(znr, zni, fnu, kode, mr, nn, cyr, cyi, nw, rl, fnul,
     * tol, elim, alim)
      if (nw.lt.0) go to 240
      nz=nw
      go to 110
   90 continue
c-----------------------------------------------------------------------
c     uniform asymptotic expansions for fnu.gt.fnul
c-----------------------------------------------------------------------
      mr = 0
      if ((znr.ge.0.0d0) .and. (znr.ne.0.0d0 .or. zni.ge.0.0d0 .or.
     * m.ne.2)) go to 100
      mr = -mm
      if (znr.ne.0.0d0 .or. zni.ge.0.0d0) go to 100
      znr = -znr
      zni = -zni
  100 continue
      call zbunk(znr, zni, fnu, kode, mr, nn, cyr, cyi, nw, tol, elim,
     * alim)
      if (nw.lt.0) go to 240
      nz = nz + nw
  110 continue
c-----------------------------------------------------------------------
c     h(m,fnu,z) = -fmm*(i/hpi)*(zt**fnu)*k(fnu,-z*zt)
c
c     zt=exp(-fmm*hpi*i) = cmplx(0.0,-fmm), fmm=3-2*m, m=1,2
c-----------------------------------------------------------------------
      sgn = dsign(hpi,-fmm)
c-----------------------------------------------------------------------
c     calculate exp(fnu*hpi*i) to minimize losses of significance
c     when fnu is large
c-----------------------------------------------------------------------
      inu = fnu
      inuh = inu/2
      ir = inu - 2*inuh
      arg = (fnu-(inu-ir))*sgn
      rhpi = 1.0d0/sgn
c     zni = rhpi*cos(arg)
c     znr = -rhpi*sin(arg)
      csgni = rhpi*cos(arg)
      csgnr = -rhpi*sin(arg)
      if (mod(inuh,2).eq.0) go to 120
c     znr = -znr
c     zni = -zni
      csgnr = -csgnr
      csgni = -csgni
  120 continue
      zti = -fmm
      rtol = 1.0d0/tol
      ascle = ufl*rtol
      do 130 i=1,nn
c       str = cyr(i)*znr - cyi(i)*zni
c       cyi(i) = cyr(i)*zni + cyi(i)*znr
c       cyr(i) = str
c       str = -zni*zti
c       zni = znr*zti
c       znr = str
        aa = cyr(i)
        bb = cyi(i)
        atol = 1.0d0
        if (max(abs(aa),abs(bb)).gt.ascle) go to 135
          aa = aa*rtol
          bb = bb*rtol
          atol = tol
  135 continue
      str = aa*csgnr - bb*csgni
      sti = aa*csgni + bb*csgnr
      cyr(i) = str*atol
      cyi(i) = sti*atol
      str = -csgni*zti
      csgni = csgnr*zti
      csgnr = str
  130 continue
      return
  140 continue
      if (znr.lt.0.0d0) go to 230
      return
  230 continue
      nz=0
      ierr=2
      return
  240 continue
      if(nw.eq.(-1)) go to 230
      nz=0
      ierr=5
      return
  260 continue
      nz=0
      ierr=4
      return
      end
