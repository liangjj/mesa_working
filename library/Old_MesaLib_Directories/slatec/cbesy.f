*deck cbesy
      subroutine cbesy (z, fnu, kode, n, cy, nz, cwrk, ierr)
c***begin prologue  cbesy
c***purpose  compute a sequence of the bessel functions y(a,z) for
c            complex argument z and real nonnegative orders a=b,b+1,
c            b+2,... where b>0.  a scaling option is available to
c            help avoid overflow.
c***library   slatec
c***category  c10a4
c***type      complex (cbesy-c, zbesy-c)
c***keywords  bessel functions of complex argument,
c             bessel functions of second kind, weber's function,
c             y bessel functions
c***author  amos, d. e., (snl)
c***description
c
c         on kode=1, cbesy computes an n member sequence of complex
c         bessel functions cy(l)=y(fnu+l-1,z) for real nonnegative
c         orders fnu+l-1, l=1,...,n and complex z in the cut plane
c         -pi<arg(z)<=pi.  on kode=2, cbesy returns the scaled
c         functions
c
c            cy(l) = exp(-abs(y))*y(fnu+l-1,z),  l=1,...,n, y=im(z)
c
c         which remove the exponential growth in both the upper and
c         lower half planes as z goes to infinity.  definitions and
c         notation are found in the nbs handbook of mathematical
c         functions (ref. 1).
c
c         input
c           z      - nonzero argument of type complex
c           fnu    - initial order of type real, fnu>=0
c           kode   - a parameter to indicate the scaling option
c                    kode=1  returns
c                            cy(l)=y(fnu+l-1,z), l=1,...,n
c                        =2  returns
c                            cy(l)=y(fnu+l-1,z)*exp(-abs(y)), l=1,...,n
c                            where y=im(z)
c           n      - number of terms in the sequence, n>=1
c           cwrk   - a work vector of type complex and dimension n
c
c         output
c           cy     - result vector of type complex
c           nz     - number of underflows set to zero
c                    nz=0    normal return
c                    nz>0    cy(l)=0 for nz values of l, usually on
c                            kode=2 (the underflows may not be in an
c                            uninterrupted sequence)
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
c            y(a,z) = (h(1,a,z) - h(2,a,z))/(2*i)
c
c         where the hankel functions are computed as described in cbesh.
c
c         for negative orders, the formula
c
c            y(-a,z) = y(a,z)*cos(a*pi) + j(a,z)*sin(a*pi)
c
c         can be used.  however, for large orders close to half odd
c         integers the function changes radically.  when a is a large
c         positive half odd integer, the magnitude of y(-a,z)=j(a,z)*
c         sin(a*pi) is a large negative power of ten.  but when a is
c         not a half odd integer, y(a,z) dominates in magnitude with a
c         large positive power of ten and the most that the second term
c         can be reduced is by unit roundoff from the coefficient.
c         thus,  wide changes can occur within unit roundoff of a large
c         half odd integer.  here, large means a>abs(z).
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
c***routines called  cbesh, i1mach, r1mach
c***revision history  (yymmdd)
c   830501  date written
c   890801  revision date from version 3.2
c   910415  prologue converted to version 4.0 format.  (bab)
c   920128  category corrected.  (wrb)
c   920811  prologue revised.  (dwl)
c***end prologue  cbesy
c
      complex cwrk, cy, c1, c2, ex, hci, z, zu, zv
      real elim, ey, fnu, r1, r2, tay, xx, yy, r1mach, r1m5, ascle,
     *  rtol, atol, tol, aa, bb
      integer i, ierr, k, kode, k1, k2, n, nz, nz1, nz2, i1mach
      dimension cy(n), cwrk(n)
c***first executable statement  cbesy
      xx = real(z)
      yy = aimag(z)
      ierr = 0
      nz=0
      if (xx.eq.0.0e0 .and. yy.eq.0.0e0) ierr=1
      if (fnu.lt.0.0e0) ierr=1
      if (kode.lt.1 .or. kode.gt.2) ierr=1
      if (n.lt.1) ierr=1
      if (ierr.ne.0) return
      hci = cmplx(0.0e0,0.5e0)
      call cbesh(z, fnu, kode, 1, n, cy, nz1, ierr)
      if (ierr.ne.0.and.ierr.ne.3) go to 170
      call cbesh(z, fnu, kode, 2, n, cwrk, nz2, ierr)
      if (ierr.ne.0.and.ierr.ne.3) go to 170
      nz = min(nz1,nz2)
      if (kode.eq.2) go to 60
      do 50 i=1,n
        cy(i) = hci*(cwrk(i)-cy(i))
   50 continue
      return
   60 continue
      tol = max(r1mach(4),1.0e-18)
      k1 = i1mach(12)
      k2 = i1mach(13)
      k = min(abs(k1),abs(k2))
      r1m5 = r1mach(5)
c-----------------------------------------------------------------------
c     elim is the approximate exponential under- and overflow limit
c-----------------------------------------------------------------------
      elim = 2.303e0*(k*r1m5-3.0e0)
      r1 = cos(xx)
      r2 = sin(xx)
      ex = cmplx(r1,r2)
      ey = 0.0e0
      tay = abs(yy+yy)
      if (tay.lt.elim) ey = exp(-tay)
      if (yy.lt.0.0e0) go to 90
      c1 = ex*cmplx(ey,0.0e0)
      c2 = conjg(ex)
   70 continue
      nz = 0
      rtol = 1.0e0/tol
      ascle = r1mach(1)*rtol*1.0e+3
      do 80 i=1,n
c       cy(i) = hci*(c2*cwrk(i)-c1*cy(i))
        zv = cwrk(i)
        aa=real(zv)
        bb=aimag(zv)
        atol=1.0e0
        if (max(abs(aa),abs(bb)).gt.ascle) go to 75
          zv = zv*cmplx(rtol,0.0e0)
          atol = tol
   75   continue
        zv = zv*c2*hci
        zv = zv*cmplx(atol,0.0e0)
        zu=cy(i)
        aa=real(zu)
        bb=aimag(zu)
        atol=1.0e0
        if (max(abs(aa),abs(bb)).gt.ascle) go to 85
          zu = zu*cmplx(rtol,0.0e0)
          atol = tol
   85   continue
        zu = zu*c1*hci
        zu = zu*cmplx(atol,0.0e0)
        cy(i) = zv - zu
        if (cy(i).eq.cmplx(0.0e0,0.0e0) .and. ey.eq.0.0e0) nz = nz + 1
   80 continue
      return
   90 continue
      c1 = ex
      c2 = conjg(ex)*cmplx(ey,0.0e0)
      go to 70
  170 continue
      nz = 0
      return
      end
