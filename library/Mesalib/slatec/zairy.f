*deck zairy
      subroutine zairy (zr, zi, id, kode, air, aii, nz, ierr)
c***begin prologue  zairy
c***purpose  compute the airy function ai(z) or its derivative dai/dz
c            for complex argument z.  a scaling option is available
c            to help avoid underflow and overflow.
c***library   slatec
c***category  c10d
c***type      complex (cairy-c, zairy-c)
c***keywords  airy function, bessel function of order one third,
c             bessel function of order two thirds
c***author  amos, d. e., (snl)
c***description
c
c                      ***a double precision routine***
c         on kode=1, zairy computes the complex airy function ai(z)
c         or its derivative dai/dz on id=0 or id=1 respectively. on
c         kode=2, a scaling option exp(zeta)*ai(z) or exp(zeta)*dai/dz
c         is provided to remove the exponential decay in -pi/3<arg(z)
c         <pi/3 and the exponential growth in pi/3<abs(arg(z))<pi where
c         zeta=(2/3)*z**(3/2).
c
c         while the airy functions ai(z) and dai/dz are analytic in
c         the whole z-plane, the corresponding scaled functions defined
c         for kode=2 have a cut along the negative real axis.
c
c         input
c           zr     - double precision real part of argument z
c           zi     - double precision imag part of argument z
c           id     - order of derivative, id=0 or id=1
c           kode   - a parameter to indicate the scaling option
c                    kode=1  returns
c                            ai=ai(z)  on id=0
c                            ai=dai/dz on id=1
c                            at z=z
c                        =2  returns
c                            ai=exp(zeta)*ai(z)  on id=0
c                            ai=exp(zeta)*dai/dz on id=1
c                            at z=z where zeta=(2/3)*z**(3/2)
c
c         output
c           air    - double precision real part of result
c           aii    - double precision imag part of result
c           nz     - underflow indicator
c                    nz=0    normal return
c                    nz=1    ai=0 due to underflow in
c                            -pi/3<arg(z)<pi/3 on kode=1
c           ierr   - error flag
c                    ierr=0  normal return     - computation completed
c                    ierr=1  input error       - no computation
c                    ierr=2  overflow          - no computation
c                            (re(z) too large with kode=1)
c                    ierr=3  precision warning - computation completed
c                            (result has less than half precision)
c                    ierr=4  precision error   - no computation
c                            (result has no precision)
c                    ierr=5  algorithmic error - no computation
c                            (termination condition not met)
c
c *long description:
c
c         ai(z) and dai/dz are computed from k bessel functions by
c
c                ai(z) =  c*sqrt(z)*k(1/3,zeta)
c               dai/dz = -c*   z   *k(2/3,zeta)
c                    c =  1/(pi*sqrt(3))
c                 zeta =  (2/3)*z**(3/2)
c
c         when abs(z)>1 and from power series when abs(z)<=1.
c
c         in most complex variable computation, one must evaluate ele-
c         mentary functions.  when the magnitude of z is large, losses
c         of significance by argument reduction occur.  consequently, if
c         the magnitude of zeta=(2/3)*z**(3/2) exceeds u1=sqrt(0.5/ur),
c         then losses exceeding half precision are likely and an error
c         flag ierr=3 is triggered where ur=max(d1mach(4),1.0d-18) is
c         double precision unit roundoff limited to 18 digits precision.
c         also, if the magnitude of zeta is larger than u2=0.5/ur, then
c         all significance is lost and ierr=4.  in order to use the int
c         function, zeta must be further restricted not to exceed
c         u3=i1mach(9)=largest integer.  thus, the magnitude of zeta
c         must be restricted by min(u2,u3).  in ieee arithmetic, u1,u2,
c         and u3 are approximately 2.0e+3, 4.2e+6, 2.1e+9 in single
c         precision and 4.7e+7, 2.3e+15, 2.1e+9 in double precision.
c         this makes u2 limiting is single precision and u3 limiting
c         in double precision.  this means that the magnitude of z
c         cannot exceed approximately 3.4e+4 in single precision and
c         2.1e+6 in double precision.  this also means that one can
c         expect to retain, in the worst cases on 32-bit machines,
c         no digits in single precision and only 6 digits in double
c         precision.
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
c         magnitude of the larger component. in these extreme cases,
c         the principal phase angle is on the order of +p, -p, pi/2-p,
c         or -pi/2+p.
c
c***references  1. m. abramowitz and i. a. stegun, handbook of mathe-
c                 matical functions, national bureau of standards
c                 applied mathematics series 55, u. s. department
c                 of commerce, tenth printing (1972) or later.
c               2. d. e. amos, computation of bessel functions of
c                 complex argument and large order, report sand83-0643,
c                 sandia national laboratories, albuquerque, nm, may
c                 1983.
c               3. d. e. amos, a subroutine package for bessel functions
c                 of a complex argument and nonnegative order, report
c                 sand85-1018, sandia national laboratory, albuquerque,
c                 nm, may 1985.
c               4. d. e. amos, a portable package for bessel functions
c                 of a complex argument and nonnegative order, acm
c                 transactions on mathematical software, 12 (september
c                 1986), pp. 265-273.
c
c***routines called  d1mach, i1mach, zabs, zacai, zbknu, zexp, zsqrt
c***revision history  (yymmdd)
c   830501  date written
c   890801  revision date from version 3.2
c   910415  prologue converted to version 4.0 format.  (bab)
c   920128  category corrected.  (wrb)
c   920811  prologue revised.  (dwl)
c   930122  added zexp and zsqrt to external statement.  (rwc)
c***end prologue  zairy
c     complex ai,cone,csq,cy,s1,s2,trm1,trm2,z,zta,z3
      double precision aa, ad, aii, air, ak, alim, atrm, az, az3, bk,
     * cc, ck, coef, conei, coner, csqi, csqr, cyi, cyr, c1, c2, dig,
     * dk, d1, d2, elim, fid, fnu, ptr, rl, r1m5, sfac, sti, str,
     * s1i, s1r, s2i, s2r, tol, trm1i, trm1r, trm2i, trm2r, tth, zeroi,
     * zeror, zi, zr, ztai, ztar, z3i, z3r, d1mach, zabs, alaz, bb
      integer id, ierr, iflag, k, kode, k1, k2, mr, nn, nz, i1mach
      dimension cyr(1), cyi(1)
      external zabs, zexp, zsqrt
      data tth, c1, c2, coef /6.66666666666666667d-01,
     * 3.55028053887817240d-01,2.58819403792806799d-01,
     * 1.83776298473930683d-01/
      data zeror, zeroi, coner, conei /0.0d0,0.0d0,1.0d0,0.0d0/
c***first executable statement  zairy
      ierr = 0
      nz=0
      if (id.lt.0 .or. id.gt.1) ierr=1
      if (kode.lt.1 .or. kode.gt.2) ierr=1
      if (ierr.ne.0) return
      az = zabs(zr,zi)
      tol = max(d1mach(4),1.0d-18)
      fid = id
      if (az.gt.1.0d0) go to 70
c-----------------------------------------------------------------------
c     power series for abs(z).le.1.
c-----------------------------------------------------------------------
      s1r = coner
      s1i = conei
      s2r = coner
      s2i = conei
      if (az.lt.tol) go to 170
      aa = az*az
      if (aa.lt.tol/az) go to 40
      trm1r = coner
      trm1i = conei
      trm2r = coner
      trm2i = conei
      atrm = 1.0d0
      str = zr*zr - zi*zi
      sti = zr*zi + zi*zr
      z3r = str*zr - sti*zi
      z3i = str*zi + sti*zr
      az3 = az*aa
      ak = 2.0d0 + fid
      bk = 3.0d0 - fid - fid
      ck = 4.0d0 - fid
      dk = 3.0d0 + fid + fid
      d1 = ak*dk
      d2 = bk*ck
      ad = min(d1,d2)
      ak = 24.0d0 + 9.0d0*fid
      bk = 30.0d0 - 9.0d0*fid
      do 30 k=1,25
        str = (trm1r*z3r-trm1i*z3i)/d1
        trm1i = (trm1r*z3i+trm1i*z3r)/d1
        trm1r = str
        s1r = s1r + trm1r
        s1i = s1i + trm1i
        str = (trm2r*z3r-trm2i*z3i)/d2
        trm2i = (trm2r*z3i+trm2i*z3r)/d2
        trm2r = str
        s2r = s2r + trm2r
        s2i = s2i + trm2i
        atrm = atrm*az3/ad
        d1 = d1 + ak
        d2 = d2 + bk
        ad = min(d1,d2)
        if (atrm.lt.tol*ad) go to 40
        ak = ak + 18.0d0
        bk = bk + 18.0d0
   30 continue
   40 continue
      if (id.eq.1) go to 50
      air = s1r*c1 - c2*(zr*s2r-zi*s2i)
      aii = s1i*c1 - c2*(zr*s2i+zi*s2r)
      if (kode.eq.1) return
      call zsqrt(zr, zi, str, sti)
      ztar = tth*(zr*str-zi*sti)
      ztai = tth*(zr*sti+zi*str)
      call zexp(ztar, ztai, str, sti)
      ptr = air*str - aii*sti
      aii = air*sti + aii*str
      air = ptr
      return
   50 continue
      air = -s2r*c2
      aii = -s2i*c2
      if (az.le.tol) go to 60
      str = zr*s1r - zi*s1i
      sti = zr*s1i + zi*s1r
      cc = c1/(1.0d0+fid)
      air = air + cc*(str*zr-sti*zi)
      aii = aii + cc*(str*zi+sti*zr)
   60 continue
      if (kode.eq.1) return
      call zsqrt(zr, zi, str, sti)
      ztar = tth*(zr*str-zi*sti)
      ztai = tth*(zr*sti+zi*str)
      call zexp(ztar, ztai, str, sti)
      ptr = str*air - sti*aii
      aii = str*aii + sti*air
      air = ptr
      return
c-----------------------------------------------------------------------
c     case for abs(z).gt.1.0
c-----------------------------------------------------------------------
   70 continue
      fnu = (1.0d0+fid)/3.0d0
c-----------------------------------------------------------------------
c     set parameters related to machine constants.
c     tol is the approximate unit roundoff limited to 1.0d-18.
c     elim is the approximate exponential over- and underflow limit.
c     exp(-elim).lt.exp(-alim)=exp(-elim)/tol    and
c     exp(elim).gt.exp(alim)=exp(elim)*tol       are intervals near
c     underflow and overflow limits where scaled arithmetic is done.
c     rl is the lower boundary of the asymptotic expansion for large z.
c     dig = number of base 10 digits in tol = 10**(-dig).
c-----------------------------------------------------------------------
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
      alaz = log(az)
c-----------------------------------------------------------------------
c     test for proper range
c-----------------------------------------------------------------------
      aa=0.5d0/tol
      bb=i1mach(9)*0.5d0
      aa=min(aa,bb)
      aa=aa**tth
      if (az.gt.aa) go to 260
      aa=sqrt(aa)
      if (az.gt.aa) ierr=3
      call zsqrt(zr, zi, csqr, csqi)
      ztar = tth*(zr*csqr-zi*csqi)
      ztai = tth*(zr*csqi+zi*csqr)
c-----------------------------------------------------------------------
c     re(zta).le.0 when re(z).lt.0, especially when im(z) is small
c-----------------------------------------------------------------------
      iflag = 0
      sfac = 1.0d0
      ak = ztai
      if (zr.ge.0.0d0) go to 80
      bk = ztar
      ck = -abs(bk)
      ztar = ck
      ztai = ak
   80 continue
      if (zi.ne.0.0d0) go to 90
      if (zr.gt.0.0d0) go to 90
      ztar = 0.0d0
      ztai = ak
   90 continue
      aa = ztar
      if (aa.ge.0.0d0 .and. zr.gt.0.0d0) go to 110
      if (kode.eq.2) go to 100
c-----------------------------------------------------------------------
c     overflow test
c-----------------------------------------------------------------------
      if (aa.gt.(-alim)) go to 100
      aa = -aa + 0.25d0*alaz
      iflag = 1
      sfac = tol
      if (aa.gt.elim) go to 270
  100 continue
c-----------------------------------------------------------------------
c     cbknu and cacon return exp(zta)*k(fnu,zta) on kode=2
c-----------------------------------------------------------------------
      mr = 1
      if (zi.lt.0.0d0) mr = -1
      call zacai(ztar, ztai, fnu, kode, mr, 1, cyr, cyi, nn, rl, tol,
     * elim, alim)
      if (nn.lt.0) go to 280
      nz = nz + nn
      go to 130
  110 continue
      if (kode.eq.2) go to 120
c-----------------------------------------------------------------------
c     underflow test
c-----------------------------------------------------------------------
      if (aa.lt.alim) go to 120
      aa = -aa - 0.25d0*alaz
      iflag = 2
      sfac = 1.0d0/tol
      if (aa.lt.(-elim)) go to 210
  120 continue
      call zbknu(ztar, ztai, fnu, kode, 1, cyr, cyi, nz, tol, elim,
     * alim)
  130 continue
      s1r = cyr(1)*coef
      s1i = cyi(1)*coef
      if (iflag.ne.0) go to 150
      if (id.eq.1) go to 140
      air = csqr*s1r - csqi*s1i
      aii = csqr*s1i + csqi*s1r
      return
  140 continue
      air = -(zr*s1r-zi*s1i)
      aii = -(zr*s1i+zi*s1r)
      return
  150 continue
      s1r = s1r*sfac
      s1i = s1i*sfac
      if (id.eq.1) go to 160
      str = s1r*csqr - s1i*csqi
      s1i = s1r*csqi + s1i*csqr
      s1r = str
      air = s1r/sfac
      aii = s1i/sfac
      return
  160 continue
      str = -(s1r*zr-s1i*zi)
      s1i = -(s1r*zi+s1i*zr)
      s1r = str
      air = s1r/sfac
      aii = s1i/sfac
      return
  170 continue
      aa = 1.0d+3*d1mach(1)
      s1r = zeror
      s1i = zeroi
      if (id.eq.1) go to 190
      if (az.le.aa) go to 180
      s1r = c2*zr
      s1i = c2*zi
  180 continue
      air = c1 - s1r
      aii = -s1i
      return
  190 continue
      air = -c2
      aii = 0.0d0
      aa = sqrt(aa)
      if (az.le.aa) go to 200
      s1r = 0.5d0*(zr*zr-zi*zi)
      s1i = zr*zi
  200 continue
      air = air + c1*s1r
      aii = aii + c1*s1i
      return
  210 continue
      nz = 1
      air = zeror
      aii = zeroi
      return
  270 continue
      nz = 0
      ierr=2
      return
  280 continue
      if(nn.eq.(-1)) go to 270
      nz=0
      ierr=5
      return
  260 continue
      ierr=4
      nz=0
      return
      end
