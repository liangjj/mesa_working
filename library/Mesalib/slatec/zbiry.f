*deck zbiry
      subroutine zbiry (zr, zi, id, kode, bir, bii, ierr)
c***begin prologue  zbiry
c***purpose  compute the airy function bi(z) or its derivative dbi/dz
c            for complex argument z.  a scaling option is available
c            to help avoid overflow.
c***library   slatec
c***category  c10d
c***type      complex (cbiry-c, zbiry-c)
c***keywords  airy function, bessel function of order one third,
c             bessel function of order two thirds
c***author  amos, d. e., (snl)
c***description
c
c                      ***a double precision routine***
c         on kode=1, zbiry computes the complex airy function bi(z)
c         or its derivative dbi/dz on id=0 or id=1 respectively.
c         on kode=2, a scaling option exp(abs(re(zeta)))*bi(z) or
c         exp(abs(re(zeta)))*dbi/dz is provided to remove the
c         exponential behavior in both the left and right half planes
c         where zeta=(2/3)*z**(3/2).
c
c         the airy functions bi(z) and dbi/dz are analytic in the
c         whole z-plane, and the scaling option does not destroy this
c         property.
c
c         input
c           zr     - double precision real part of argument z
c           zi     - double precision imag part of argument z
c           id     - order of derivative, id=0 or id=1
c           kode   - a parameter to indicate the scaling option
c                    kode=1  returns
c                            bi=bi(z)  on id=0
c                            bi=dbi/dz on id=1
c                            at z=z
c                        =2  returns
c                            bi=exp(abs(re(zeta)))*bi(z)  on id=0
c                            bi=exp(abs(re(zeta)))*dbi/dz on id=1
c                            at z=z where zeta=(2/3)*z**(3/2)
c
c         output
c           bir    - double precision real part of result
c           bii    - double precision imag part of result
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
c         bi(z) and dbi/dz are computed from i bessel functions by
c
c                bi(z) =  c*sqrt(z)*( i(-1/3,zeta) + i(1/3,zeta) )
c               dbi/dz =  c*   z   *( i(-2/3,zeta) + i(2/3,zeta) )
c                    c =  1/sqrt(3)
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
c***routines called  d1mach, i1mach, zabs, zbinu, zdiv, zsqrt
c***revision history  (yymmdd)
c   830501  date written
c   890801  revision date from version 3.2
c   910415  prologue converted to version 4.0 format.  (bab)
c   920128  category corrected.  (wrb)
c   920811  prologue revised.  (dwl)
c   930122  added zsqrt to external statement.  (rwc)
c***end prologue  zbiry
c     complex bi,cone,csq,cy,s1,s2,trm1,trm2,z,zta,z3
      double precision aa, ad, ak, alim, atrm, az, az3, bb, bii, bir,
     * bk, cc, ck, coef, conei, coner, csqi, csqr, cyi, cyr, c1, c2,
     * dig, dk, d1, d2, eaa, elim, fid, fmr, fnu, fnul, pi, rl, r1m5,
     * sfac, sti, str, s1i, s1r, s2i, s2r, tol, trm1i, trm1r, trm2i,
     * trm2r, tth, zi, zr, ztai, ztar, z3i, z3r, d1mach, zabs
      integer id, ierr, k, kode, k1, k2, nz, i1mach
      dimension cyr(2), cyi(2)
      external zabs, zsqrt
      data tth, c1, c2, coef, pi /6.66666666666666667d-01,
     * 6.14926627446000736d-01,4.48288357353826359d-01,
     * 5.77350269189625765d-01,3.14159265358979324d+00/
      data coner, conei /1.0d0,0.0d0/
c***first executable statement  zbiry
      ierr = 0
      nz=0
      if (id.lt.0 .or. id.gt.1) ierr=1
      if (kode.lt.1 .or. kode.gt.2) ierr=1
      if (ierr.ne.0) return
      az = zabs(zr,zi)
      tol = max(d1mach(4),1.0d-18)
      fid = id
      if (az.gt.1.0e0) go to 70
c-----------------------------------------------------------------------
c     power series for abs(z).le.1.
c-----------------------------------------------------------------------
      s1r = coner
      s1i = conei
      s2r = coner
      s2i = conei
      if (az.lt.tol) go to 130
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
      bir = c1*s1r + c2*(zr*s2r-zi*s2i)
      bii = c1*s1i + c2*(zr*s2i+zi*s2r)
      if (kode.eq.1) return
      call zsqrt(zr, zi, str, sti)
      ztar = tth*(zr*str-zi*sti)
      ztai = tth*(zr*sti+zi*str)
      aa = ztar
      aa = -abs(aa)
      eaa = exp(aa)
      bir = bir*eaa
      bii = bii*eaa
      return
   50 continue
      bir = s2r*c2
      bii = s2i*c2
      if (az.le.tol) go to 60
      cc = c1/(1.0d0+fid)
      str = s1r*zr - s1i*zi
      sti = s1r*zi + s1i*zr
      bir = bir + cc*(str*zr-sti*zi)
      bii = bii + cc*(str*zi+sti*zr)
   60 continue
      if (kode.eq.1) return
      call zsqrt(zr, zi, str, sti)
      ztar = tth*(zr*str-zi*sti)
      ztai = tth*(zr*sti+zi*str)
      aa = ztar
      aa = -abs(aa)
      eaa = exp(aa)
      bir = bir*eaa
      bii = bii*eaa
      return
c-----------------------------------------------------------------------
c     case for abs(z).gt.1.0
c-----------------------------------------------------------------------
   70 continue
      fnu = (1.0d0+fid)/3.0d0
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
c     test for range
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
      sfac = 1.0d0
      ak = ztai
      if (zr.ge.0.0d0) go to 80
      bk = ztar
      ck = -abs(bk)
      ztar = ck
      ztai = ak
   80 continue
      if (zi.ne.0.0d0 .or. zr.gt.0.0d0) go to 90
      ztar = 0.0d0
      ztai = ak
   90 continue
      aa = ztar
      if (kode.eq.2) go to 100
c-----------------------------------------------------------------------
c     overflow test
c-----------------------------------------------------------------------
      bb = abs(aa)
      if (bb.lt.alim) go to 100
      bb = bb + 0.25d0*log(az)
      sfac = tol
      if (bb.gt.elim) go to 190
  100 continue
      fmr = 0.0d0
      if (aa.ge.0.0d0 .and. zr.gt.0.0d0) go to 110
      fmr = pi
      if (zi.lt.0.0d0) fmr = -pi
      ztar = -ztar
      ztai = -ztai
  110 continue
c-----------------------------------------------------------------------
c     aa=factor for analytic continuation of i(fnu,zta)
c     kode=2 returns exp(-abs(xzta))*i(fnu,zta) from cbesi
c-----------------------------------------------------------------------
      call zbinu(ztar, ztai, fnu, kode, 1, cyr, cyi, nz, rl, fnul, tol,
     * elim, alim)
      if (nz.lt.0) go to 200
      aa = fmr*fnu
      z3r = sfac
      str = cos(aa)
      sti = sin(aa)
      s1r = (str*cyr(1)-sti*cyi(1))*z3r
      s1i = (str*cyi(1)+sti*cyr(1))*z3r
      fnu = (2.0d0-fid)/3.0d0
      call zbinu(ztar, ztai, fnu, kode, 2, cyr, cyi, nz, rl, fnul, tol,
     * elim, alim)
      cyr(1) = cyr(1)*z3r
      cyi(1) = cyi(1)*z3r
      cyr(2) = cyr(2)*z3r
      cyi(2) = cyi(2)*z3r
c-----------------------------------------------------------------------
c     backward recur one step for orders -1/3 or -2/3
c-----------------------------------------------------------------------
      call zdiv(cyr(1), cyi(1), ztar, ztai, str, sti)
      s2r = (fnu+fnu)*str + cyr(2)
      s2i = (fnu+fnu)*sti + cyi(2)
      aa = fmr*(fnu-1.0d0)
      str = cos(aa)
      sti = sin(aa)
      s1r = coef*(s1r+s2r*str-s2i*sti)
      s1i = coef*(s1i+s2r*sti+s2i*str)
      if (id.eq.1) go to 120
      str = csqr*s1r - csqi*s1i
      s1i = csqr*s1i + csqi*s1r
      s1r = str
      bir = s1r/sfac
      bii = s1i/sfac
      return
  120 continue
      str = zr*s1r - zi*s1i
      s1i = zr*s1i + zi*s1r
      s1r = str
      bir = s1r/sfac
      bii = s1i/sfac
      return
  130 continue
      aa = c1*(1.0d0-fid) + fid*c2
      bir = aa
      bii = 0.0d0
      return
  190 continue
      ierr=2
      nz=0
      return
  200 continue
      if(nz.eq.(-1)) go to 190
      nz=0
      ierr=5
      return
  260 continue
      ierr=4
      nz=0
      return
      end
