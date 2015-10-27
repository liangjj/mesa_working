*deck cbiry
      subroutine cbiry (z, id, kode, bi, ierr)
c***begin prologue  cbiry
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
c         on kode=1, cbiry computes the complex airy function bi(z)
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
c           z      - argument of type complex
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
c           bi     - result of type complex
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
c         flag ierr=3 is triggered where ur=r1mach(4)=unit roundoff.
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
c***routines called  cbinu, i1mach, r1mach
c***revision history  (yymmdd)
c   830501  date written
c   890801  revision date from version 3.2
c   910415  prologue converted to version 4.0 format.  (bab)
c   920128  category corrected.  (wrb)
c   920811  prologue revised.  (dwl)
c***end prologue  cbiry
      complex bi, cone, csq, cy, s1, s2, trm1, trm2, z, zta, z3
      real aa, ad, ak, alim, atrm, az, az3, bb, bk, ck, coef, c1, c2,
     * dig, dk, d1, d2, elim, fid, fmr, fnu, fnul, pi, rl, r1m5, sfac,
     * tol, tth, zi, zr, z3i, z3r, r1mach
      integer id, ierr, k, kode, k1, k2, nz, i1mach
      dimension cy(2)
      data tth, c1, c2, coef, pi /6.66666666666666667e-01,
     * 6.14926627446000736e-01,4.48288357353826359e-01,
     * 5.77350269189625765e-01,3.14159265358979324e+00/
      data  cone / (1.0e0,0.0e0) /
c***first executable statement  cbiry
      ierr = 0
      nz=0
      if (id.lt.0 .or. id.gt.1) ierr=1
      if (kode.lt.1 .or. kode.gt.2) ierr=1
      if (ierr.ne.0) return
      az = abs(z)
      tol = max(r1mach(4),1.0e-18)
      fid = id
      if (az.gt.1.0e0) go to 60
c-----------------------------------------------------------------------
c     power series for abs(z).le.1.
c-----------------------------------------------------------------------
      s1 = cone
      s2 = cone
      if (az.lt.tol) go to 110
      aa = az*az
      if (aa.lt.tol/az) go to 40
      trm1 = cone
      trm2 = cone
      atrm = 1.0e0
      z3 = z*z*z
      az3 = az*aa
      ak = 2.0e0 + fid
      bk = 3.0e0 - fid - fid
      ck = 4.0e0 - fid
      dk = 3.0e0 + fid + fid
      d1 = ak*dk
      d2 = bk*ck
      ad = min(d1,d2)
      ak = 24.0e0 + 9.0e0*fid
      bk = 30.0e0 - 9.0e0*fid
      z3r = real(z3)
      z3i = aimag(z3)
      do 30 k=1,25
        trm1 = trm1*cmplx(z3r/d1,z3i/d1)
        s1 = s1 + trm1
        trm2 = trm2*cmplx(z3r/d2,z3i/d2)
        s2 = s2 + trm2
        atrm = atrm*az3/ad
        d1 = d1 + ak
        d2 = d2 + bk
        ad = min(d1,d2)
        if (atrm.lt.tol*ad) go to 40
        ak = ak + 18.0e0
        bk = bk + 18.0e0
   30 continue
   40 continue
      if (id.eq.1) go to 50
      bi = s1*cmplx(c1,0.0e0) + z*s2*cmplx(c2,0.0e0)
      if (kode.eq.1) return
      zta = z*csqrt(z)*cmplx(tth,0.0e0)
      aa = real(zta)
      aa = -abs(aa)
      bi = bi*cmplx(exp(aa),0.0e0)
      return
   50 continue
      bi = s2*cmplx(c2,0.0e0)
      if (az.gt.tol) bi = bi + z*z*s1*cmplx(c1/(1.0e0+fid),0.0e0)
      if (kode.eq.1) return
      zta = z*csqrt(z)*cmplx(tth,0.0e0)
      aa = real(zta)
      aa = -abs(aa)
      bi = bi*cmplx(exp(aa),0.0e0)
      return
c-----------------------------------------------------------------------
c     case for abs(z).gt.1.0
c-----------------------------------------------------------------------
   60 continue
      fnu = (1.0e0+fid)/3.0e0
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
c-----------------------------------------------------------------------
c     test for range
c-----------------------------------------------------------------------
      aa=0.5e0/tol
      bb=i1mach(9)*0.5e0
      aa=min(aa,bb)
      aa=aa**tth
      if (az.gt.aa) go to 190
      aa=sqrt(aa)
      if (az.gt.aa) ierr=3
      csq=csqrt(z)
      zta=z*csq*cmplx(tth,0.0e0)
c-----------------------------------------------------------------------
c     re(zta).le.0 when re(z).lt.0, especially when im(z) is small
c-----------------------------------------------------------------------
      sfac = 1.0e0
      zi = aimag(z)
      zr = real(z)
      ak = aimag(zta)
      if (zr.ge.0.0e0) go to 70
      bk = real(zta)
      ck = -abs(bk)
      zta = cmplx(ck,ak)
   70 continue
      if (zi.eq.0.0e0 .and. zr.le.0.0e0) zta = cmplx(0.0e0,ak)
      aa = real(zta)
      if (kode.eq.2) go to 80
c-----------------------------------------------------------------------
c     overflow test
c-----------------------------------------------------------------------
      bb = abs(aa)
      if (bb.lt.alim) go to 80
      bb = bb + 0.25e0*alog(az)
      sfac = tol
      if (bb.gt.elim) go to 170
   80 continue
      fmr = 0.0e0
      if (aa.ge.0.0e0 .and. zr.gt.0.0e0) go to 90
      fmr = pi
      if (zi.lt.0.0e0) fmr = -pi
      zta = -zta
   90 continue
c-----------------------------------------------------------------------
c     aa=factor for analytic continuation of i(fnu,zta)
c     kode=2 returns exp(-abs(xzta))*i(fnu,zta) from cbinu
c-----------------------------------------------------------------------
      call cbinu(zta, fnu, kode, 1, cy, nz, rl, fnul, tol, elim, alim)
      if (nz.lt.0) go to 180
      aa = fmr*fnu
      z3 = cmplx(sfac,0.0e0)
      s1 = cy(1)*cmplx(cos(aa),sin(aa))*z3
      fnu = (2.0e0-fid)/3.0e0
      call cbinu(zta, fnu, kode, 2, cy, nz, rl, fnul, tol, elim, alim)
      cy(1) = cy(1)*z3
      cy(2) = cy(2)*z3
c-----------------------------------------------------------------------
c     backward recur one step for orders -1/3 or -2/3
c-----------------------------------------------------------------------
      s2 = cy(1)*cmplx(fnu+fnu,0.0e0)/zta + cy(2)
      aa = fmr*(fnu-1.0e0)
      s1 = (s1+s2*cmplx(cos(aa),sin(aa)))*cmplx(coef,0.0e0)
      if (id.eq.1) go to 100
      s1 = csq*s1
      bi = s1*cmplx(1.0e0/sfac,0.0e0)
      return
  100 continue
      s1 = z*s1
      bi = s1*cmplx(1.0e0/sfac,0.0e0)
      return
  110 continue
      aa = c1*(1.0e0-fid) + fid*c2
      bi = cmplx(aa,0.0e0)
      return
  170 continue
      nz=0
      ierr=2
      return
  180 continue
      if(nz.eq.(-1)) go to 170
      nz=0
      ierr=5
      return
  190 continue
      ierr=4
      nz=0
      return
      end
