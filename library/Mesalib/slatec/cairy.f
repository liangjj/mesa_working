*deck cairy
      subroutine cairy (z, id, kode, ai, nz, ierr)
c***begin prologue  cairy
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
c         on kode=1, cairy computes the complex airy function ai(z)
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
c           z      - argument of type complex
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
c           ai     - result of type complex
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
c***routines called  cacai, cbknu, i1mach, r1mach
c***revision history  (yymmdd)
c   830501  date written
c   890801  revision date from version 3.2
c   910415  prologue converted to version 4.0 format.  (bab)
c   920128  category corrected.  (wrb)
c   920811  prologue revised.  (dwl)
c***end prologue  cairy
      complex ai, cone, csq, cy, s1, s2, trm1, trm2, z, zta, z3
      real aa, ad, ak, alim, atrm, az, az3, bk, ck, coef, c1, c2, dig,
     * dk, d1, d2, elim, fid, fnu, rl, r1m5, sfac, tol, tth, zi, zr,
     * z3i, z3r, r1mach, bb, alaz
      integer id, ierr, iflag, k, kode, k1, k2, mr, nn, nz, i1mach
      dimension cy(1)
      data tth, c1, c2, coef /6.66666666666666667e-01,
     * 3.55028053887817240e-01,2.58819403792806799e-01,
     * 1.83776298473930683e-01/
      data  cone / (1.0e0,0.0e0) /
c***first executable statement  cairy
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
      if (az.lt.tol) go to 160
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
      ai = s1*cmplx(c1,0.0e0) - z*s2*cmplx(c2,0.0e0)
      if (kode.eq.1) return
      zta = z*csqrt(z)*cmplx(tth,0.0e0)
      ai = ai*cexp(zta)
      return
   50 continue
      ai = -s2*cmplx(c2,0.0e0)
      if (az.gt.tol) ai = ai + z*z*s1*cmplx(c1/(1.0e0+fid),0.0e0)
      if (kode.eq.1) return
      zta = z*csqrt(z)*cmplx(tth,0.0e0)
      ai = ai*cexp(zta)
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
      alaz=alog(az)
c-----------------------------------------------------------------------
c     test for range
c-----------------------------------------------------------------------
      aa=0.5e0/tol
      bb=i1mach(9)*0.5e0
      aa=min(aa,bb)
      aa=aa**tth
      if (az.gt.aa) go to 260
      aa=sqrt(aa)
      if (az.gt.aa) ierr=3
      csq=csqrt(z)
      zta=z*csq*cmplx(tth,0.0e0)
c-----------------------------------------------------------------------
c     re(zta).le.0 when re(z).lt.0, especially when im(z) is small
c-----------------------------------------------------------------------
      iflag = 0
      sfac = 1.0e0
      zi = aimag(z)
      zr = real(z)
      ak = aimag(zta)
      if (zr.ge.0.0e0) go to 70
      bk = real(zta)
      ck = -abs(bk)
      zta = cmplx(ck,ak)
   70 continue
      if (zi.ne.0.0e0) go to 80
      if (zr.gt.0.0e0) go to 80
      zta = cmplx(0.0e0,ak)
   80 continue
      aa = real(zta)
      if (aa.ge.0.0e0 .and. zr.gt.0.0e0) go to 100
      if (kode.eq.2) go to 90
c-----------------------------------------------------------------------
c     overflow test
c-----------------------------------------------------------------------
      if (aa.gt.(-alim)) go to 90
      aa = -aa + 0.25e0*alaz
      iflag = 1
      sfac = tol
      if (aa.gt.elim) go to 240
   90 continue
c-----------------------------------------------------------------------
c     cbknu and cacai return exp(zta)*k(fnu,zta) on kode=2
c-----------------------------------------------------------------------
      mr = 1
      if (zi.lt.0.0e0) mr = -1
      call cacai(zta, fnu, kode, mr, 1, cy, nn, rl, tol, elim, alim)
      if (nn.lt.0) go to 250
      nz = nz + nn
      go to 120
  100 continue
      if (kode.eq.2) go to 110
c-----------------------------------------------------------------------
c     underflow test
c-----------------------------------------------------------------------
      if (aa.lt.alim) go to 110
      aa = -aa - 0.25e0*alaz
      iflag = 2
      sfac = 1.0e0/tol
      if (aa.lt.(-elim)) go to 180
  110 continue
      call cbknu(zta, fnu, kode, 1, cy, nz, tol, elim, alim)
  120 continue
      s1 = cy(1)*cmplx(coef,0.0e0)
      if (iflag.ne.0) go to 140
      if (id.eq.1) go to 130
      ai = csq*s1
      return
  130 ai = -z*s1
      return
  140 continue
      s1 = s1*cmplx(sfac,0.0e0)
      if (id.eq.1) go to 150
      s1 = s1*csq
      ai = s1*cmplx(1.0e0/sfac,0.0e0)
      return
  150 continue
      s1 = -s1*z
      ai = s1*cmplx(1.0e0/sfac,0.0e0)
      return
  160 continue
      aa = 1.0e+3*r1mach(1)
      s1 = cmplx(0.0e0,0.0e0)
      if (id.eq.1) go to 170
      if (az.gt.aa) s1 = cmplx(c2,0.0e0)*z
      ai = cmplx(c1,0.0e0) - s1
      return
  170 continue
      ai = -cmplx(c2,0.0e0)
      aa = sqrt(aa)
      if (az.gt.aa) s1 = z*z*cmplx(0.5e0,0.0e0)
      ai = ai + s1*cmplx(c1,0.0e0)
      return
  180 continue
      nz = 1
      ai = cmplx(0.0e0,0.0e0)
      return
  240 continue
      nz = 0
      ierr=2
      return
  250 continue
      if(nn.eq.(-1)) go to 240
      nz=0
      ierr=5
      return
  260 continue
      ierr=4
      nz=0
      return
      end
