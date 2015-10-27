*deck drj
      double precision function drj (x, y, z, p, ier)
c***begin prologue  drj
c***purpose  compute the incomplete or complete (x or y or z is zero)
c            elliptic integral of the 3rd kind.  for x, y, and z non-
c            negative, at most one of them zero, and p positive,
c             rj(x,y,z,p) = integral from zero to infinity of
c                              -1/2     -1/2     -1/2     -1
c                    (3/2)(t+x)    (t+y)    (t+z)    (t+p)  dt.
c***library   slatec
c***category  c14
c***type      double precision (rj-s, drj-d)
c***keywords  complete elliptic integral, duplication theorem,
c             incomplete elliptic integral, integral of the third kind,
c             taylor series
c***author  carlson, b. c.
c             ames laboratory-doe
c             iowa state university
c             ames, ia  50011
c           notis, e. m.
c             ames laboratory-doe
c             iowa state university
c             ames, ia  50011
c           pexton, r. l.
c             lawrence livermore national laboratory
c             livermore, ca  94550
c***description
c
c   1.     drj
c          standard fortran function routine
c          double precision version
c          the routine calculates an approximation result to
c          drj(x,y,z,p) = integral from zero to infinity of
c
c                                -1/2     -1/2     -1/2     -1
c                      (3/2)(t+x)    (t+y)    (t+z)    (t+p)  dt,
c
c          where x, y, and z are nonnegative, at most one of them is
c          zero, and p is positive.  if x or y or z is zero, the
c          integral is complete.  the duplication theorem is iterated
c          until the variables are nearly equal, and the function is
c          then expanded in taylor series to fifth order.
c
c
c   2.     calling sequence
c          drj( x, y, z, p, ier )
c
c          parameters on entry
c          values assigned by the calling routine
c
c          x      - double precision, nonnegative variable
c
c          y      - double precision, nonnegative variable
c
c          z      - double precision, nonnegative variable
c
c          p      - double precision, positive variable
c
c
c          on  return    (values assigned by the drj routine)
c
c          drj     - double precision approximation to the integral
c
c          ier    - integer
c
c                   ier = 0 normal and reliable termination of the
c                           routine. it is assumed that the requested
c                           accuracy has been achieved.
c
c                   ier >  0 abnormal termination of the routine
c
c
c          x, y, z, p are unaltered.
c
c
c   3.    error messages
c
c         value of ier assigned by the drj routine
c
c              value assigned         error message printed
c              ier = 1                min(x,y,z) .lt. 0.0d0
c                  = 2                min(x+y,x+z,y+z,p) .lt. lolim
c                  = 3                max(x,y,z,p) .gt. uplim
c
c
c
c   4.     control parameters
c
c                  values of lolim, uplim, and errtol are set by the
c                  routine.
c
c
c          lolim and uplim determine the valid range of x, y, z, and p
c
c          lolim is not less than the cube root of the value
c          of lolim used in the routine for drc.
c
c          uplim is not greater than 0.3 times the cube root of
c          the value of uplim used in the routine for drc.
c
c
c                     acceptable values for:   lolim      uplim
c                     ibm 360/370 series   :   2.0d-26     3.0d+24
c                     cdc 6000/7000 series :   5.0d-98     3.0d+106
c                     univac 1100 series   :   5.0d-103    6.0d+101
c                     cray                 :   1.32d-822   1.4d+821
c                     vax 11 series        :   2.5d-13     9.0d+11
c
c
c
c          errtol determines the accuracy of the answer
c
c                 the value assigned by the routine will result
c                 in solution precision within 1-2 decimals of
c                 "machine precision".
c
c
c
c
c          relative error due to truncation of the series for drj
c          is less than 3 * errtol ** 6 / (1 - errtol) ** 3/2.
c
c
c
c        the accuracy of the computed approximation to the integral
c        can be controlled by choosing the value of errtol.
c        truncation of a taylor series after terms of fifth order
c        introduces an error less than the amount shown in the
c        second column of the following table for each value of
c        errtol in the first column.  in addition to the truncation
c        error there will be round-off error, but in practice the
c        total error from both sources is usually less than the
c        amount given in the table.
c
c
c
c          sample choices:  errtol   relative truncation
c                                    error less than
c                           1.0d-3    4.0d-18
c                           3.0d-3    3.0d-15
c                           1.0d-2    4.0d-12
c                           3.0d-2    3.0d-9
c                           1.0d-1    4.0d-6
c
c                    decreasing errtol by a factor of 10 yields six more
c                    decimal digits of accuracy at the expense of one or
c                    two more iterations of the duplication theorem.
c
c *long description:
c
c   drj special comments
c
c
c     check by addition theorem: drj(x,x+z,x+w,x+p)
c     + drj(y,y+z,y+w,y+p) + (a-b) * drj(a,b,b,a) + 3.0d0 / sqrt(a)
c     = drj(0,z,w,p), where x,y,z,w,p are positive and x * y
c     = z * w,  a = p * p * (x+y+z+w),  b = p * (p+x) * (p+y),
c     and b - a = p * (p-z) * (p-w).  the sum of the third and
c     fourth terms on the left side is 3.0d0 * drc(a,b).
c
c
c          on input:
c
c     x, y, z, and p are the variables in the integral drj(x,y,z,p).
c
c
c          on output:
c
c
c          x, y, z, p are unaltered.
c
c          ********************************************************
c
c          warning: changes in the program may improve speed at the
c                   expense of robustness.
c
c    -------------------------------------------------------------------
c
c
c   special double precision functions via drj and drf
c
c
c                  legendre form of elliptic integral of 3rd kind
c                  -----------------------------------------
c
c
c                          phi         2         -1
c             p(phi,k,n) = int (1+n sin (theta) )   *
c                           0
c
c
c                                  2    2         -1/2
c                             *(1-k  sin (theta) )     d theta
c
c
c                                           2          2   2
c                        = sin (phi) drf(cos (phi), 1-k sin (phi),1)
c
c                                   3             2         2   2
c                         -(n/3) sin (phi) drj(cos (phi),1-k sin (phi),
c
c                                  2
c                         1,1+n sin (phi))
c
c
c
c                  bulirsch form of elliptic integral of 3rd kind
c                  -----------------------------------------
c
c
c                                            2 2    2
c                  el3(x,kc,p) = x drf(1,1+kc x ,1+x ) +
c
c                                            3           2 2    2     2
c                               +(1/3)(1-p) x  drj(1,1+kc x ,1+x ,1+px )
c
c
c                                           2
c                  cel(kc,p,a,b) = a rf(0,kc ,1) +
c
c
c                                                      2
c                                 +(1/3)(b-pa) drj(0,kc ,1,p)
c
c
c                  heuman's lambda function
c                  -----------------------------------------
c
c
c                                2                      2      2    1/2
c                  l(a,b,p) =(cos (a)sin(b)cos(b)/(1-cos (a)sin (b))   )
c
c                                            2         2       2
c                            *(sin(p) drf(cos (p),1-sin (a) sin (p),1)
c
c                                 2       3            2       2
c                            +(sin (a) sin (p)/(3(1-cos (a) sin (b))))
c
c                                    2         2       2
c                            *drj(cos (p),1-sin (a) sin (p),1,1-
c
c                                2       2          2       2
c                            -sin (a) sin (p)/(1-cos (a) sin (b))))
c
c
c
c                  (pi/2) lambda0(a,b) =l(a,b,pi/2) =
c
c                   2                         2       2    -1/2
c              = cos (a)  sin(b) cos(b) (1-cos (a) sin (b))
c
c                           2                  2       2
c                 *drf(0,cos (a),1) + (1/3) sin (a) cos (a)
c
c                                      2       2    -3/2
c                 *sin(b) cos(b) (1-cos (a) sin (b))
c
c                           2         2       2          2       2
c                 *drj(0,cos (a),1,cos (a) cos (b)/(1-cos (a) sin (b)))
c
c
c                  jacobi zeta function
c                  -----------------------------------------
c
c                        2                     2   2    1/2
c             z(b,k) = (k/3) sin(b) cos(b) (1-k sin (b))
c
c
c                                  2      2   2                 2
c                        *drj(0,1-k ,1,1-k sin (b)) / drf (0,1-k ,1)
c
c
c  ---------------------------------------------------------------------
c
c***references  b. c. carlson and e. m. notis, algorithms for incomplete
c                 elliptic integrals, acm transactions on mathematical
c                 software 7, 3 (september 1981), pp. 398-403.
c               b. c. carlson, computing elliptic integrals by
c                 duplication, numerische mathematik 33, (1979),
c                 pp. 1-16.
c               b. c. carlson, elliptic integrals of the first kind,
c                 siam journal of mathematical analysis 8, (1977),
c                 pp. 231-242.
c***routines called  d1mach, drc, xermsg
c***revision history  (yymmdd)
c   790801  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891009  removed unreferenced statement labels.  (wrb)
c   891009  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c   900510  changed calls to xermsg to standard form, and some
c           editorial changes.  (rwc)).
c   920501  reformatted the references section.  (wrb)
c***end prologue  drj
      integer ier
      character*16 xern3, xern4, xern5, xern6, xern7
      double precision alfa, beta, c1, c2, c3, c4, ea, eb, ec, e2, e3
      double precision lolim, uplim, epslon, errtol, d1mach
      double precision lamda, mu, p, pn, pndev
      double precision power4, drc, sigma, s1, s2, s3, x, xn, xndev
      double precision xnroot, y, yn, yndev, ynroot, z, zn, zndev,
     * znroot
      logical first
      save errtol,lolim,uplim,c1,c2,c3,c4,first
      data first /.true./
c
c***first executable statement  drj
      if (first) then
         errtol = (d1mach(3)/3.0d0)**(1.0d0/6.0d0)
         lolim  = (5.0d0 * d1mach(1))**(1.0d0/3.0d0)
         uplim  = 0.30d0*( d1mach(2) / 5.0d0)**(1.0d0/3.0d0)
c
         c1 = 3.0d0/14.0d0
         c2 = 1.0d0/3.0d0
         c3 = 3.0d0/22.0d0
         c4 = 3.0d0/26.0d0
      endif
      first = .false.
c
c         call error handler if necessary.
c
      drj = 0.0d0
      if (min(x,y,z).lt.0.0d0) then
         ier = 1
         write (xern3, '(1pe15.6)') x
         write (xern4, '(1pe15.6)') y
         write (xern5, '(1pe15.6)') z
         call xermsg ('slatec', 'drj',
     *      'min(x,y,z).lt.0 where x = ' // xern3 // ' y = ' // xern4 //
     *      ' and z = ' // xern5, 1, 1)
         return
      endif
c
      if (max(x,y,z,p).gt.uplim) then
         ier = 3
         write (xern3, '(1pe15.6)') x
         write (xern4, '(1pe15.6)') y
         write (xern5, '(1pe15.6)') z
         write (xern6, '(1pe15.6)') p
         write (xern7, '(1pe15.6)') uplim
         call xermsg ('slatec', 'drj',
     *      'max(x,y,z,p).gt.uplim where x = ' // xern3 // ' y = ' //
     *      xern4 // ' z = ' // xern5 // ' p = ' // xern6 //
     *      ' and uplim = ' // xern7, 3, 1)
         return
      endif
c
      if (min(x+y,x+z,y+z,p).lt.lolim) then
         ier = 2
         write (xern3, '(1pe15.6)') x
         write (xern4, '(1pe15.6)') y
         write (xern5, '(1pe15.6)') z
         write (xern6, '(1pe15.6)') p
         write (xern7, '(1pe15.6)') lolim
         call xermsg ('slatec', 'rj',
     *      'min(x+y,x+z,y+z,p).lt.lolim where x = ' // xern3 //
     *      ' y = ' // xern4 // ' z = '  // xern5 // ' p = ' // xern6 //
     *      ' and lolim = ', 2, 1)
         return
      endif
c
      ier = 0
      xn = x
      yn = y
      zn = z
      pn = p
      sigma = 0.0d0
      power4 = 1.0d0
c
   30 mu = (xn+yn+zn+pn+pn)*0.20d0
      xndev = (mu-xn)/mu
      yndev = (mu-yn)/mu
      zndev = (mu-zn)/mu
      pndev = (mu-pn)/mu
      epslon = max(abs(xndev), abs(yndev), abs(zndev), abs(pndev))
      if (epslon.lt.errtol) go to 40
      xnroot =  sqrt(xn)
      ynroot =  sqrt(yn)
      znroot =  sqrt(zn)
      lamda = xnroot*(ynroot+znroot) + ynroot*znroot
      alfa = pn*(xnroot+ynroot+znroot) + xnroot*ynroot*znroot
      alfa = alfa*alfa
      beta = pn*(pn+lamda)*(pn+lamda)
      sigma = sigma + power4*drc(alfa,beta,ier)
      power4 = power4*0.250d0
      xn = (xn+lamda)*0.250d0
      yn = (yn+lamda)*0.250d0
      zn = (zn+lamda)*0.250d0
      pn = (pn+lamda)*0.250d0
      go to 30
c
   40 ea = xndev*(yndev+zndev) + yndev*zndev
      eb = xndev*yndev*zndev
      ec = pndev*pndev
      e2 = ea - 3.0d0*ec
      e3 = eb + 2.0d0*pndev*(ea-ec)
      s1 = 1.0d0 + e2*(-c1+0.750d0*c3*e2-1.50d0*c4*e3)
      s2 = eb*(0.50d0*c2+pndev*(-c3-c3+pndev*c4))
      s3 = pndev*ea*(c2-pndev*c3) - c2*pndev*ec
      drj = 3.0d0*sigma + power4*(s1+s2+s3)/(mu* sqrt(mu))
      return
      end
