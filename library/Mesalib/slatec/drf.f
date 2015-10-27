*deck drf
      double precision function drf (x, y, z, ier)
c***begin prologue  drf
c***purpose  compute the incomplete or complete elliptic integral of the
c            1st kind.  for x, y, and z non-negative and at most one of
c            them zero, rf(x,y,z) = integral from zero to infinity of
c                                -1/2     -1/2     -1/2
c                      (1/2)(t+x)    (t+y)    (t+z)    dt.
c            if x, y or z is zero, the integral is complete.
c***library   slatec
c***category  c14
c***type      double precision (rf-s, drf-d)
c***keywords  complete elliptic integral, duplication theorem,
c             incomplete elliptic integral, integral of the first kind,
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
c   1.     drf
c          evaluate an incomplete (or complete) elliptic integral
c          of the first kind
c          standard fortran function routine
c          double precision version
c          the routine calculates an approximation result to
c          drf(x,y,z) = integral from zero to infinity of
c
c                               -1/2     -1/2     -1/2
c                     (1/2)(t+x)    (t+y)    (t+z)    dt,
c
c          where x, y, and z are nonnegative and at most one of them
c          is zero.  if one of them  is zero, the integral is complete.
c          the duplication theorem is iterated until the variables are
c          nearly equal, and the function is then expanded in taylor
c          series to fifth order.
c
c   2.     calling sequence
c          drf( x, y, z, ier )
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
c
c
c          on return    (values assigned by the drf routine)
c
c          drf     - double precision approximation to the integral
c
c          ier    - integer
c
c                   ier = 0 normal and reliable termination of the
c                           routine. it is assumed that the requested
c                           accuracy has been achieved.
c
c                   ier >  0 abnormal termination of the routine
c
c          x, y, z are unaltered.
c
c
c   3.    error messages
c
c
c         value of ier assigned by the drf routine
c
c                  value assigned         error message printed
c                  ier = 1                min(x,y,z) .lt. 0.0d0
c                      = 2                min(x+y,x+z,y+z) .lt. lolim
c                      = 3                max(x,y,z) .gt. uplim
c
c
c
c   4.     control parameters
c
c                  values of lolim, uplim, and errtol are set by the
c                  routine.
c
c          lolim and uplim determine the valid range of x, y and z
c
c          lolim  - lower limit of valid arguments
c
c                   not less than 5 * (machine minimum).
c
c          uplim  - upper limit of valid arguments
c
c                   not greater than (machine maximum) / 5.
c
c
c                     acceptable values for:   lolim      uplim
c                     ibm 360/370 series   :   3.0d-78     1.0d+75
c                     cdc 6000/7000 series :   1.0d-292    1.0d+321
c                     univac 1100 series   :   1.0d-307    1.0d+307
c                     cray                 :   2.3d-2466   1.09d+2465
c                     vax 11 series        :   1.5d-38     3.0d+37
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
c          errtol - relative error due to truncation is less than
c                   errtol ** 6 / (4 * (1-errtol)  .
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
c
c
c          sample choices:  errtol   relative truncation
c                                    error less than
c                           1.0d-3    3.0d-19
c                           3.0d-3    2.0d-16
c                           1.0d-2    3.0d-13
c                           3.0d-2    2.0d-10
c                           1.0d-1    3.0d-7
c
c
c                    decreasing errtol by a factor of 10 yields six more
c                    decimal digits of accuracy at the expense of one or
c                    two more iterations of the duplication theorem.
c
c *long description:
c
c   drf special comments
c
c
c
c          check by addition theorem: drf(x,x+z,x+w) + drf(y,y+z,y+w)
c          = drf(0,z,w), where x,y,z,w are positive and x * y = z * w.
c
c
c          on input:
c
c          x, y, and z are the variables in the integral drf(x,y,z).
c
c
c          on output:
c
c
c          x, y, z are unaltered.
c
c
c
c          ********************************************************
c
c          warning: changes in the program may improve speed at the
c                   expense of robustness.
c
c
c
c   special double precision functions via drf
c
c
c
c
c                  legendre form of elliptic integral of 1st kind
c
c                  -----------------------------------------
c
c
c
c                                             2         2   2
c                  f(phi,k) = sin(phi) drf(cos (phi),1-k sin (phi),1)
c
c
c                                  2
c                  k(k) = drf(0,1-k ,1)
c
c
c                         pi/2     2   2      -1/2
c                       = int  (1-k sin (phi) )   d phi
c                          0
c
c
c
c                  bulirsch form of elliptic integral of 1st kind
c
c                  -----------------------------------------
c
c
c                                          2 2    2
c                  el1(x,kc) = x drf(1,1+kc x ,1+x )
c
c
c                  lemniscate constant a
c
c                  -----------------------------------------
c
c
c                       1      4 -1/2
c                  a = int (1-s )    ds = drf(0,1,2) = drf(0,2,1)
c                       0
c
c
c
c    -------------------------------------------------------------------
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
c***routines called  d1mach, xermsg
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
c           editorial changes.  (rwc))
c   920501  reformatted the references section.  (wrb)
c***end prologue  drf
      character*16 xern3, xern4, xern5, xern6
      integer ier
      double precision lolim, uplim, epslon, errtol, d1mach
      double precision c1, c2, c3, e2, e3, lamda
      double precision mu, s, x, xn, xndev
      double precision xnroot, y, yn, yndev, ynroot, z, zn, zndev,
     * znroot
      logical first
      save errtol,lolim,uplim,c1,c2,c3,first
      data first /.true./
c
c***first executable statement  drf
c
      if (first) then
         errtol = (4.0d0*d1mach(3))**(1.0d0/6.0d0)
         lolim  = 5.0d0 * d1mach(1)
         uplim  = d1mach(2)/5.0d0
c
         c1 = 1.0d0/24.0d0
         c2 = 3.0d0/44.0d0
         c3 = 1.0d0/14.0d0
      endif
      first = .false.
c
c         call error handler if necessary.
c
      drf = 0.0d0
      if (min(x,y,z).lt.0.0d0) then
         ier = 1
         write (xern3, '(1pe15.6)') x
         write (xern4, '(1pe15.6)') y
         write (xern5, '(1pe15.6)') z
         call xermsg ('slatec', 'drf',
     *      'min(x,y,z).lt.0 where x = ' // xern3 // ' y = ' // xern4 //
     *      ' and z = ' // xern5, 1, 1)
         return
      endif
c
      if (max(x,y,z).gt.uplim) then
         ier = 3
         write (xern3, '(1pe15.6)') x
         write (xern4, '(1pe15.6)') y
         write (xern5, '(1pe15.6)') z
         write (xern6, '(1pe15.6)') uplim
         call xermsg ('slatec', 'drf',
     *      'max(x,y,z).gt.uplim where x = '  // xern3 // ' y = ' //
     *      xern4 // ' z = ' // xern5 // ' and uplim = ' // xern6, 3, 1)
         return
      endif
c
      if (min(x+y,x+z,y+z).lt.lolim) then
         ier = 2
         write (xern3, '(1pe15.6)') x
         write (xern4, '(1pe15.6)') y
         write (xern5, '(1pe15.6)') z
         write (xern6, '(1pe15.6)') lolim
         call xermsg ('slatec', 'drf',
     *      'min(x+y,x+z,y+z).lt.lolim where x = ' // xern3 //
     *      ' y = ' // xern4 // ' z = ' // xern5 // ' and lolim = ' //
     *      xern6, 2, 1)
         return
      endif
c
      ier = 0
      xn = x
      yn = y
      zn = z
c
   30 mu = (xn+yn+zn)/3.0d0
      xndev = 2.0d0 - (mu+xn)/mu
      yndev = 2.0d0 - (mu+yn)/mu
      zndev = 2.0d0 - (mu+zn)/mu
      epslon = max(abs(xndev),abs(yndev),abs(zndev))
      if (epslon.lt.errtol) go to 40
      xnroot = sqrt(xn)
      ynroot = sqrt(yn)
      znroot = sqrt(zn)
      lamda = xnroot*(ynroot+znroot) + ynroot*znroot
      xn = (xn+lamda)*0.250d0
      yn = (yn+lamda)*0.250d0
      zn = (zn+lamda)*0.250d0
      go to 30
c
   40 e2 = xndev*yndev - zndev*zndev
      e3 = xndev*yndev*zndev
      s  = 1.0d0 + (c1*e2-0.10d0-c2*e3)*e2 + c3*e3
      drf = s/sqrt(mu)
c
      return
      end
