*deck rf
      real function rf (x, y, z, ier)
c***begin prologue  rf
c***purpose  compute the incomplete or complete elliptic integral of the
c            1st kind.  for x, y, and z non-negative and at most one of
c            them zero, rf(x,y,z) = integral from zero to infinity of
c                                -1/2     -1/2     -1/2
c                      (1/2)(t+x)    (t+y)    (t+z)    dt.
c            if x, y or z is zero, the integral is complete.
c***library   slatec
c***category  c14
c***type      single precision (rf-s, drf-d)
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
c   1.     rf
c          evaluate an incomplete (or complete) elliptic integral
c          of the first kind
c          standard fortran function routine
c          single precision version
c          the routine calculates an approximation result to
c          rf(x,y,z) = integral from zero to infinity of
c
c                               -1/2     -1/2     -1/2
c                     (1/2)(t+x)    (t+y)    (t+z)    dt,
c
c          where x, y, and z are nonnegative and at most one of them
c          is zero.  if one of them is zero, the integral is complete.
c          the duplication theorem is iterated until the variables are
c          nearly equal, and the function is then expanded in taylor
c          series to fifth order.
c
c   2.     calling sequence
c          rf( x, y, z, ier )
c
c          parameters on entry
c          values assigned by the calling routine
c
c          x      - single precision, nonnegative variable
c
c          y      - single precision, nonnegative variable
c
c          z      - single precision, nonnegative variable
c
c
c
c          on return     (values assigned by the rf routine)
c
c          rf     - single precision approximation to the integral
c
c          ier    - integer
c
c                   ier = 0 normal and reliable termination of the
c                           routine.  it is assumed that the requested
c                           accuracy has been achieved.
c
c                   ier >  0 abnormal termination of the routine
c
c          x, y, z are unaltered.
c
c
c   3.    error messages
c
c         value of ier assigned by the rf routine
c
c                  value assigned         error message printed
c                  ier = 1                min(x,y,z) .lt. 0.0e0
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
c                     ibm 360/370 series   :   3.0e-78     1.0e+75
c                     cdc 6000/7000 series :   1.0e-292    1.0e+321
c                     univac 1100 series   :   1.0e-37     1.0e+37
c                     cray                 :   2.3e-2466   1.09e+2465
c                     vax 11 series        :   1.5e-38     3.0e+37
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
c              the accuracy of the computed approximation to the inte-
c              gral can be controlled by choosing the value of errtol.
c              truncation of a taylor series after terms of fifth order
c              introduces an error less than the amount shown in the
c              second column of the following table for each value of
c              errtol in the first column.  in addition to the trunca-
c              tion error there will be round-off error, but in prac-
c              tice the total error from both sources is usually less
c              than the amount given in the table.
c
c
c
c
c
c          sample choices:  errtol   relative truncation
c                                    error less than
c                           1.0e-3    3.0e-19
c                           3.0e-3    2.0e-16
c                           1.0e-2    3.0e-13
c                           3.0e-2    2.0e-10
c                           1.0e-1    3.0e-7
c
c
c                    decreasing errtol by a factor of 10 yields six more
c                    decimal digits of accuracy at the expense of one or
c                    two more iterations of the duplication theorem.
c
c *long description:
c
c   rf special comments
c
c
c
c          check by addition theorem: rf(x,x+z,x+w) + rf(y,y+z,y+w)
c          = rf(0,z,w), where x,y,z,w are positive and x * y = z * w.
c
c
c          on input:
c
c          x, y, and z are the variables in the integral rf(x,y,z).
c
c
c          on output:
c
c
c          x, y, and z are unaltered.
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
c   special functions via rf
c
c
c                  legendre form of elliptic integral of 1st kind
c                  ----------------------------------------------
c
c
c                                            2         2   2
c                  f(phi,k) = sin(phi) rf(cos (phi),1-k sin (phi),1)
c
c
c                                 2
c                  k(k) = rf(0,1-k ,1)
c
c                         pi/2     2   2      -1/2
c                       = int  (1-k sin (phi) )   d phi
c                          0
c
c
c
c
c
c                  bulirsch form of elliptic integral of 1st kind
c                  ----------------------------------------------
c
c
c                                         2 2    2
c                  el1(x,kc) = x rf(1,1+kc x ,1+x )
c
c
c
c
c                  lemniscate constant a
c                  ---------------------
c
c
c                       1      4 -1/2
c                  a = int (1-s )    ds = rf(0,1,2) = rf(0,2,1)
c                       0
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
c***routines called  r1mach, xermsg
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
c***end prologue  rf
      character*16 xern3, xern4, xern5, xern6
      integer ier
      real lolim, uplim, epslon, errtol
      real c1, c2, c3, e2, e3, lamda
      real mu, s, x, xn, xndev
      real xnroot, y, yn, yndev, ynroot, z, zn, zndev,
     * znroot
      logical first
      save errtol,lolim,uplim,c1,c2,c3,first
      data first /.true./
c
c***first executable statement  rf
c
      if (first) then
         errtol = (4.0e0*r1mach(3))**(1.0e0/6.0e0)
         lolim  = 5.0e0 * r1mach(1)
         uplim  = r1mach(2)/5.0e0
c
         c1 = 1.0e0/24.0e0
         c2 = 3.0e0/44.0e0
         c3 = 1.0e0/14.0e0
      endif
      first = .false.
c
c         call error handler if necessary.
c
      rf = 0.0e0
      if (min(x,y,z).lt.0.0e0) then
         ier = 1
         write (xern3, '(1pe15.6)') x
         write (xern4, '(1pe15.6)') y
         write (xern5, '(1pe15.6)') z
         call xermsg ('slatec', 'rf',
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
         call xermsg ('slatec', 'rf',
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
         call xermsg ('slatec', 'rf',
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
   30 mu = (xn+yn+zn)/3.0e0
      xndev = 2.0e0 - (mu+xn)/mu
      yndev = 2.0e0 - (mu+yn)/mu
      zndev = 2.0e0 - (mu+zn)/mu
      epslon = max(abs(xndev), abs(yndev), abs(zndev))
      if (epslon.lt.errtol) go to 40
      xnroot =  sqrt(xn)
      ynroot =  sqrt(yn)
      znroot =  sqrt(zn)
      lamda = xnroot*(ynroot+znroot) + ynroot*znroot
      xn = (xn+lamda)*0.250e0
      yn = (yn+lamda)*0.250e0
      zn = (zn+lamda)*0.250e0
      go to 30
c
   40 e2 = xndev*yndev - zndev*zndev
      e3 = xndev*yndev*zndev
      s  = 1.0e0 + (c1*e2-0.10e0-c2*e3)*e2 + c3*e3
      rf = s/sqrt(mu)
c
      return
      end
