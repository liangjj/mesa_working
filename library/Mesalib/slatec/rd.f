*deck rd
      real function rd (x, y, z, ier)
c***begin prologue  rd
c***purpose  compute the incomplete or complete elliptic integral of the
c            2nd kind.  for x and y nonnegative, x+y and z positive,
c             rd(x,y,z) = integral from zero to infinity of
c                                -1/2     -1/2     -3/2
c                      (3/2)(t+x)    (t+y)    (t+z)    dt.
c            if x or y is zero, the integral is complete.
c***library   slatec
c***category  c14
c***type      single precision (rd-s, drd-d)
c***keywords  complete elliptic integral, duplication theorem,
c             incomplete elliptic integral, integral of the second kind,
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
c   1.     rd
c          evaluate an incomplete (or complete) elliptic integral
c          of the second kind
c          standard fortran function routine
c          single precision version
c          the routine calculates an approximation result to
c          rd(x,y,z) = integral from zero to infinity of
c                              -1/2     -1/2     -3/2
c                    (3/2)(t+x)    (t+y)    (t+z)    dt,
c          where x and y are nonnegative, x + y is positive, and z is
c          positive.  if x or y is zero, the integral is complete.
c          the duplication theorem is iterated until the variables are
c          nearly equal, and the function is then expanded in taylor
c          series to fifth order.
c
c   2.     calling sequence
c
c          rd( x, y, z, ier )
c
c          parameters on entry
c          values assigned by the calling routine
c
c          x      - single precision, nonnegative variable
c
c          y      - single precision, nonnegative variable
c
c                   x + y is positive
c
c          z      - real, positive variable
c
c
c
c          on return     (values assigned by the rd routine)
c
c          rd     - real approximation to the integral
c
c
c          ier    - integer
c
c                   ier = 0 normal and reliable termination of the
c                           routine.  it is assumed that the requested
c                           accuracy has been achieved.
c
c                   ier >  0 abnormal termination of the routine
c
c
c          x, y, z are unaltered.
c
c   3.    error messages
c
c         value of ier assigned by the rd routine
c
c                  value assigned         error message printed
c                  ier = 1                min(x,y) .lt. 0.0e0
c                      = 2                min(x + y, z ) .lt. lolim
c                      = 3                max(x,y,z) .gt. uplim
c
c
c   4.     control parameters
c
c                  values of lolim, uplim, and errtol are set by the
c                  routine.
c
c          lolim and uplim determine the valid range of x, y, and z
c
c          lolim  - lower limit of valid arguments
c
c                    not less  than 2 / (machine maximum) ** (2/3).
c
c          uplim  - upper limit of valid arguments
c
c                    not greater than (0.1e0 * errtol / machine
c                    minimum) ** (2/3), where errtol is described below.
c                    in the following table it is assumed that errtol
c                    will never be chosen smaller than 1.0e-5.
c
c
c                    acceptable values for:   lolim      uplim
c                    ibm 360/370 series   :   6.0e-51     1.0e+48
c                    cdc 6000/7000 series :   5.0e-215    2.0e+191
c                    univac 1100 series   :   1.0e-25     2.0e+21
c                    cray                 :   3.0e-1644   1.69e+1640
c                    vax 11 series        :   1.0e-25     4.5e+21
c
c
c          errtol determines the accuracy of the answer
c
c                 the value assigned by the routine will result
c                 in solution precision within 1-2 decimals of
c                 "machine precision".
c
c          errtol    relative error due to truncation is less than
c                    3 * errtol ** 6 / (1-errtol) ** 3/2.
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
c          sample choices:  errtol   relative truncation
c                                    error less than
c                           1.0e-3    4.0e-18
c                           3.0e-3    3.0e-15
c                           1.0e-2    4.0e-12
c                           3.0e-2    3.0e-9
c                           1.0e-1    4.0e-6
c
c
c                    decreasing errtol by a factor of 10 yields six more
c                    decimal digits of accuracy at the expense of one or
c                    two more iterations of the duplication theorem.
c
c *long description:
c
c   rd special comments
c
c
c
c          check: rd(x,y,z) + rd(y,z,x) + rd(z,x,y)
c          = 3 /  sqrt(x * y * z), where x, y, and z are positive.
c
c
c          on input:
c
c          x, y, and z are the variables in the integral rd(x,y,z).
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
c           warning: changes in the program may improve speed at the
c                    expense of robustness.
c
c
c
c    -------------------------------------------------------------------
c
c
c   special functions via rd and rf
c
c
c                  legendre form of elliptic integral of 2nd kind
c                  ----------------------------------------------
c
c
c                                            2         2   2
c                  e(phi,k) = sin(phi) rf(cos (phi),1-k sin (phi),1) -
c
c                     2      3            2         2   2
c                  -(k/3) sin (phi) rd(cos (phi),1-k sin (phi),1)
c
c
c                                 2        2           2
c                  e(k) = rf(0,1-k ,1) - (k/3) rd(3,1-k ,1)
c
c
c                         pi/2     2   2      1/2
c                       = int  (1-k sin (phi) )  d phi
c                          0
c
c
c
c                  bulirsch form of elliptic integral of 2nd kind
c                  ----------------------------------------------
c
c                                              2 2    2
c                  el2(x,kc,a,b) = ax rf(1,1+kc x ,1+x ) +
c
c                                              3         2 2    2
c                                 +(1/3)(b-a) x rd(1,1+kc x ,1+x )
c
c
c
c                  legendre form of alternative elliptic integral of 2nd
c                  -----------------------------------------------------
c                        kind
c                        ----
c
c                            q     2       2   2  -1/2
c                  d(q,k) = int sin p  (1-k sin p)     dp
c                            0
c
c
c
c                                   3          2     2   2
c                  d(q,k) =(1/3)(sin q)  rd(cos q,1-k sin q,1)
c
c
c
c
c
c                  lemniscate constant b
c                  ---------------------
c
c
c
c                       1    2    4 -1/2
c                  b = int  s (1-s )    ds
c                       0
c
c
c                  b =(1/3)rd (0,2,1)
c
c
c
c
c                  heuman's lambda function
c                  ------------------------
c
c
c
c                  (pi/2) lambda0(a,b) =
c
c                                       2                2
c                     = sin(b) (rf(0,cos (a),1)-(1/3) sin (a) *
c
c                               2              2         2       2
c                      *rd(0,cos (a),1)) rf(cos (b),1-cos (a) sin (b),1)
c
c                               2       3            2
c                     -(1/3) cos (a) sin (b) rf(0,cos (a),1) *
c
c                             2         2       2
c                      *rd(cos (b),1-cos (a) sin (b),1)
c
c
c
c                  jacobi zeta function
c                  --------------------
c
c
c                             2                2       2   2
c                  z(b,k) = (k/3) sin(b) rf(cos (b),1-k sin (b),1)
c
c
c                                      2            2
c                             *rd(0,1-k ,1)/rf(0,1-k ,1)
c
c                               2       3          2       2   2
c                            -(k /3) sin (b) rd(cos (b),1-k sin (b),1)
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
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c   900510  modify calls to xermsg to put in standard form.  (rwc)
c   920501  reformatted the references section.  (wrb)
c***end prologue  rd
      character*16 xern3, xern4, xern5, xern6
      integer ier
      real lolim, uplim, epslon, errtol
      real c1, c2, c3, c4, ea, eb, ec, ed, ef, lamda
      real mu, power4, sigma, s1, s2, x, xn, xndev
      real xnroot, y, yn, yndev, ynroot, z, zn, zndev, znroot
      logical first
      save errtol, lolim, uplim, c1, c2, c3, c4, first
      data first /.true./
c
c***first executable statement  rd
      if (first) then
         errtol = (r1mach(3)/3.0e0)**(1.0e0/6.0e0)
         lolim  = 2.0e0/(r1mach(2))**(2.0e0/3.0e0)
         tuplim = r1mach(1)**(1.0e0/3.0e0)
         tuplim = (0.10e0*errtol)**(1.0e0/3.0e0)/tuplim
         uplim  = tuplim**2.0e0
c
         c1 = 3.0e0/14.0e0
         c2 = 1.0e0/6.0e0
         c3 = 9.0e0/22.0e0
         c4 = 3.0e0/26.0e0
      endif
      first = .false.
c
c         call error handler if necessary.
c
      rd = 0.0e0
      if( min(x,y).lt.0.0e0) then
         ier = 1
         write (xern3, '(1pe15.6)') x
         write (xern4, '(1pe15.6)') y
         call xermsg ('slatec', 'rd',
     *      'min(x,y).lt.0 where x = ' // xern3 // ' and y = ' //
     *      xern4, 1, 1)
         return
      endif
c
      if (max(x,y,z).gt.uplim) then
         ier = 3
         write (xern3, '(1pe15.6)') x
         write (xern4, '(1pe15.6)') y
         write (xern5, '(1pe15.6)') z
         write (xern6, '(1pe15.6)') uplim
         call xermsg ('slatec', 'rd',
     *      'max(x,y,z).gt.uplim where x = ' // xern3 // ' y = ' //
     *      xern4 // ' z = ' // xern5 // ' and uplim = ' // xern6,
     *      3, 1)
         return
      endif
c
      if (min(x+y,z).lt.lolim) then
         ier = 2
         write (xern3, '(1pe15.6)') x
         write (xern4, '(1pe15.6)') y
         write (xern5, '(1pe15.6)') z
         write (xern6, '(1pe15.6)') lolim
         call xermsg ('slatec', 'rd',
     *      'min(x+y,z).lt.lolim where x = ' // xern3 // ' y = ' //
     *      xern4 // ' z = ' // xern5 // ' and lolim = ' // xern6,
     *      2, 1)
         return
      endif
c
      ier = 0
      xn = x
      yn = y
      zn = z
      sigma = 0.0e0
      power4 = 1.0e0
c
   30 mu = (xn+yn+3.0e0*zn)*0.20e0
      xndev = (mu-xn)/mu
      yndev = (mu-yn)/mu
      zndev = (mu-zn)/mu
      epslon = max(abs(xndev), abs(yndev), abs(zndev))
      if (epslon.lt.errtol) go to 40
      xnroot = sqrt(xn)
      ynroot = sqrt(yn)
      znroot = sqrt(zn)
      lamda = xnroot*(ynroot+znroot) + ynroot*znroot
      sigma = sigma + power4/(znroot*(zn+lamda))
      power4 = power4*0.250e0
      xn = (xn+lamda)*0.250e0
      yn = (yn+lamda)*0.250e0
      zn = (zn+lamda)*0.250e0
      go to 30
c
   40 ea = xndev*yndev
      eb = zndev*zndev
      ec = ea - eb
      ed = ea - 6.0e0*eb
      ef = ed + ec + ec
      s1 = ed*(-c1+0.250e0*c3*ed-1.50e0*c4*zndev*ef)
      s2 = zndev*(c2*ef+zndev*(-c3*ec+zndev*c4*ea))
      rd = 3.0e0*sigma + power4*(1.0e0+s1+s2)/(mu* sqrt(mu))
c
      return
      end
