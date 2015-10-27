*deck drc
      double precision function drc (x, y, ier)
c***begin prologue  drc
c***purpose  calculate a double precision approximation to
c             drc(x,y) = integral from zero to infinity of
c                              -1/2     -1
c                    (1/2)(t+x)    (t+y)  dt,
c            where x is nonnegative and y is positive.
c***library   slatec
c***category  c14
c***type      double precision (rc-s, drc-d)
c***keywords  duplication theorem, elementary functions,
c             elliptic integral, taylor series
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
c   1.     drc
c          standard fortran function routine
c          double precision version
c          the routine calculates an approximation result to
c          drc(x,y) = integral from zero to infinity of
c
c                              -1/2     -1
c                    (1/2)(t+x)    (t+y)  dt,
c
c          where x is nonnegative and y is positive.  the duplication
c          theorem is iterated until the variables are nearly equal,
c          and the function is then expanded in taylor series to fifth
c          order.  logarithmic, inverse circular, and inverse hyper-
c          bolic functions can be expressed in terms of drc.
c
c   2.     calling sequence
c          drc( x, y, ier )
c
c          parameters on entry
c          values assigned by the calling routine
c
c          x      - double precision, nonnegative variable
c
c          y      - double precision, positive variable
c
c
c
c          on return  (values assigned by the drc routine)
c
c          drc    - double precision approximation to the integral
c
c          ier    - integer to indicate normal or abnormal termination.
c
c                     ier = 0 normal and reliable termination of the
c                             routine.  it is assumed that the requested
c                             accuracy has been achieved.
c
c                     ier > 0 abnormal termination of the routine
c
c          x and y are unaltered.
c
c   3.    error messages
c
c         value of ier assigned by the drc routine
c
c                  value assigned         error message printed
c                  ier = 1                x.lt.0.0d0.or.y.le.0.0d0
c                      = 2                x+y.lt.lolim
c                      = 3                max(x,y) .gt. uplim
c
c   4.     control parameters
c
c                  values of lolim, uplim, and errtol are set by the
c                  routine.
c
c          lolim and uplim determine the valid range of x and y
c
c          lolim  - lower limit of valid arguments
c
c                   not less  than 5 * (machine minimum)  .
c
c          uplim  - upper limit of valid arguments
c
c                   not greater than (machine maximum) / 5 .
c
c
c                     acceptable values for:   lolim       uplim
c                     ibm 360/370 series   :   3.0d-78     1.0d+75
c                     cdc 6000/7000 series :   1.0d-292    1.0d+321
c                     univac 1100 series   :   1.0d-307    1.0d+307
c                     cray                 :   2.3d-2466   1.0d+2465
c                     vax 11 series        :   1.5d-38     3.0d+37
c
c          errtol determines the accuracy of the answer
c
c                 the value assigned by the routine will result
c                 in solution precision within 1-2 decimals of
c                 "machine precision".
c
c
c          errtol  - relative error due to truncation is less than
c                    16 * errtol ** 6 / (1 - 2 * errtol).
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
c          sample choices:  errtol   relative truncation
c                                    error less than
c                           1.0d-3    2.0d-17
c                           3.0d-3    2.0d-14
c                           1.0d-2    2.0d-11
c                           3.0d-2    2.0d-8
c                           1.0d-1    2.0d-5
c
c
c                    decreasing errtol by a factor of 10 yields six more
c                    decimal digits of accuracy at the expense of one or
c                    two more iterations of the duplication theorem.
c
c *long description:
c
c   drc special comments
c
c
c
c
c                  check: drc(x,x+z) + drc(y,y+z) = drc(0,z)
c
c                  where x, y, and z are positive and x * y = z * z
c
c
c          on input:
c
c          x, and y are the variables in the integral drc(x,y).
c
c          on output:
c
c          x and y are unaltered.
c
c
c
c                    drc(0,1/4)=drc(1/16,1/8)=pi=3.14159...
c
c                    drc(9/4,2)=ln(2)
c
c
c
c          ********************************************************
c
c          warning: changes in the program may improve speed at the
c                   expense of robustness.
c
c
c   --------------------------------------------------------------------
c
c   special functions via drc
c
c
c
c                  ln x                x .gt. 0
c
c                                             2
c                  ln(x) = (x-1) drc(((1+x)/2)  , x )
c
c
c   --------------------------------------------------------------------
c
c                  arcsin x            -1 .le. x .le. 1
c
c                                       2
c                  arcsin x = x drc (1-x  ,1 )
c
c   --------------------------------------------------------------------
c
c                  arccos x            0 .le. x .le. 1
c
c
c                                     2       2
c                  arccos x = sqrt(1-x ) drc(x  ,1 )
c
c   --------------------------------------------------------------------
c
c                  arctan x            -inf .lt. x .lt. +inf
c
c                                        2
c                  arctan x = x drc(1,1+x  )
c
c   --------------------------------------------------------------------
c
c                  arccot x            0 .le. x .lt. inf
c
c                                  2   2
c                  arccot x = drc(x  ,x +1 )
c
c   --------------------------------------------------------------------
c
c                  arcsinh x           -inf .lt. x .lt. +inf
c
c                                       2
c                  arcsinh x = x drc(1+x  ,1 )
c
c   --------------------------------------------------------------------
c
c                  arccosh x           x .ge. 1
c
c                                    2         2
c                  arccosh x = sqrt(x -1) drc(x  ,1 )
c
c   --------------------------------------------------------------------
c
c                  arctanh x           -1 .lt. x .lt. 1
c
c                                         2
c                  arctanh x = x drc(1,1-x  )
c
c   --------------------------------------------------------------------
c
c                  arccoth x           x .gt. 1
c
c                                   2   2
c                  arccoth x = drc(x  ,x -1 )
c
c   --------------------------------------------------------------------
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
c***end prologue  drc
      character*16 xern3, xern4, xern5
      integer ier
      double precision c1, c2, errtol, lamda, lolim, d1mach
      double precision mu, s, sn, uplim, x, xn, y, yn
      logical first
      save errtol,lolim,uplim,c1,c2,first
      data first /.true./
c
c***first executable statement  drc
      if (first) then
         errtol = (d1mach(3)/16.0d0)**(1.0d0/6.0d0)
         lolim  = 5.0d0 * d1mach(1)
         uplim  = d1mach(2) / 5.0d0
c
         c1 = 1.0d0/7.0d0
         c2 = 9.0d0/22.0d0
      endif
      first = .false.
c
c         call error handler if necessary.
c
      drc = 0.0d0
      if (x.lt.0.0d0.or.y.le.0.0d0) then
         ier = 1
         write (xern3, '(1pe15.6)') x
         write (xern4, '(1pe15.6)') y
         call xermsg ('slatec', 'drc',
     *      'x.lt.0 .or. y.le.0 where x = ' // xern3 // ' and y = ' //
     *      xern4, 1, 1)
         return
      endif
c
      if (max(x,y).gt.uplim) then
         ier = 3
         write (xern3, '(1pe15.6)') x
         write (xern4, '(1pe15.6)') y
         write (xern5, '(1pe15.6)') uplim
         call xermsg ('slatec', 'drc',
     *      'max(x,y).gt.uplim where x = '  // xern3 // ' y = ' //
     *      xern4 // ' and uplim = ' // xern5, 3, 1)
         return
      endif
c
      if (x+y.lt.lolim) then
         ier = 2
         write (xern3, '(1pe15.6)') x
         write (xern4, '(1pe15.6)') y
         write (xern5, '(1pe15.6)') lolim
         call xermsg ('slatec', 'drc',
     *      'x+y.lt.lolim where x = ' // xern3 // ' y = ' // xern4 //
     *      ' and lolim = ' // xern5, 2, 1)
         return
      endif
c
      ier = 0
      xn = x
      yn = y
c
   30 mu = (xn+yn+yn)/3.0d0
      sn = (yn+mu)/mu - 2.0d0
      if (abs(sn).lt.errtol) go to 40
      lamda = 2.0d0*sqrt(xn)*sqrt(yn) + yn
      xn = (xn+lamda)*0.250d0
      yn = (yn+lamda)*0.250d0
      go to 30
c
   40 s = sn*sn*(0.30d0+sn*(c1+sn*(0.3750d0+sn*c2)))
      drc = (1.0d0+s)/sqrt(mu)
      return
      end
