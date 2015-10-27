*deck davint
      subroutine davint (x, y, n, xlo, xup, ans, ierr)
c***begin prologue  davint
c***purpose  integrate a function tabulated at arbitrarily spaced
c            abscissas using overlapping parabolas.
c***library   slatec
c***category  h2a1b2
c***type      double precision (avint-s, davint-d)
c***keywords  integration, quadrature, tabulated data
c***author  jones, r. e., (snla)
c***description
c
c     abstract
c         davint integrates a function tabulated at arbitrarily spaced
c         abscissas.  the limits of integration need not coincide
c         with the tabulated abscissas.
c
c         a method of overlapping parabolas fitted to the data is used
c         provided that there are at least 3 abscissas between the
c         limits of integration.  davint also handles two special cases.
c         if the limits of integration are equal, davint returns a
c         result of zero regardless of the number of tabulated values.
c         if there are only two function values, davint uses the
c         trapezoid rule.
c
c     description of parameters
c         the user must dimension all arrays appearing in the call list
c              x(n), y(n)
c
c         input--
c      x    - double precision array of abscissas, which must be in
c             increasing order.
c      y    - double precision array of function values. i.e.,
c                y(i)=func(x(i))
c      n    - the integer number of function values supplied.
c                n .ge. 2 unless xlo = xup.
c      xlo  - double precision lower limit of integration
c      xup  - double precision upper limit of integration.  must have
c              xlo.le.xup
c
c         output--
c      ans  - double precision computed approximate value of integral
c      ierr - a status code
c           --normal code
c                =1 means the requested integration was performed.
c           --abnormal codes
c                =2 means xup was less than xlo.
c                =3 means the number of x(i) between xlo and xup
c                   (inclusive) was less than 3 and neither of the two
c                   special cases described in the abstract occurred.
c                   no integration was performed.
c                =4 means the restriction x(i+1).gt.x(i) was violated.
c                =5 means the number n of function values was .lt. 2.
c                   ans is set to zero if ierr=2,3,4,or 5.
c
c    davint is documented completely in sc-m-69-335
c    original program from *numerical integration* by davis & rabinowitz
c    adaptation and modifications by rondall e jones.
c
c***references  r. e. jones, approximate integrator of functions
c                 tabulated at arbitrarily spaced abscissas,
c                 report sc-m-69-335, sandia laboratories, 1969.
c***routines called  xermsg
c***revision history  (yymmdd)
c   690901  date written
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   920501  reformatted the references section.  (wrb)
c***end prologue  davint
c
      integer i, ierr, inlft, inrt, istart, istop, n
      double precision a, ans, b, c, ca, cb, cc, fl, fr, r3, rp5,
     1     slope, sum, syl, syl2, syl3, syu, syu2, syu3, term1, term2,
     2     term3, x, x1, x12, x13, x2, x23, x3, xlo, xup, y
      dimension x(*),y(*)
c     begin block permitting ...exits to 190
c        begin block permitting ...exits to 180
c***first executable statement  davint
            ierr = 1
            ans = 0.0d0
            if (xlo .gt. xup) go to 160
               if (xlo .eq. xup) go to 150
                  if (n .ge. 2) go to 10
                     ierr = 5
                     call xermsg ('slatec', 'davint',
     +                  'less than two function values were supplied.',
     +                  4, 1)
c     ...............exit
                     go to 190
   10             continue
                  do 20 i = 2, n
c        ............exit
                     if (x(i) .le. x(i-1)) go to 180
c                 ...exit
                     if (x(i) .gt. xup) go to 30
   20             continue
   30             continue
                  if (n .ge. 3) go to 40
c
c                    special n=2 case
                     slope = (y(2) - y(1))/(x(2) - x(1))
                     fl = y(1) + slope*(xlo - x(1))
                     fr = y(2) + slope*(xup - x(2))
                     ans = 0.5d0*(fl + fr)*(xup - xlo)
c     ...............exit
                     go to 190
   40             continue
                  if (x(n-2) .ge. xlo) go to 50
                     ierr = 3
                     call xermsg ('slatec', 'davint',
     +                  'there were less than three function values ' //
     +                  'between the limits of integration.', 4, 1)
c     ...............exit
                     go to 190
   50             continue
                  if (x(3) .le. xup) go to 60
                     ierr = 3
                     call xermsg ('slatec', 'davint',
     +                  'there were less than three function values ' //
     +                  'between the limits of integration.', 4, 1)
c     ...............exit
                     go to 190
   60             continue
                  i = 1
   70             if (x(i) .ge. xlo) go to 80
                     i = i + 1
                  go to 70
   80             continue
                  inlft = i
                  i = n
   90             if (x(i) .le. xup) go to 100
                     i = i - 1
                  go to 90
  100             continue
                  inrt = i
                  if ((inrt - inlft) .ge. 2) go to 110
                     ierr = 3
                     call xermsg ('slatec', 'davint',
     +                  'there were less than three function values ' //
     +                  'between the limits of integration.', 4, 1)
c     ...............exit
                     go to 190
  110             continue
                  istart = inlft
                  if (inlft .eq. 1) istart = 2
                  istop = inrt
                  if (inrt .eq. n) istop = n - 1
c
                  r3 = 3.0d0
                  rp5 = 0.5d0
                  sum = 0.0d0
                  syl = xlo
                  syl2 = syl*syl
                  syl3 = syl2*syl
c
                  do 140 i = istart, istop
                     x1 = x(i-1)
                     x2 = x(i)
                     x3 = x(i+1)
                     x12 = x1 - x2
                     x13 = x1 - x3
                     x23 = x2 - x3
                     term1 = y(i-1)/(x12*x13)
                     term2 = -y(i)/(x12*x23)
                     term3 = y(i+1)/(x13*x23)
                     a = term1 + term2 + term3
                     b = -(x2 + x3)*term1 - (x1 + x3)*term2
     1                   - (x1 + x2)*term3
                     c = x2*x3*term1 + x1*x3*term2 + x1*x2*term3
                     if (i .gt. istart) go to 120
                        ca = a
                        cb = b
                        cc = c
                     go to 130
  120                continue
                        ca = 0.5d0*(a + ca)
                        cb = 0.5d0*(b + cb)
                        cc = 0.5d0*(c + cc)
  130                continue
                     syu = x2
                     syu2 = syu*syu
                     syu3 = syu2*syu
                     sum = sum + ca*(syu3 - syl3)/r3
     1                     + cb*rp5*(syu2 - syl2) + cc*(syu - syl)
                     ca = a
                     cb = b
                     cc = c
                     syl = syu
                     syl2 = syu2
                     syl3 = syu3
  140             continue
                  syu = xup
                  ans = sum + ca*(syu**3 - syl3)/r3
     1                  + cb*rp5*(syu**2 - syl2) + cc*(syu - syl)
  150          continue
            go to 170
  160       continue
               ierr = 2
               call xermsg ('slatec', 'davint',
     +            'the upper limit of integration was not greater ' //
     +            'than the lower limit.', 4, 1)
  170       continue
c     ......exit
            go to 190
  180    continue
         ierr = 4
         call xermsg ('slatec', 'davint',
     +      'the abscissas were not strictly increasing.  must have ' //
     +      'x(i-1) .lt. x(i) for all i.', 4, 1)
  190 continue
      return
      end
