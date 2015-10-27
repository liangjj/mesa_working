*deck avint
      subroutine avint (x, y, n, xlo, xup, ans, ierr)
c***begin prologue  avint
c***purpose  integrate a function tabulated at arbitrarily spaced
c            abscissas using overlapping parabolas.
c***library   slatec
c***category  h2a1b2
c***type      single precision (avint-s, davint-d)
c***keywords  integration, quadrature, tabulated data
c***author  jones, r. e., (snla)
c***description
c
c     abstract
c         avint integrates a function tabulated at arbitrarily spaced
c         abscissas.  the limits of integration need not coincide
c         with the tabulated abscissas.
c
c         a method of overlapping parabolas fitted to the data is used
c         provided that there are at least 3 abscissas between the
c         limits of integration.  avint also handles two special cases.
c         if the limits of integration are equal, avint returns a result
c         of zero regardless of the number of tabulated values.
c         if there are only two function values, avint uses the
c         trapezoid rule.
c
c     description of parameters
c         the user must dimension all arrays appearing in the call list
c              x(n), y(n).
c
c         input--
c         x    - real array of abscissas, which must be in increasing
c                order.
c         y    - real array of functional values. i.e., y(i)=func(x(i)).
c         n    - the integer number of function values supplied.
c                n .ge. 2 unless xlo = xup.
c         xlo  - real lower limit of integration.
c         xup  - real upper limit of integration.
c                must have xlo .le. xup.
c
c         output--
c         ans  - computed approximate value of integral
c         ierr - a status code
c              --normal code
c                =1 means the requested integration was performed.
c              --abnormal codes
c                =2 means xup was less than xlo.
c                =3 means the number of x(i) between xlo and xup
c                   (inclusive) was less than 3 and neither of the two
c                   special cases described in the abstract occurred.
c                   no integration was performed.
c                =4 means the restriction x(i+1) .gt. x(i) was violated.
c                =5 means the number n of function values was .lt. 2.
c                ans is set to zero if ierr=2,3,4,or 5.
c
c     avint is documented completely in sc-m-69-335
c     original program from "numerical integration" by davis &
c     rabinowitz.
c     adaptation and modifications for sandia mathematical program
c     library by rondall e. jones.
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
c   900326  removed duplicate information from description section.
c           (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  avint
c
      double precision r3,rp5,sum,syl,syl2,syl3,syu,syu2,syu3,x1,x2,x3
     1,x12,x13,x23,term1,term2,term3,a,b,c,ca,cb,cc
      dimension x(*),y(*)
c***first executable statement  avint
      ierr=1
      ans =0.0
      if (xlo-xup) 3,100,200
    3 if (n.lt.2) go to 215
      do 5 i=2,n
      if (x(i).le.x(i-1)) go to 210
      if (x(i).gt.xup) go to 6
    5 continue
    6 continue
      if (n.ge.3) go to 9
c
c     special n=2 case
      slope = (y(2)-y(1))/(x(2)-x(1))
      fl = y(1) + slope*(xlo-x(1))
      fr = y(2) + slope*(xup-x(2))
      ans = 0.5*(fl+fr)*(xup-xlo)
      return
    9 continue
      if (x(n-2).lt.xlo)  go to 205
      if (x(3).gt.xup)    go to 205
      i = 1
   10 if (x(i).ge.xlo) go to 15
      i = i+1
      go to 10
   15 inlft = i
      i = n
   20 if (x(i).le.xup) go to 25
      i = i-1
      go to 20
   25 inrt = i
      if ((inrt-inlft).lt.2) go to 205
      istart = inlft
      if (inlft.eq.1) istart = 2
      istop  = inrt
      if (inrt.eq.n)  istop  = n-1
c
      r3 = 3.0d0
      rp5= 0.5d0
      sum = 0.0
      syl = xlo
      syl2= syl*syl
      syl3= syl2*syl
c
      do 50 i=istart,istop
      x1 = x(i-1)
      x2 = x(i)
      x3 = x(i+1)
      x12 = x1-x2
      x13 = x1-x3
      x23 = x2-x3
      term1 = dble(y(i-1))/(x12*x13)
      term2 =-dble(y(i)) /(x12*x23)
      term3 = dble(y(i+1))/(x13*x23)
      a = term1+term2+term3
      b = -(x2+x3)*term1 - (x1+x3)*term2 - (x1+x2)*term3
      c = x2*x3*term1 + x1*x3*term2 + x1*x2*term3
      if (i-istart) 30,30,35
   30 ca = a
      cb = b
      cc = c
      go to 40
   35 ca = 0.5*(a+ca)
      cb = 0.5*(b+cb)
      cc = 0.5*(c+cc)
   40 syu = x2
      syu2= syu*syu
      syu3= syu2*syu
      sum = sum + ca*(syu3-syl3)/r3  + cb*rp5*(syu2-syl2) + cc*(syu-syl)
      ca  = a
      cb  = b
      cc  = c
      syl = syu
      syl2= syu2
      syl3= syu3
   50 continue
      syu = xup
      ans = sum + ca*(syu**3-syl3)/r3 + cb*rp5*(syu**2-syl2)
     1  + cc*(syu-syl)
  100 return
  200 ierr=2
      call xermsg ('slatec', 'avint',
     +   'the upper limit of integration was not greater than the ' //
     +   'lower limit.', 4, 1)
      return
  205 ierr=3
      call xermsg ('slatec', 'avint',
     +   'there were less than three function values between the ' //
     +   'limits of integration.', 4, 1)
      return
  210 ierr=4
      call xermsg ('slatec', 'avint',
     +   'the abscissas were not strictly increasing.  must have ' //
     +   'x(i-1) .lt. x(i) for all i.', 4, 1)
      return
  215 ierr=5
      call xermsg ('slatec', 'avint',
     +   'less than two function values were supplied.', 4, 1)
      return
      end
