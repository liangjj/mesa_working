*deck numer
      integer function numer(ngrp)
c**   begin prologue     numer
c**   date written       850601  yymmdd
c**   revision date      yymmdd  yymmdd
c**   keywords           digits
c**   author             martin, richard (lanl)
c**   source             @(#)util.f     3.2   11/9/87
c**   purpose            forms a number from hollerith daata.
c**   description
c     numer is an integer function used as:
c     number=numer(ngrp)
c     ngrp    a vector containing the point group of the molecule.
c     the number is extracted from the second and third elements.
c**   references
c**   routines called    (none)
c**   end prologue       numer
      implicit integer(a-z)
      dimension ngrp(1), n(10)
      data n/1h0, 1h1, 1h2, 1h3, 1h4, 1h5, 1h6, 1h7, 1h8, 1h9/
c
      idigit = ngrp(2)
      jdigit = ngrp(3)
      numer = 0
      n1 = -1
      n2 = -1
      do 20 i=1,10
         if (idigit .eq. n(i)) n1 = i - 1
         if (jdigit .eq. n(i)) n2 = i - 1
 20   continue
      numer = n1
      if (n2 .ne. -1) numer = 10*n1 + n2
      if (numer .lt. 0) numer = 0
c
      return
      end
