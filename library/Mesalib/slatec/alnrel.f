*deck alnrel
      function alnrel (x)
c***begin prologue  alnrel
c***purpose  evaluate ln(1+x) accurate in the sense of relative error.
c***library   slatec (fnlib)
c***category  c4b
c***type      single precision (alnrel-s, dlnrel-d, clnrel-c)
c***keywords  elementary functions, fnlib, logarithm
c***author  fullerton, w., (lanl)
c***description
c
c alnrel(x) evaluates ln(1+x) accurately in the sense of relative
c error when x is very small.  this routine must be used to
c maintain relative error accuracy whenever x is small and
c accurately known.
c
c series for alnr       on the interval -3.75000d-01 to  3.75000d-01
c                                        with weighted error   1.93e-17
c                                         log weighted error  16.72
c                               significant figures required  16.44
c                                    decimal places required  17.40
c
c***references  (none)
c***routines called  csevl, inits, r1mach, xermsg
c***revision history  (yymmdd)
c   770401  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c***end prologue  alnrel
      dimension alnrcs(23)
      logical first
      save alnrcs, nlnrel, xmin, first
      data alnrcs( 1) /   1.0378693562 743770e0 /
      data alnrcs( 2) /   -.1336430150 4908918e0 /
      data alnrcs( 3) /    .0194082491 35520563e0 /
      data alnrcs( 4) /   -.0030107551 12753577e0 /
      data alnrcs( 5) /    .0004869461 47971548e0 /
      data alnrcs( 6) /   -.0000810548 81893175e0 /
      data alnrcs( 7) /    .0000137788 47799559e0 /
      data alnrcs( 8) /   -.0000023802 21089435e0 /
      data alnrcs( 9) /    .0000004164 04162138e0 /
      data alnrcs(10) /   -.0000000735 95828378e0 /
      data alnrcs(11) /    .0000000131 17611876e0 /
      data alnrcs(12) /   -.0000000023 54670931e0 /
      data alnrcs(13) /    .0000000004 25227732e0 /
      data alnrcs(14) /   -.0000000000 77190894e0 /
      data alnrcs(15) /    .0000000000 14075746e0 /
      data alnrcs(16) /   -.0000000000 02576907e0 /
      data alnrcs(17) /    .0000000000 00473424e0 /
      data alnrcs(18) /   -.0000000000 00087249e0 /
      data alnrcs(19) /    .0000000000 00016124e0 /
      data alnrcs(20) /   -.0000000000 00002987e0 /
      data alnrcs(21) /    .0000000000 00000554e0 /
      data alnrcs(22) /   -.0000000000 00000103e0 /
      data alnrcs(23) /    .0000000000 00000019e0 /
      data first /.true./
c***first executable statement  alnrel
      if (first) then
         nlnrel = inits (alnrcs, 23, 0.1*r1mach(3))
         xmin = -1.0 + sqrt(r1mach(4))
      endif
      first = .false.
c
      if (x .le. (-1.0)) call xermsg ('slatec', 'alnrel', 'x is le -1',
     +   2, 2)
      if (x .lt. xmin) call xermsg ('slatec', 'alnrel',
     +   'answer lt half precision because x too near -1', 1, 1)
c
      if (abs(x).le.0.375) alnrel = x*(1. -
     1  x*csevl (x/.375, alnrcs, nlnrel))
      if (abs(x).gt.0.375) alnrel = log (1.0+x)
c
      return
      end
