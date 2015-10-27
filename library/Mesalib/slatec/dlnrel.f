*deck dlnrel
      double precision function dlnrel (x)
c***begin prologue  dlnrel
c***purpose  evaluate ln(1+x) accurate in the sense of relative error.
c***library   slatec (fnlib)
c***category  c4b
c***type      double precision (alnrel-s, dlnrel-d, clnrel-c)
c***keywords  elementary functions, fnlib, logarithm
c***author  fullerton, w., (lanl)
c***description
c
c dlnrel(x) calculates the double precision natural logarithm of
c (1.0+x) for double precision argument x.  this routine should
c be used when x is small and accurate to calculate the logarithm
c accurately (in the relative error sense) in the neighborhood
c of 1.0.
c
c series for alnr       on the interval -3.75000e-01 to  3.75000e-01
c                                        with weighted error   6.35e-32
c                                         log weighted error  31.20
c                               significant figures required  30.93
c                                    decimal places required  32.01
c
c***references  (none)
c***routines called  d1mach, dcsevl, initds, xermsg
c***revision history  (yymmdd)
c   770601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c***end prologue  dlnrel
      double precision alnrcs(43), x, xmin,  dcsevl, d1mach
      logical first
      save alnrcs, nlnrel, xmin, first
      data alnrcs(  1) / +.1037869356 2743769800 6862677190 98 d+1     /
      data alnrcs(  2) / -.1336430150 4908918098 7660415531 33 d+0     /
      data alnrcs(  3) / +.1940824913 5520563357 9261993747 50 d-1     /
      data alnrcs(  4) / -.3010755112 7535777690 3765377765 92 d-2     /
      data alnrcs(  5) / +.4869461479 7154850090 4563665091 37 d-3     /
      data alnrcs(  6) / -.8105488189 3175356066 8099430086 22 d-4     /
      data alnrcs(  7) / +.1377884779 9559524782 9382514960 59 d-4     /
      data alnrcs(  8) / -.2380221089 4358970251 3699929149 35 d-5     /
      data alnrcs(  9) / +.4164041621 3865183476 3918599019 89 d-6     /
      data alnrcs( 10) / -.7359582837 8075994984 2668370319 98 d-7     /
      data alnrcs( 11) / +.1311761187 6241674949 1522943450 11 d-7     /
      data alnrcs( 12) / -.2354670931 7742425136 6960923301 75 d-8     /
      data alnrcs( 13) / +.4252277327 6034997775 6380529625 67 d-9     /
      data alnrcs( 14) / -.7719089413 4840796826 1081074933 00 d-10    /
      data alnrcs( 15) / +.1407574648 1359069909 2153564721 91 d-10    /
      data alnrcs( 16) / -.2576907205 8024680627 5370786275 84 d-11    /
      data alnrcs( 17) / +.4734240666 6294421849 1543950059 38 d-12    /
      data alnrcs( 18) / -.8724901267 4742641745 3012632926 75 d-13    /
      data alnrcs( 19) / +.1612461490 2740551465 7398331191 15 d-13    /
      data alnrcs( 20) / -.2987565201 5665773006 7107924168 15 d-14    /
      data alnrcs( 21) / +.5548070120 9082887983 0413216972 79 d-15    /
      data alnrcs( 22) / -.1032461915 8271569595 1413339619 32 d-15    /
      data alnrcs( 23) / +.1925023920 3049851177 8785032448 68 d-16    /
      data alnrcs( 24) / -.3595507346 5265150011 1897078442 66 d-17    /
      data alnrcs( 25) / +.6726454253 7876857892 1945742267 73 d-18    /
      data alnrcs( 26) / -.1260262416 8735219252 0824256375 46 d-18    /
      data alnrcs( 27) / +.2364488440 8606210044 9161589555 19 d-19    /
      data alnrcs( 28) / -.4441937705 0807936898 8783891797 33 d-20    /
      data alnrcs( 29) / +.8354659446 4034259016 2412939946 66 d-21    /
      data alnrcs( 30) / -.1573155941 6479562574 8992535210 66 d-21    /
      data alnrcs( 31) / +.2965312874 0247422686 1543697066 66 d-22    /
      data alnrcs( 32) / -.5594958348 1815947292 1560132266 66 d-23    /
      data alnrcs( 33) / +.1056635426 8835681048 1872841386 66 d-23    /
      data alnrcs( 34) / -.1997248368 0670204548 3149994666 66 d-24    /
      data alnrcs( 35) / +.3778297781 8839361421 0498559999 99 d-25    /
      data alnrcs( 36) / -.7153158688 9081740345 0381653333 33 d-26    /
      data alnrcs( 37) / +.1355248846 3674213646 5020245333 33 d-26    /
      data alnrcs( 38) / -.2569467304 8487567430 0798293333 33 d-27    /
      data alnrcs( 39) / +.4874775606 6216949076 4595199999 99 d-28    /
      data alnrcs( 40) / -.9254211253 0849715321 1323733333 33 d-29    /
      data alnrcs( 41) / +.1757859784 1760239233 2697600000 00 d-29    /
      data alnrcs( 42) / -.3341002667 7731010351 3770666666 66 d-30    /
      data alnrcs( 43) / +.6353393618 0236187354 1802666666 66 d-31    /
      data first /.true./
c***first executable statement  dlnrel
      if (first) then
         nlnrel = initds (alnrcs, 43, 0.1*real(d1mach(3)))
         xmin = -1.0d0 + sqrt(d1mach(4))
      endif
      first = .false.
c
      if (x .le. (-1.d0)) call xermsg ('slatec', 'dlnrel', 'x is le -1'
     +   , 2, 2)
      if (x .lt. xmin) call xermsg ('slatec', 'dlnrel',
     +   'answer lt half precision because x too near -1', 1, 1)
c
      if (abs(x).le.0.375d0) dlnrel = x*(1.d0 -
     1  x*dcsevl (x/.375d0, alnrcs, nlnrel))
c
      if (abs(x).gt.0.375d0) dlnrel = log (1.0d0+x)
c
      return
      end
