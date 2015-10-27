*deck dai
      double precision function dai (x)
c***begin prologue  dai
c***purpose  evaluate the airy function.
c***library   slatec (fnlib)
c***category  c10d
c***type      double precision (ai-s, dai-d)
c***keywords  airy function, fnlib, special functions
c***author  fullerton, w., (lanl)
c***description
c
c dai(x) calculates the double precision airy function for double
c precision argument x.
c
c series for aif        on the interval -1.00000e+00 to  1.00000e+00
c                                        with weighted error   8.37e-33
c                                         log weighted error  32.08
c                               significant figures required  30.87
c                                    decimal places required  32.63
c
c series for aig        on the interval -1.00000e+00 to  1.00000e+00
c                                        with weighted error   7.47e-34
c                                         log weighted error  33.13
c                               significant figures required  31.50
c                                    decimal places required  33.68
c
c***references  (none)
c***routines called  d1mach, d9aimp, daie, dcsevl, initds, xermsg
c***revision history  (yymmdd)
c   770701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   920618  removed space from variable names.  (rwc, wrb)
c***end prologue  dai
      double precision x, aifcs(13), aigcs(13), theta, xm, xmax, x3sml,
     1  z, d1mach, dcsevl, daie, xmaxt
      logical first
      save aifcs, aigcs, naif, naig, x3sml, xmax, first
      data aifcs(  1) / -.3797135849 6669997496 1970894694 14 d-1     /
      data aifcs(  2) / +.5919188853 7263638574 3197280137 77 d-1     /
      data aifcs(  3) / +.9862928057 7279975365 6038910440 60 d-3     /
      data aifcs(  4) / +.6848843819 0765667554 8548301824 12 d-5     /
      data aifcs(  5) / +.2594202596 2194713019 4892790814 03 d-7     /
      data aifcs(  6) / +.6176612774 0813750329 4457496972 36 d-10    /
      data aifcs(  7) / +.1009245417 2466117901 4295562246 01 d-12    /
      data aifcs(  8) / +.1201479251 1179938141 2880332253 33 d-15    /
      data aifcs(  9) / +.1088294558 8716991878 5252954666 66 d-18    /
      data aifcs( 10) / +.7751377219 6684887039 2384000000 00 d-22    /
      data aifcs( 11) / +.4454811203 7175638391 4666666666 66 d-25    /
      data aifcs( 12) / +.2109284523 1692343466 6666666666 66 d-28    /
      data aifcs( 13) / +.8370173591 0741333333 3333333333 33 d-32    /
      data aigcs(  1) / +.1815236558 1161273011 5562099578 64 d-1     /
      data aigcs(  2) / +.2157256316 6010755534 0306388199 68 d-1     /
      data aigcs(  3) / +.2567835698 7483249659 0524280901 33 d-3     /
      data aigcs(  4) / +.1426521411 9792403898 8294969217 21 d-5     /
      data aigcs(  5) / +.4572114920 0180426070 4340975581 91 d-8     /
      data aigcs(  6) / +.9525170843 5647098607 3922788405 92 d-11    /
      data aigcs(  7) / +.1392563460 5771399051 1504206861 90 d-13    /
      data aigcs(  8) / +.1507099914 2762379592 3069911386 66 d-16    /
      data aigcs(  9) / +.1255914831 2567778822 7032053333 33 d-19    /
      data aigcs( 10) / +.8306307377 0821340343 8293333333 33 d-23    /
      data aigcs( 11) / +.4465753849 3718567445 3333333333 33 d-26    /
      data aigcs( 12) / +.1990085503 4518869333 3333333333 33 d-29    /
      data aigcs( 13) / +.7470288525 6533333333 3333333333 33 d-33    /
      data first /.true./
c***first executable statement  dai
      if (first) then
         naif = initds (aifcs, 13, 0.1*real(d1mach(3)))
         naig = initds (aigcs, 13, 0.1*real(d1mach(3)))
c
         x3sml = d1mach(3)**0.3334d0
         xmaxt = (-1.5d0*log(d1mach(1)))**0.6667d0
         xmax = xmaxt - xmaxt*log(xmaxt)/(4.0d0*sqrt(xmaxt)+1.0d0)
     *           - 0.01d0
      endif
      first = .false.
c
      if (x.ge.(-1.d0)) go to 20
      call d9aimp (x, xm, theta)
      dai = xm * cos(theta)
      return
c
 20   if (x.gt.1.0d0) go to 30
      z = 0.0d0
      if (abs(x).gt.x3sml) z = x**3
      dai = 0.375d0 + (dcsevl (z, aifcs, naif) - x*(0.25d0 +
     1  dcsevl (z, aigcs, naig)) )
      return
c
 30   if (x.gt.xmax) go to 40
      dai = daie(x) * exp(-2.0d0*x*sqrt(x)/3.0d0)
      return
c
 40   dai = 0.0d0
      call xermsg ('slatec', 'dai', 'x so big ai underflows', 1, 1)
      return
c
      end
