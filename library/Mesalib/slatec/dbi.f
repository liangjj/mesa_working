*deck dbi
      double precision function dbi (x)
c***begin prologue  dbi
c***purpose  evaluate the bairy function (the airy function of the
c            second kind).
c***library   slatec (fnlib)
c***category  c10d
c***type      double precision (bi-s, dbi-d)
c***keywords  bairy function, fnlib, special functions
c***author  fullerton, w., (lanl)
c***description
c
c dbi(x) calculates the double precision airy function of the
c second kind for double precision argument x.
c
c series for bif        on the interval -1.00000e+00 to  1.00000e+00
c                                        with weighted error   1.45e-32
c                                         log weighted error  31.84
c                               significant figures required  30.85
c                                    decimal places required  32.40
c
c series for big        on the interval -1.00000e+00 to  1.00000e+00
c                                        with weighted error   1.29e-33
c                                         log weighted error  32.89
c                               significant figures required  31.48
c                                    decimal places required  33.45
c
c series for bif2       on the interval  1.00000e+00 to  8.00000e+00
c                                        with weighted error   6.08e-32
c                                         log weighted error  31.22
c                        approx significant figures required  30.8
c                                    decimal places required  31.80
c
c series for big2       on the interval  1.00000e+00 to  8.00000e+00
c                                        with weighted error   4.91e-33
c                                         log weighted error  32.31
c                        approx significant figures required  31.6
c                                    decimal places required  32.90
c
c***references  (none)
c***routines called  d1mach, d9aimp, dbie, dcsevl, initds, xermsg
c***revision history  (yymmdd)
c   770701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c***end prologue  dbi
      double precision x, bifcs(13), bigcs(13), bif2cs(15), big2cs(15),
     1  theta, xm, xmax, x3sml, z,  d1mach, dcsevl, dbie
      logical first
      save bifcs, bigcs, bif2cs, big2cs, nbif, nbig,
     1 nbif2, nbig2, x3sml, xmax, first
      data bifcs(  1) / -.1673021647 1986649483 5374239281 76 d-1     /
      data bifcs(  2) / +.1025233583 4249445611 4263627777 57 d+0     /
      data bifcs(  3) / +.1708309250 7381516539 4296502420 13 d-2     /
      data bifcs(  4) / +.1186254546 7744681179 2164592100 40 d-4     /
      data bifcs(  5) / +.4493290701 7792133694 5318879272 42 d-7     /
      data bifcs(  6) / +.1069820714 3387889067 5677676636 28 d-9     /
      data bifcs(  7) / +.1748064339 9771824706 0105176285 73 d-12    /
      data bifcs(  8) / +.2081023107 1761711025 8818918343 99 d-15    /
      data bifcs(  9) / +.1884981469 5665416509 9279717333 33 d-18    /
      data bifcs( 10) / +.1342577917 3097804625 8826666666 66 d-21    /
      data bifcs( 11) / +.7715959342 9658887893 3333333333 33 d-25    /
      data bifcs( 12) / +.3653387961 7478566399 9999999999 99 d-28    /
      data bifcs( 13) / +.1449756592 7953066666 6666666666 66 d-31    /
      data bigcs(  1) / +.2246622324 8574522283 4682201390 24 d-1     /
      data bigcs(  2) / +.3736477545 3019545441 7275616667 52 d-1     /
      data bigcs(  3) / +.4447621895 7212285696 2152943266 39 d-3     /
      data bigcs(  4) / +.2470807563 6329384245 4945919488 82 d-5     /
      data bigcs(  5) / +.7919135339 5149635134 8624262855 96 d-8     /
      data bigcs(  6) / +.1649807985 1827779880 8878724027 06 d-10    /
      data bigcs(  7) / +.2411990666 4835455909 2475011228 41 d-13    /
      data bigcs(  8) / +.2610373623 6091436985 1847812693 33 d-16    /
      data bigcs(  9) / +.2175308297 7160323853 1237920000 00 d-19    /
      data bigcs( 10) / +.1438694640 0390433219 4837333333 33 d-22    /
      data bigcs( 11) / +.7734912561 2083468629 3333333333 33 d-26    /
      data bigcs( 12) / +.3446929203 3849002666 6666666666 66 d-29    /
      data bigcs( 13) / +.1293891927 3216000000 0000000000 00 d-32    /
      data bif2cs(  1) / +.0998457269 3816041044 6828425799 3 d+0      /
      data bif2cs(  2) / +.4786249778 6300553772 2114673182 31 d+0     /
      data bif2cs(  3) / +.2515521196 0433011771 3244154366 75 d-1     /
      data bif2cs(  4) / +.5820693885 2326456396 5156978722 16 d-3     /
      data bif2cs(  5) / +.7499765964 4377865943 8614573782 17 d-5     /
      data bif2cs(  6) / +.6134602870 3493836681 4030103564 74 d-7     /
      data bif2cs(  7) / +.3462753885 1480632900 4342687333 59 d-9     /
      data bif2cs(  8) / +.1428891008 0270254287 7708467489 31 d-11    /
      data bif2cs(  9) / +.4496270429 8334641895 0564721792 00 d-14    /
      data bif2cs( 10) / +.1114232306 5833011708 4283001066 66 d-16    /
      data bif2cs( 11) / +.2230479106 6175002081 5178666666 66 d-19    /
      data bif2cs( 12) / +.3681577873 6393142842 9226666666 66 d-22    /
      data bif2cs( 13) / +.5096086844 9338261333 3333333333 33 d-25    /
      data bif2cs( 14) / +.6000338692 6288554666 6666666666 66 d-28    /
      data bif2cs( 15) / +.6082749744 6570666666 6666666666 66 d-31    /
      data big2cs(  1) / +.0333056621 4551434046 5176188111 647 d+0    /
      data big2cs(  2) / +.1613092151 2319706761 3287532084 943 d+0    /
      data big2cs(  3) / +.6319007309 6134286912 1615634921 173 d-2    /
      data big2cs(  4) / +.1187904568 1625173638 9780192304 567 d-3    /
      data big2cs(  5) / +.1304534588 6200265614 7116485012 843 d-5    /
      data big2cs(  6) / +.9374125995 5352172954 6809615508 936 d-8    /
      data big2cs(  7) / +.4745801886 7472515378 8510169834 595 d-10   /
      data big2cs(  8) / +.1783107265 0948139980 0065667560 946 d-12   /
      data big2cs(  9) / +.5167591927 8495818037 4276356640 000 d-15   /
      data big2cs( 10) / +.1190045083 8682712512 9496251733 333 d-17   /
      data big2cs( 11) / +.2229828806 6640351727 7063466666 666 d-20   /
      data big2cs( 12) / +.3465519230 2768941972 2666666666 666 d-23   /
      data big2cs( 13) / +.4539263363 2050451413 3333333333 333 d-26   /
      data big2cs( 14) / +.5078849965 1352234666 6666666666 666 d-29   /
      data big2cs( 15) / +.4910206746 9653333333 3333333333 333 d-32   /
      data first /.true./
c***first executable statement  dbi
      if (first) then
         eta = 0.1*real(d1mach(3))
         nbif = initds (bifcs, 13, eta)
         nbig = initds (bigcs, 13, eta)
         nbif2 = initds (bif2cs, 15, eta)
         nbig2 = initds (big2cs, 15, eta)
c
         x3sml = eta**0.3333
         xmax = (1.5*log(d1mach(2)))**0.6666d0
      endif
      first = .false.
c
      if (x.ge.(-1.0d0)) go to 20
      call d9aimp (x, xm, theta)
      dbi = xm * sin(theta)
      return
c
 20   if (x.gt.1.0d0) go to 30
      z = 0.d0
      if (abs(x).gt.x3sml) z = x**3
      dbi = 0.625 + dcsevl (z, bifcs, nbif) + x*(0.4375d0 +
     1  dcsevl (z, bigcs, nbig))
      return
c
 30   if (x.gt.2.0d0) go to 40
      z = (2.0d0*x**3 - 9.0d0)/7.d0
      dbi = 1.125d0 + dcsevl (z, bif2cs, nbif2) + x*(0.625d0 +
     1  dcsevl (z, big2cs, nbig2))
      return
c
 40   if (x .gt. xmax) call xermsg ('slatec', 'dbi',
     +   'x so big that bi overflows', 1, 2)
c
      dbi = dbie(x) * exp(2.0d0*x*sqrt(x)/3.0d0)
      return
c
      end
