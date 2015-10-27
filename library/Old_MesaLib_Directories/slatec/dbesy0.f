*deck dbesy0
      double precision function dbesy0 (x)
c***begin prologue  dbesy0
c***purpose  compute the bessel function of the second kind of order
c            zero.
c***library   slatec (fnlib)
c***category  c10a1
c***type      double precision (besy0-s, dbesy0-d)
c***keywords  bessel function, fnlib, order zero, second kind,
c             special functions
c***author  fullerton, w., (lanl)
c***description
c
c dbesy0(x) calculates the double precision bessel function of the
c second kind of order zero for double precision argument x.
c
c series for by0        on the interval  0.          to  1.60000e+01
c                                        with weighted error   8.14e-32
c                                         log weighted error  31.09
c                               significant figures required  30.31
c                                    decimal places required  31.73
c
c***references  (none)
c***routines called  d1mach, d9b0mp, dbesj0, dcsevl, initds, xermsg
c***revision history  (yymmdd)
c   770701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c***end prologue  dbesy0
      double precision x, by0cs(19), ampl, theta, twodpi, xsml,
     1  y, d1mach, dcsevl, dbesj0
      logical first
      save by0cs, twodpi, nty0, xsml, first
      data by0cs(  1) / -.1127783939 2865573217 9398054602 8 d-1      /
      data by0cs(  2) / -.1283452375 6042034604 8088453183 8 d+0      /
      data by0cs(  3) / -.1043788479 9794249365 8176227661 8 d+0      /
      data by0cs(  4) / +.2366274918 3969695409 2415926461 3 d-1      /
      data by0cs(  5) / -.2090391647 7004862391 9622395034 2 d-2      /
      data by0cs(  6) / +.1039754539 3905725209 9924657638 1 d-3      /
      data by0cs(  7) / -.3369747162 4239720967 1877534503 7 d-5      /
      data by0cs(  8) / +.7729384267 6706671585 2136721637 1 d-7      /
      data by0cs(  9) / -.1324976772 6642595914 4347606896 4 d-8      /
      data by0cs( 10) / +.1764823261 5404527921 0038936315 8 d-10     /
      data by0cs( 11) / -.1881055071 5801962006 0282301206 9 d-12     /
      data by0cs( 12) / +.1641865485 3661495027 9223718574 9 d-14     /
      data by0cs( 13) / -.1195659438 6046060857 4599100672 0 d-16     /
      data by0cs( 14) / +.7377296297 4401858424 9411242666 6 d-19     /
      data by0cs( 15) / -.3906843476 7104373307 4090666666 6 d-21     /
      data by0cs( 16) / +.1795503664 4361579498 2912000000 0 d-23     /
      data by0cs( 17) / -.7229627125 4480104789 3333333333 3 d-26     /
      data by0cs( 18) / +.2571727931 6351685973 3333333333 3 d-28     /
      data by0cs( 19) / -.8141268814 1636949333 3333333333 3 d-31     /
      data twodpi / 0.6366197723 6758134307 5535053490 057 d0 /
      data first /.true./
c***first executable statement  dbesy0
      if (first) then
         nty0 = initds (by0cs, 19, 0.1*real(d1mach(3)))
         xsml = sqrt(4.0d0*d1mach(3))
      endif
      first = .false.
c
      if (x .le. 0.d0) call xermsg ('slatec', 'dbesy0',
     +   'x is zero or negative', 1, 2)
      if (x.gt.4.0d0) go to 20
c
      y = 0.d0
      if (x.gt.xsml) y = x*x
      dbesy0 = twodpi*log(0.5d0*x)*dbesj0(x) + .375d0 + dcsevl (
     1  .125d0*y-1.d0, by0cs, nty0)
      return
c
 20   call d9b0mp (x, ampl, theta)
      dbesy0 = ampl * sin(theta)
      return
c
      end
