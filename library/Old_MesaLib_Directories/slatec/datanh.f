*deck datanh
      double precision function datanh (x)
c***begin prologue  datanh
c***purpose  compute the arc hyperbolic tangent.
c***library   slatec (fnlib)
c***category  c4c
c***type      double precision (atanh-s, datanh-d, catanh-c)
c***keywords  arc hyperbolic tangent, atanh, elementary functions,
c             fnlib, inverse hyperbolic tangent
c***author  fullerton, w., (lanl)
c***description
c
c datanh(x) calculates the double precision arc hyperbolic
c tangent for double precision argument x.
c
c series for atnh       on the interval  0.          to  2.50000e-01
c                                        with weighted error   6.86e-32
c                                         log weighted error  31.16
c                               significant figures required  30.00
c                                    decimal places required  31.88
c
c***references  (none)
c***routines called  d1mach, dcsevl, initds, xermsg
c***revision history  (yymmdd)
c   770601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c***end prologue  datanh
      double precision x, atnhcs(27), dxrel, sqeps, y, dcsevl, d1mach
      logical first
      save atnhcs, nterms, dxrel, sqeps, first
      data atnhcs(  1) / +.9439510239 3195492308 4289221863 3 d-1      /
      data atnhcs(  2) / +.4919843705 5786159472 0003457666 8 d-1      /
      data atnhcs(  3) / +.2102593522 4554327634 7932733175 2 d-2      /
      data atnhcs(  4) / +.1073554449 7761165846 4073104527 6 d-3      /
      data atnhcs(  5) / +.5978267249 2930314786 4278751787 2 d-5      /
      data atnhcs(  6) / +.3505062030 8891348459 6683488620 0 d-6      /
      data atnhcs(  7) / +.2126374343 7653403508 9621931443 1 d-7      /
      data atnhcs(  8) / +.1321694535 7155271921 2980172305 5 d-8      /
      data atnhcs(  9) / +.8365875501 1780703646 2360405295 9 d-10     /
      data atnhcs( 10) / +.5370503749 3110021638 8143458777 2 d-11     /
      data atnhcs( 11) / +.3486659470 1571079229 7124578429 0 d-12     /
      data atnhcs( 12) / +.2284549509 6034330155 2402411972 2 d-13     /
      data atnhcs( 13) / +.1508407105 9447930448 7422906755 8 d-14     /
      data atnhcs( 14) / +.1002418816 8041091261 3699572283 7 d-15     /
      data atnhcs( 15) / +.6698674738 1650695397 1552688298 6 d-17     /
      data atnhcs( 16) / +.4497954546 4949310830 8332762453 3 d-18     /
      data atnhcs( 17) / +.3032954474 2794535416 8236714666 6 d-19     /
      data atnhcs( 18) / +.2052702064 1909368264 6386141866 6 d-20     /
      data atnhcs( 19) / +.1393848977 0538377131 9301461333 3 d-21     /
      data atnhcs( 20) / +.9492580637 2245769719 5895466666 6 d-23     /
      data atnhcs( 21) / +.6481915448 2423076049 8244266666 6 d-24     /
      data atnhcs( 22) / +.4436730205 7236152726 3232000000 0 d-25     /
      data atnhcs( 23) / +.3043465618 5431616389 1200000000 0 d-26     /
      data atnhcs( 24) / +.2091881298 7923934740 4799999999 9 d-27     /
      data atnhcs( 25) / +.1440445411 2340505613 6533333333 3 d-28     /
      data atnhcs( 26) / +.9935374683 1416404650 6666666666 6 d-30     /
      data atnhcs( 27) / +.6863462444 3582600533 3333333333 3 d-31     /
      data first /.true./
c***first executable statement  datanh
      if (first) then
         nterms = initds (atnhcs, 27, 0.1*real(d1mach(3)) )
         dxrel = sqrt(d1mach(4))
         sqeps = sqrt(3.0d0*d1mach(3))
      endif
      first = .false.
c
      y = abs(x)
      if (y .ge. 1.d0) call xermsg ('slatec', 'datanh', 'abs(x) ge 1',
     +   2, 2)
c
      if (1.d0-y .lt. dxrel) call xermsg ('slatec', 'datanh',
     +   'answer lt half precision because abs(x) too near 1', 1, 1)
c
      datanh = x
      if (y.gt.sqeps .and. y.le.0.5d0) datanh = x*(1.0d0 +
     1  dcsevl (8.d0*x*x-1.d0, atnhcs, nterms) )
      if (y.gt.0.5d0) datanh = 0.5d0*log ((1.0d0+x)/(1.0d0-x))
c
      return
      end
