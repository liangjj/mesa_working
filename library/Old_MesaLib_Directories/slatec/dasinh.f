*deck dasinh
      double precision function dasinh (x)
c***begin prologue  dasinh
c***purpose  compute the arc hyperbolic sine.
c***library   slatec (fnlib)
c***category  c4c
c***type      double precision (asinh-s, dasinh-d, casinh-c)
c***keywords  arc hyperbolic sine, asinh, elementary functions, fnlib,
c             inverse hyperbolic sine
c***author  fullerton, w., (lanl)
c***description
c
c dasinh(x) calculates the double precision arc hyperbolic
c sine for double precision argument x.
c
c***references  (none)
c***routines called  d1mach, dcsevl, initds
c***revision history  (yymmdd)
c   770601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c***end prologue  dasinh
      double precision x, asnhcs(39), aln2, sqeps, xmax, y,
     1  dcsevl, d1mach
      logical first
      save asnhcs, aln2, nterms, xmax, sqeps, first
      data asnhcs(  1) / -.1282003991 1738186343 3721273592 68 d+0     /
      data asnhcs(  2) / -.5881176118 9951767565 2117571383 62 d-1     /
      data asnhcs(  3) / +.4727465432 2124815640 7252497560 29 d-2     /
      data asnhcs(  4) / -.4938363162 6536172101 3601747902 73 d-3     /
      data asnhcs(  5) / +.5850620705 8557412287 4948352593 21 d-4     /
      data asnhcs(  6) / -.7466998328 9313681354 7550692171 88 d-5     /
      data asnhcs(  7) / +.1001169358 3558199265 9661920158 12 d-5     /
      data asnhcs(  8) / -.1390354385 8708333608 6164722588 86 d-6     /
      data asnhcs(  9) / +.1982316948 3172793547 3173602371 48 d-7     /
      data asnhcs( 10) / -.2884746841 7848843612 7472728003 17 d-8     /
      data asnhcs( 11) / +.4267296546 7159937953 4575149959 07 d-9     /
      data asnhcs( 12) / -.6397608465 4366357868 7526323096 81 d-10    /
      data asnhcs( 13) / +.9699168608 9064704147 8782931311 79 d-11    /
      data asnhcs( 14) / -.1484427697 2043770830 2466583656 96 d-11    /
      data asnhcs( 15) / +.2290373793 9027447988 0401843789 83 d-12    /
      data asnhcs( 16) / -.3558839513 2732645159 9789426513 10 d-13    /
      data asnhcs( 17) / +.5563969408 0056789953 3745390885 54 d-14    /
      data asnhcs( 18) / -.8746250959 9624678045 6665935201 62 d-15    /
      data asnhcs( 19) / +.1381524884 4526692155 8688022981 29 d-15    /
      data asnhcs( 20) / -.2191668828 2900363984 9551422641 49 d-16    /
      data asnhcs( 21) / +.3490465852 4827565638 3139237068 80 d-17    /
      data asnhcs( 22) / -.5578578840 0895742439 6301570321 06 d-18    /
      data asnhcs( 23) / +.8944514661 7134012551 0508827989 33 d-19    /
      data asnhcs( 24) / -.1438342634 6571317305 5518452394 66 d-19    /
      data asnhcs( 25) / +.2319181187 2169963036 3261446826 66 d-20    /
      data asnhcs( 26) / -.3748700795 3314343674 5706045439 99 d-21    /
      data asnhcs( 27) / +.6073210982 2064279404 5492428800 00 d-22    /
      data asnhcs( 28) / -.9859940276 4633583177 3701734400 00 d-23    /
      data asnhcs( 29) / +.1603921745 2788496315 2326382933 33 d-23    /
      data asnhcs( 30) / -.2613884735 0287686596 7161343999 99 d-24    /
      data asnhcs( 31) / +.4267084960 6857390833 3581653333 33 d-25    /
      data asnhcs( 32) / -.6977021703 9185243299 7307733333 33 d-26    /
      data asnhcs( 33) / +.1142508833 6806858659 8126933333 33 d-26    /
      data asnhcs( 34) / -.1873529207 8860968933 0210133333 33 d-27    /
      data asnhcs( 35) / +.3076358441 4464922794 0659200000 00 d-28    /
      data asnhcs( 36) / -.5057736403 1639824787 0463999999 99 d-29    /
      data asnhcs( 37) / +.8325075471 2689142224 2133333333 33 d-30    /
      data asnhcs( 38) / -.1371845728 2501044163 9253333333 33 d-30    /
      data asnhcs( 39) / +.2262986842 6552784104 1066666666 66 d-31    /
      data aln2 / 0.6931471805 5994530941 7232121458 18d0 /
      data first /.true./
c***first executable statement  dasinh
      if (first) then
         nterms = initds (asnhcs, 39, 0.1*real(d1mach(3)) )
         sqeps = sqrt(d1mach(3))
         xmax = 1.0d0/sqeps
      endif
      first = .false.
c
      y = abs(x)
      if (y.gt.1.0d0) go to 20
c
      dasinh = x
      if (y.gt.sqeps) dasinh = x*(1.0d0 + dcsevl (2.d0*x*x-1.d0,
     1  asnhcs, nterms) )
      return
 20   if (y.lt.xmax) dasinh = log (y+sqrt(y*y+1.d0))
      if (y.ge.xmax) dasinh = aln2 + log(y)
      dasinh = sign (dasinh, x)
      return
c
      end
