*deck dfac
      double precision function dfac (n)
c***begin prologue  dfac
c***purpose  compute the factorial function.
c***library   slatec (fnlib)
c***category  c1
c***type      double precision (fac-s, dfac-d)
c***keywords  factorial, fnlib, special functions
c***author  fullerton, w., (lanl)
c***description
c
c dfac(n) calculates the double precision factorial for integer
c argument n.
c
c***references  (none)
c***routines called  d9lgmc, dgamlm, xermsg
c***revision history  (yymmdd)
c   770601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c***end prologue  dfac
      double precision facn(31), sq2pil, x, xmax, xmin,  d9lgmc
      save facn, sq2pil, nmax
      data facn  (  1) / +.1000000000 0000000000 0000000000 000 d+1    /
      data facn  (  2) / +.1000000000 0000000000 0000000000 000 d+1    /
      data facn  (  3) / +.2000000000 0000000000 0000000000 000 d+1    /
      data facn  (  4) / +.6000000000 0000000000 0000000000 000 d+1    /
      data facn  (  5) / +.2400000000 0000000000 0000000000 000 d+2    /
      data facn  (  6) / +.1200000000 0000000000 0000000000 000 d+3    /
      data facn  (  7) / +.7200000000 0000000000 0000000000 000 d+3    /
      data facn  (  8) / +.5040000000 0000000000 0000000000 000 d+4    /
      data facn  (  9) / +.4032000000 0000000000 0000000000 000 d+5    /
      data facn  ( 10) / +.3628800000 0000000000 0000000000 000 d+6    /
      data facn  ( 11) / +.3628800000 0000000000 0000000000 000 d+7    /
      data facn  ( 12) / +.3991680000 0000000000 0000000000 000 d+8    /
      data facn  ( 13) / +.4790016000 0000000000 0000000000 000 d+9    /
      data facn  ( 14) / +.6227020800 0000000000 0000000000 000 d+10   /
      data facn  ( 15) / +.8717829120 0000000000 0000000000 000 d+11   /
      data facn  ( 16) / +.1307674368 0000000000 0000000000 000 d+13   /
      data facn  ( 17) / +.2092278988 8000000000 0000000000 000 d+14   /
      data facn  ( 18) / +.3556874280 9600000000 0000000000 000 d+15   /
      data facn  ( 19) / +.6402373705 7280000000 0000000000 000 d+16   /
      data facn  ( 20) / +.1216451004 0883200000 0000000000 000 d+18   /
      data facn  ( 21) / +.2432902008 1766400000 0000000000 000 d+19   /
      data facn  ( 22) / +.5109094217 1709440000 0000000000 000 d+20   /
      data facn  ( 23) / +.1124000727 7776076800 0000000000 000 d+22   /
      data facn  ( 24) / +.2585201673 8884976640 0000000000 000 d+23   /
      data facn  ( 25) / +.6204484017 3323943936 0000000000 000 d+24   /
      data facn  ( 26) / +.1551121004 3330985984 0000000000 000 d+26   /
      data facn  ( 27) / +.4032914611 2660563558 4000000000 000 d+27   /
      data facn  ( 28) / +.1088886945 0418352160 7680000000 000 d+29   /
      data facn  ( 29) / +.3048883446 1171386050 1504000000 000 d+30   /
      data facn  ( 30) / +.8841761993 7397019545 4361600000 000 d+31   /
      data facn  ( 31) / +.2652528598 1219105863 6308480000 000 d+33   /
      data sq2pil / 0.9189385332 0467274178 0329736405 62 d0 /
      data nmax / 0 /
c***first executable statement  dfac
      if (nmax.ne.0) go to 10
      call dgamlm (xmin, xmax)
      nmax = xmax - 1.d0
c
 10   if (n .lt. 0) call xermsg ('slatec', 'dfac',
     +   'factorial of negative integer undefined', 1, 2)
c
      if (n.le.30) dfac = facn(n+1)
      if (n.le.30) return
c
      if (n .gt. nmax) call xermsg ('slatec', 'dfac',
     +   'n so big factorial(n) overflows', 2, 2)
c
      x = n + 1
      dfac = exp ((x-0.5d0)*log(x) - x + sq2pil + d9lgmc(x) )
c
      return
      end
