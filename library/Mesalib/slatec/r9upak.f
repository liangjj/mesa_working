*deck r9upak
      subroutine r9upak (x, y, n)
c***begin prologue  r9upak
c***purpose  unpack a floating point number x so that x = y*2**n.
c***library   slatec (fnlib)
c***category  a6b
c***type      single precision (r9upak-s, d9upak-d)
c***keywords  fnlib, unpack
c***author  fullerton, w., (lanl)
c***description
c
c   unpack a floating point number x so that x = y*2.0**n, where
c   0.5 .le. abs(y) .lt. 1.0.
c
c***references  (none)
c***routines called  (none)
c***revision history  (yymmdd)
c   780701  date written
c   861211  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c***end prologue  r9upak
c***first executable statement  r9upak
      absx = abs(x)
      n = 0
      if (x.eq.0.0e0) go to 30
c
   10 if (absx.ge.0.5e0) go to 20
      n = n-1
      absx = absx*2.0e0
      go to 10
c
   20 if (absx.lt.1.0e0) go to 30
      n = n+1
      absx = absx*0.5e0
      go to 20
c
   30 y = sign(absx,x)
      return
c
      end
