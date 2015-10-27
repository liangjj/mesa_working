*deck dei
      double precision function dei (x)
c***begin prologue  dei
c***purpose  compute the exponential integral ei(x).
c***library   slatec (fnlib)
c***category  c5
c***type      double precision (ei-s, dei-d)
c***keywords  ei function, exponential integral, fnlib,
c             special functions
c***author  fullerton, w., (lanl)
c***description
c
c dei calculates the double precision exponential integral, ei(x), for
c positive double precision argument x and the cauchy principal value
c for negative x.  if principal values are used everywhere, then, for
c all x,
c
c    ei(x) = -e1(-x)
c or
c    e1(x) = -ei(-x).
c
c***references  (none)
c***routines called  de1
c***revision history  (yymmdd)
c   770701  date written
c   891115  modified prologue description.  (wrb)
c   891115  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c***end prologue  dei
      double precision x, de1
c***first executable statement  dei
      dei = -de1(-x)
c
      return
      end
