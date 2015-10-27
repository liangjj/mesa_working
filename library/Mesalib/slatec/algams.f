*deck algams
      subroutine algams (x, algam, sgngam)
c***begin prologue  algams
c***purpose  compute the logarithm of the absolute value of the gamma
c            function.
c***library   slatec (fnlib)
c***category  c7a
c***type      single precision (algams-s, dlgams-d)
c***keywords  absolute value of the logarithm of the gamma function,
c             fnlib, special functions
c***author  fullerton, w., (lanl)
c***description
c
c evaluates the logarithm of the absolute value of the gamma
c function.
c     x           - input argument
c     algam       - result
c     sgngam      - is set to the sign of gamma(x) and will
c                   be returned at +1.0 or -1.0.
c
c***references  (none)
c***routines called  alngam
c***revision history  (yymmdd)
c   770701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c***end prologue  algams
c***first executable statement  algams
      algam = alngam(x)
      sgngam = 1.0
      if (x.gt.0.0) return
c
      int = mod (-aint(x), 2.0) + 0.1
      if (int.eq.0) sgngam = -1.0
c
      return
      end
