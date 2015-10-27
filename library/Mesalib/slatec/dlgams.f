*deck dlgams
      subroutine dlgams (x, dlgam, sgngam)
c***begin prologue  dlgams
c***purpose  compute the logarithm of the absolute value of the gamma
c            function.
c***library   slatec (fnlib)
c***category  c7a
c***type      double precision (algams-s, dlgams-d)
c***keywords  absolute value of the logarithm of the gamma function,
c             fnlib, special functions
c***author  fullerton, w., (lanl)
c***description
c
c dlgams(x,dlgam,sgngam) calculates the double precision natural
c logarithm of the absolute value of the gamma function for
c double precision argument x and stores the result in double
c precision argument dlgam.
c
c***references  (none)
c***routines called  dlngam
c***revision history  (yymmdd)
c   770701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c***end prologue  dlgams
      double precision x, dlgam, sgngam, dlngam
c***first executable statement  dlgams
      dlgam = dlngam(x)
      sgngam = 1.0d0
      if (x.gt.0.d0) return
c
      int = mod (-aint(x), 2.0d0) + 0.1d0
      if (int.eq.0) sgngam = -1.0d0
c
      return
      end
