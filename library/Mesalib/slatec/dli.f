*deck dli
      double precision function dli (x)
c***begin prologue  dli
c***purpose  compute the logarithmic integral.
c***library   slatec (fnlib)
c***category  c5
c***type      double precision (ali-s, dli-d)
c***keywords  fnlib, logarithmic integral, special functions
c***author  fullerton, w., (lanl)
c***description
c
c dli(x) calculates the double precision logarithmic integral
c for double precision argument x.
c
c***references  (none)
c***routines called  dei, xermsg
c***revision history  (yymmdd)
c   770701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c***end prologue  dli
      double precision x, dei
c***first executable statement  dli
      if (x .le. 0.d0) call xermsg ('slatec', 'dli',
     +   'log integral undefined for x le 0', 1, 2)
      if (x .eq. 1.d0) call xermsg ('slatec', 'dli',
     +   'log integral undefined for x = 0', 2, 2)
c
      dli = dei (log(x))
c
      return
      end
