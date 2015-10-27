*deck cgamr
      complex function cgamr (z)
c***begin prologue  cgamr
c***purpose  compute the reciprocal of the gamma function.
c***library   slatec (fnlib)
c***category  c7a
c***type      complex (gamr-s, dgamr-d, cgamr-c)
c***keywords  fnlib, reciprocal gamma function, special functions
c***author  fullerton, w., (lanl)
c***description
c
c cgamr(z) calculates the reciprocal gamma function for complex
c argument z.  this is a preliminary version that is not accurate.
c
c***references  (none)
c***routines called  clngam, xerclr, xgetf, xsetf
c***revision history  (yymmdd)
c   770701  date written
c   861211  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c***end prologue  cgamr
      complex z, clngam
c***first executable statement  cgamr
      cgamr = (0.0, 0.0)
      x = real (z)
      if (x.le.0.0 .and. aint(x).eq.x .and. aimag(z).eq.0.0) return
c
      call xgetf (irold)
      call xsetf (1)
      cgamr = clngam(z)
      call xerclr
      call xsetf (irold)
      cgamr = exp (-cgamr)
c
      return
      end
