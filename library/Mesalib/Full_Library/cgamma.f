      function cgamma(z)
c***begin prologue  cgamma
c***date written   770701   (yymmdd)
c***revision date  820801   (yymmdd)
c***category no.  c7a
c***keywords  complete gamma function,complex,gamma function,
c             special function
c***author  fullerton, w., (lanl)
c***purpose  computes the gamma function of complex argument.
c***description
c
c cgamma(z) calculates the complete gamma function for complex
c argument z.  this is a preliminary version that is portable
c but not accurate.
c***references  (none)
c***routines called  clngam
c***end prologue  cgamma
      implicit real*8(a-h,o-z)
      complex*16 cgamma,  z, clngam
c***first executable statement  cgamma
      cgamma = exp (clngam(z))
c
      return
      end
