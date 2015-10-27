      function c_arg(z)
c***begin prologue  c_arg
c***date written   770401   (yymmdd)
c***revision date  820801   (yymmdd)
c***category no.  a4a
c***keywords  argument,complex,complex number,elementary function
c***author  fullerton, w., (lanl)
c***purpose  computes the argument of a complex number.
c***description
c
c c_arg(z) calculates the argument of the complex number z.  note
c that c_arg returns a real result.  if z = x+iy, then c_arg is atan(y/x),
c except when both x and y are zero, in which case the result
c will be zero.
c***references  (none)
c***routines called  (none)
c***end prologue  c_arg
      implicit real*8 (a-h,o-z)
      complex*16 z
c***first executable statement  c_arg
      c_arg = 0.0d0
      a1=real(z)
      a2=imag(z)
      if (a1.ne.0.d0 .or. a2.ne.0.d0) c_arg =
     1  atan2 (a2,a1)
c
      return
      end
