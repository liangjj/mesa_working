*deck prvec
      function prvec (m, u, v)
c***begin prologue  prvec
c***subsidiary
c***purpose  subsidiary to bvsup
c***library   slatec
c***type      single precision (prvec-s, dprvec-d)
c***author  watts, h. a., (snla)
c***description
c
c  this subroutine computes the inner product of a vector u
c  with the imaginary product or mate vector corresponding to v
c
c***see also  bvsup
c***routines called  sdot
c***revision history  (yymmdd)
c   750601  date written
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910722  updated author section.  (als)
c***end prologue  prvec
c
      dimension u(*),v(*)
c***first executable statement  prvec
      n=m/2
      np=n+1
      vp=sdot(n,u(1),1,v(np),1)
      prvec=sdot(n,u(np),1,v(1),1) - vp
      return
      end
