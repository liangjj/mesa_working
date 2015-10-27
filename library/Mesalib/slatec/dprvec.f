*deck dprvec
      double precision function dprvec (m, u, v)
c***begin prologue  dprvec
c***subsidiary
c***purpose  subsidiary to dbvsup
c***library   slatec
c***type      double precision (prvec-s, dprvec-d)
c***author  watts, h. a., (snla)
c***description
c
c  this subroutine computes the inner product of a vector u
c  with the imaginary product or mate vector corresponding to v.
c
c***see also  dbvsup
c***routines called  ddot
c***revision history  (yymmdd)
c   750601  date written
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910722  updated author section.  (als)
c***end prologue  dprvec
c
      double precision ddot
      integer m, n, np
      double precision u(*), v(*), vp
c***first executable statement  dprvec
      n = m/2
      np = n + 1
      vp = ddot(n,u(1),1,v(np),1)
      dprvec = ddot(n,u(np),1,v(1),1) - vp
      return
      end
