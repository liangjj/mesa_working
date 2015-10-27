*deck dqwgtc
      double precision function dqwgtc (x, c, p2, p3, p4, kp)
c***begin prologue  dqwgtc
c***subsidiary
c***purpose  this function subprogram is used together with the
c            routine dqawc and defines the weight function.
c***library   slatec
c***type      double precision (qwgtc-s, dqwgtc-d)
c***keywords  cauchy principal value, weight function
c***author  piessens, robert
c             applied mathematics and programming division
c             k. u. leuven
c           de doncker, elise
c             applied mathematics and programming division
c             k. u. leuven
c***see also  dqk15w
c***routines called  (none)
c***revision history  (yymmdd)
c   810101  date written
c   830518  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c***end prologue  dqwgtc
c
      double precision c,p2,p3,p4,x
      integer kp
c***first executable statement  dqwgtc
      dqwgtc = 0.1d+01/(x-c)
      return
      end
