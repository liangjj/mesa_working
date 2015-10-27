*deck qwgtf
      real function qwgtf (x, omega, p2, p3, p4, integr)
c***begin prologue  qwgtf
c***subsidiary
c***purpose  this function subprogram is used together with the
c            routine qawf and defines the weight function.
c***library   slatec
c***type      single precision (qwgtf-s, dqwgtf-d)
c***keywords  cos or sin in weight function
c***author  piessens, robert
c             applied mathematics and programming division
c             k. u. leuven
c           de doncker, elise
c             applied mathematics and programming division
c             k. u. leuven
c***see also  qk15w
c***routines called  (none)
c***revision history  (yymmdd)
c   810101  date written
c   830518  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c***end prologue  qwgtf
c
      real omega,omx,p2,p3,p4,x
      integer integr
c***first executable statement  qwgtf
      omx = omega*x
      go to(10,20),integr
   10 qwgtf = cos(omx)
      go to 30
   20 qwgtf = sin(omx)
   30 return
      end
