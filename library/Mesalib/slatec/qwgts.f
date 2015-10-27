*deck qwgts
      real function qwgts (x, a, b, alfa, beta, integr)
c***begin prologue  qwgts
c***subsidiary
c***purpose  this function subprogram is used together with the
c            routine qaws and defines the weight function.
c***library   slatec
c***type      single precision (qwgts-s, dqwgts-d)
c***keywords  algebraico-logarithmic, end point singularities,
c             weight function
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
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c***end prologue  qwgts
c
      real a,alfa,b,beta,bmx,x,xma
      integer integr
c***first executable statement  qwgts
      xma = x-a
      bmx = b-x
      qwgts = xma**alfa*bmx**beta
      go to (40,10,20,30),integr
   10 qwgts = qwgts*log(xma)
      go to 40
   20 qwgts = qwgts*log(bmx)
      go to 40
   30 qwgts = qwgts*log(xma)*log(bmx)
   40 return
      end
