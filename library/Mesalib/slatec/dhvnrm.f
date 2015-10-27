*deck dhvnrm
      double precision function dhvnrm (v, ncomp)
c***begin prologue  dhvnrm
c***subsidiary
c***purpose  subsidiary to ddeabm, ddebdf and dderkf
c***library   slatec
c***type      double precision (hvnrm-s, dhvnrm-d)
c***author  watts, h. a., (snla)
c***description
c
c     compute the maximum norm of the vector v(*) of length ncomp and
c     return the result as dhvnrm
c
c***see also  ddeabm, ddebdf, dderkf
c***routines called  (none)
c***revision history  (yymmdd)
c   820301  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   891024  changed references from dvnorm to dhvnrm.  (wrb)
c   891024  changed routine name from dvnorm to dhvnrm.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910722  updated author section.  (als)
c***end prologue  dhvnrm
c
      integer k, ncomp
      double precision v
      dimension v(*)
c***first executable statement  dhvnrm
      dhvnrm = 0.0d0
      do 10 k = 1, ncomp
         dhvnrm = max(dhvnrm,abs(v(k)))
   10 continue
      return
      end
