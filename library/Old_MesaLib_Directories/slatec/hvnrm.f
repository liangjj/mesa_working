*deck hvnrm
      function hvnrm (v, ncomp)
c***begin prologue  hvnrm
c***subsidiary
c***purpose  subsidiary to deabm, debdf and derkf
c***library   slatec
c***type      single precision (hvnrm-s, dhvnrm-d)
c***author  watts, h. a., (snla)
c***description
c
c     compute the maximum norm of the vector v(*) of length ncomp and
c     return the result as hvnrm.
c
c***see also  deabm, debdf, derkf
c***routines called  (none)
c***revision history  (yymmdd)
c   800501  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   891024  changed routine name from vnorm to hvnrm.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910722  updated author section.  (als)
c***end prologue  hvnrm
      dimension v(*)
c***first executable statement  hvnrm
      hvnrm=0.
      do 10 k=1,ncomp
   10   hvnrm=max(hvnrm,abs(v(k)))
      return
      end
