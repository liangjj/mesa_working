*deck dvnrms
      double precision function dvnrms (n, v, w)
c***begin prologue  dvnrms
c***subsidiary
c***purpose  subsidiary to ddebdf
c***library   slatec
c***type      double precision (vnwrms-s, dvnrms-d)
c***author  (unknown)
c***description
c
c   dvnrms computes a weighted root-mean-square vector norm for the
c   integrator package ddebdf.
c
c***see also  ddebdf
c***routines called  (none)
c***revision history  (yymmdd)
c   820301  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   890911  removed unnecessary intrinsics.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c***end prologue  dvnrms
      integer i, n
      double precision sum, v, w
      dimension v(*),w(*)
c***first executable statement  dvnrms
      sum = 0.0d0
      do 10 i = 1, n
         sum = sum + (v(i)/w(i))**2
   10 continue
      dvnrms = sqrt(sum/n)
      return
c     ----------------------- end of function dvnrms
c     ------------------------
      end
