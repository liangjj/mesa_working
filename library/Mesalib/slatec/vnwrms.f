*deck vnwrms
      real function vnwrms (n, v, w)
c***begin prologue  vnwrms
c***subsidiary
c***purpose  subsidiary to debdf
c***library   slatec
c***type      single precision (vnwrms-s, dvnrms-d)
c***author  (unknown)
c***description
c
c   vnwrms computes a weighted root-mean-square vector norm for the
c   integrator package debdf.
c
c***see also  debdf
c***routines called  (none)
c***revision history  (yymmdd)
c   800901  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c***end prologue  vnwrms
c
c
clll. optimize
c-----------------------------------------------------------------------
c this function routine computes the weighted root-mean-square norm
c of the vector of length n contained in the array v, with weights
c contained in the array w of length n..
c   vnwrms = sqrt( (1/n) * sum( v(i)/w(i) )**2 )
c-----------------------------------------------------------------------
      integer n, i
      real v, w, sum
      dimension v(*), w(*)
c***first executable statement  vnwrms
      sum = 0.0e0
      do 10 i = 1,n
 10     sum = sum + (v(i)/w(i))**2
      vnwrms = sqrt(sum/n)
      return
c----------------------- end of function vnwrms ------------------------
      end
