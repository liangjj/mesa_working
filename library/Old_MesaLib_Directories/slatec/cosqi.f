*deck cosqi
      subroutine cosqi (n, wsave)
c***begin prologue  cosqi
c***purpose  initialize a work array for cosqf and cosqb.
c***library   slatec (fftpack)
c***category  j1a3
c***type      single precision (cosqi-s)
c***keywords  cosine fourier transform, fftpack
c***author  swarztrauber, p. n., (ncar)
c***description
c
c  subroutine cosqi initializes the work array wsave which is used in
c  both cosqf1 and cosqb1.  the prime factorization of n together with
c  a tabulation of the trigonometric functions are computed and
c  stored in wsave.
c
c  input parameter
c
c  n       the length of the array to be transformed.  the method
c          is most efficient when n is a product of small primes.
c
c  output parameter
c
c  wsave   a work array which must be dimensioned at least 3*n+15.
c          the same work array can be used for both cosqf1 and cosqb1
c          as long as n remains unchanged.  different wsave arrays
c          are required for different values of n.  the contents of
c          wsave must not be changed between calls of cosqf1 or cosqb1.
c
c***references  p. n. swarztrauber, vectorizing the ffts, in parallel
c                 computations (g. rodrigue, ed.), academic press,
c                 1982, pp. 51-83.
c***routines called  rffti
c***revision history  (yymmdd)
c   790601  date written
c   830401  modified to use slatec library source file format.
c   860115  modified by ron boisvert to adhere to fortran 77 by
c           (a) changing dummy array size declarations (1) to (*),
c           (b) changing references to intrinsic function float
c               to real, and
c           (c) changing definition of variable pih by using
c               fortran intrinsic function atan instead of a data
c               statement.
c   881128  modified by dick valent to meet prologue standards.
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  cosqi
      dimension wsave(*)
c***first executable statement  cosqi
      pih = 2.*atan(1.)
      dt = pih/n
      fk = 0.
      do 101 k=1,n
         fk = fk+1.
         wsave(k) = cos(fk*dt)
  101 continue
      call rffti (n,wsave(n+1))
      return
      end
