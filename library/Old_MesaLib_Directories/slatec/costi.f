*deck costi
      subroutine costi (n, wsave)
c***begin prologue  costi
c***purpose  initialize a work array for cost.
c***library   slatec (fftpack)
c***category  j1a3
c***type      single precision (costi-s)
c***keywords  cosine fourier transform, fftpack
c***author  swarztrauber, p. n., (ncar)
c***description
c
c  subroutine costi initializes the array wsave which is used in
c  subroutine cost.  the prime factorization of n together with
c  a tabulation of the trigonometric functions are computed and
c  stored in wsave.
c
c  input parameter
c
c  n       the length of the sequence to be transformed.  the method
c          is most efficient when n-1 is a product of small primes.
c
c  output parameter
c
c  wsave   a work array which must be dimensioned at least 3*n+15.
c          different wsave arrays are required for different values
c          of n.  the contents of wsave must not be changed between
c          calls of cost.
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
c           (c) changing definition of variable pi by using
c               fortran intrinsic function atan instead of a data
c               statement.
c   881128  modified by dick valent to meet prologue standards.
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  costi
      dimension wsave(*)
c***first executable statement  costi
      if (n .le. 3) return
      pi = 4.*atan(1.)
      nm1 = n-1
      np1 = n+1
      ns2 = n/2
      dt = pi/nm1
      fk = 0.
      do 101 k=2,ns2
         kc = np1-k
         fk = fk+1.
         wsave(k) = 2.*sin(fk*dt)
         wsave(kc) = 2.*cos(fk*dt)
  101 continue
      call rffti (nm1,wsave(n+1))
      return
      end
