*deck sinti
      subroutine sinti (n, wsave)
c***begin prologue  sinti
c***purpose  initialize a work array for sint.
c***library   slatec (fftpack)
c***category  j1a3
c***type      single precision (sinti-s)
c***keywords  fftpack, fourier transform
c***author  swarztrauber, p. n., (ncar)
c***description
c
c  subroutine sinti initializes the array wsave which is used in
c  subroutine sint.  the prime factorization of n together with
c  a tabulation of the trigonometric functions are computed and
c  stored in wsave.
c
c  input parameter
c
c  n       the length of the sequence to be transformed.  the method
c          is most efficient when n+1 is a product of small primes.
c
c  output parameter
c
c  wsave   a work array with at least int(3.5*n+16) locations.
c          different wsave arrays are required for different values
c          of n.  the contents of wsave must not be changed between
c          calls of sint.
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
c***end prologue  sinti
      dimension wsave(*)
c***first executable statement  sinti
      if (n .le. 1) return
      pi = 4.*atan(1.)
      np1 = n+1
      ns2 = n/2
      dt = pi/np1
      ks = n+2
      kf = ks+ns2-1
      fk = 0.
      do 101 k=ks,kf
         fk = fk+1.
         wsave(k) = 2.*sin(fk*dt)
  101 continue
      call rffti (np1,wsave(kf+1))
      return
      end
