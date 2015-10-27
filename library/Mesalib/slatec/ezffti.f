*deck ezffti
      subroutine ezffti (n, wsave)
c***begin prologue  ezffti
c***purpose  initialize a work array for ezfftf and ezfftb.
c***library   slatec (fftpack)
c***category  j1a1
c***type      single precision (ezffti-s)
c***keywords  fftpack, fourier transform
c***author  swarztrauber, p. n., (ncar)
c***description
c
c  subroutine ezffti initializes the work array wsave which is used in
c  both ezfftf and ezfftb.  the prime factorization of n together with
c  a tabulation of the trigonometric functions are computed and
c  stored in wsave.
c
c  input parameter
c
c  n       the length of the sequence to be transformed.
c
c  output parameter
c
c  wsave   a work array which must be dimensioned at least 3*n+15.
c          the same work array can be used for both ezfftf and ezfftb
c          as long as n remains unchanged.  different wsave arrays
c          are required for different values of n.
c
c***references  p. n. swarztrauber, vectorizing the ffts, in parallel
c                 computations (g. rodrigue, ed.), academic press,
c                 1982, pp. 51-83.
c***routines called  ezfft1
c***revision history  (yymmdd)
c   790601  date written
c   830401  modified to use slatec library source file format.
c   860115  modified by ron boisvert to adhere to fortran 77 by
c           changing dummy array size declarations (1) to (*).
c   861211  revision date from version 3.2
c   881128  modified by dick valent to meet prologue standards.
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  ezffti
      dimension wsave(*)
c***first executable statement  ezffti
      if (n .eq. 1) return
      call ezfft1 (n,wsave(2*n+1),wsave(3*n+1))
      return
      end
