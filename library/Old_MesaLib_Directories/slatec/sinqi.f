*deck sinqi
      subroutine sinqi (n, wsave)
c***begin prologue  sinqi
c***purpose  initialize a work array for sinqf and sinqb.
c***library   slatec (fftpack)
c***category  j1a3
c***type      single precision (sinqi-s)
c***keywords  fftpack, fourier transform
c***author  swarztrauber, p. n., (ncar)
c***description
c
c  subroutine sinqi initializes the array wsave which is used in
c  both sinqf and sinqb.  the prime factorization of n together with
c  a tabulation of the trigonometric functions are computed and
c  stored in wsave.
c
c  input parameter
c
c  n       the length of the sequence to be transformed.  the method
c          is most efficient when n is a product of small primes.
c
c  output parameter
c
c  wsave   a work array which must be dimensioned at least 3*n+15.
c          the same work array can be used for both sinqf and sinqb
c          as long as n remains unchanged.  different wsave arrays
c          are required for different values of n.  the contents of
c          wsave must not be changed between calls of sinqf or sinqb.
c
c***references  p. n. swarztrauber, vectorizing the ffts, in parallel
c                 computations (g. rodrigue, ed.), academic press,
c                 1982, pp. 51-83.
c***routines called  cosqi
c***revision history  (yymmdd)
c   790601  date written
c   830401  modified to use slatec library source file format.
c   860115  modified by ron boisvert to adhere to fortran 77 by
c           changing dummy array size declarations (1) to (*)
c   861211  revision date from version 3.2
c   881128  modified by dick valent to meet prologue standards.
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  sinqi
      dimension wsave(*)
c***first executable statement  sinqi
      call cosqi (n,wsave)
      return
      end
