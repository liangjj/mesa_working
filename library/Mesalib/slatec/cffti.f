*deck cffti
      subroutine cffti (n, wsave)
c***begin prologue  cffti
c***subsidiary
c***purpose  initialize a work array for cfftf and cfftb.
c***library   slatec (fftpack)
c***category  j1a2
c***type      complex (rffti-s, cffti-c)
c***keywords  fftpack, fourier transform
c***author  swarztrauber, p. n., (ncar)
c***description
c
c  ********************************************************************
c  *   notice   notice   notice   notice   notice   notice   notice   *
c  ********************************************************************
c  *                                                                  *
c  *   this routine uses non-standard fortran 77 constructs and will  *
c  *   be removed from the library at a future date.  you are         *
c  *   requested to use cffti1.                                       *
c  *                                                                  *
c  ********************************************************************
c
c  subroutine cffti initializes the array wsave which is used in
c  both cfftf and cfftb.  the prime factorization of n together with
c  a tabulation of the trigonometric functions are computed and
c  stored in wsave.
c
c  input parameter
c
c  n       the length of the sequence to be transformed
c
c  output parameter
c
c  wsave   a work array which must be dimensioned at least 4*n+15.
c          the same work array can be used for both cfftf and cfftb
c          as long as n remains unchanged.  different wsave arrays
c          are required for different values of n.  the contents of
c          wsave must not be changed between calls of cfftf or cfftb.
c
c***references  p. n. swarztrauber, vectorizing the ffts, in parallel
c                 computations (g. rodrigue, ed.), academic press,
c                 1982, pp. 51-83.
c***routines called  cffti1
c***revision history  (yymmdd)
c   790601  date written
c   830401  modified to use slatec library source file format.
c   860115  modified by ron boisvert to adhere to fortran 77 by
c           changing dummy array size declarations (1) to (*).
c   861211  revision date from version 3.2
c   881128  modified by dick valent to meet prologue standards.
c   891214  prologue converted to version 4.0 format.  (bab)
c   900131  routine changed from user-callable to subsidiary
c           because of non-standard fortran 77 arguments in the
c           call to cfftb1.  (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  cffti
      dimension wsave(*)
c***first executable statement  cffti
      if (n .eq. 1) return
      iw1 = n+n+1
      iw2 = iw1+n+n
      call cffti1 (n,wsave(iw1),wsave(iw2))
      return
      end
