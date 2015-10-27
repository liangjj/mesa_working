*deck rfftf
      subroutine rfftf (n, r, wsave)
c***begin prologue  rfftf
c***subsidiary
c***purpose  compute the forward transform of a real, periodic sequence.
c***library   slatec (fftpack)
c***category  j1a1
c***type      single precision (rfftf-s, cfftf-c)
c***keywords  fftpack, fourier transform
c***author  swarztrauber, p. n., (ncar)
c***description
c
c   ********************************************************************
c   *   notice   notice   notice   notice   notice   notice   notice   *
c   ********************************************************************
c   *                                                                  *
c   *   this routine uses non-standard fortran 77 constructs and will  *
c   *   be removed from the library at a future date.  you are         *
c   *   requested to use rfftf1.                                       *
c   *                                                                  *
c   ********************************************************************
c
c   subroutine rfftf computes the fourier coefficients of a real
c   periodic sequence (fourier analysis).  the transform is defined
c   below at output parameter r.
c
c   input arguments
c
c   n       the length of the array r to be transformed.  the method
c           is most efficient when n is a product of small primes.
c           n may change so long as different work arrays are provided.
c
c   r       a real array of length n which contains the sequence
c           to be transformed.
c
c   wsave   a work array which must be dimensioned at least 2*n+15
c           in the program that calls rfftf.  the wsave array must be
c           initialized by calling subroutine rffti, and a different
c           wsave array must be used for each different value of n.
c           this initialization does not have to be repeated so long as
c           remains unchanged.  thus subsequent transforms can be
c           obtained faster than the first.  moreover, the same wsave
c           array can be used by rfftf and rfftb as long as n remains
c           unchanged.
c
c   output argument
c
c   r       r(1) = the sum from i=1 to i=n of r(i)
c
c           if n is even set l = n/2; if n is odd set l = (n+1)/2
c
c             then for k = 2,...,l
c
c                r(2*k-2) = the sum from i = 1 to i = n of
c
c                     r(i)*cos((k-1)*(i-1)*2*pi/n)
c
c                r(2*k-1) = the sum from i = 1 to i = n of
c
c                    -r(i)*sin((k-1)*(i-1)*2*pi/n)
c
c           if n is even
c
c                r(n) = the sum from i = 1 to i = n of
c
c                     (-1)**(i-1)*r(i)
c
c   note:  this transform is unnormalized since a call of rfftf
c          followed by a call of rfftb will multiply the input
c          sequence by n.
c
c   wsave  contains results which must not be destroyed between
c          calls of rfftf or rfftb.
c
c***references  p. n. swarztrauber, vectorizing the ffts, in parallel
c                 computations (g. rodrigue, ed.), academic press,
c                 1982, pp. 51-83.
c***routines called  rfftf1
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
c***end prologue  rfftf
      dimension r(*), wsave(*)
c***first executable statement  rfftf
      if (n .eq. 1) return
      call rfftf1 (n,r,wsave,wsave(n+1),wsave(2*n+1))
      return
      end
