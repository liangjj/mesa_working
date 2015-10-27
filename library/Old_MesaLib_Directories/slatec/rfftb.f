*deck rfftb
      subroutine rfftb (n, r, wsave)
c***begin prologue  rfftb
c***subsidiary
c***purpose  compute the backward fast fourier transform of a real
c            coefficient array.
c***library   slatec (fftpack)
c***category  j1a1
c***type      single precision (rfftb-s, cfftb-c)
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
c   *   requested to use rfftb1.                                       *
c   *                                                                  *
c   ********************************************************************
c
c   subroutine rfftb computes the real periodic sequence from its
c   fourier coefficients (fourier synthesis).  the transform is defined
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
c           in the program that calls rfftb.  the wsave array must be
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
c   r       for n even and for i = 1,...,n
c
c                r(i) = r(1)+(-1)**(i-1)*r(n)
c
c                     plus the sum from k=2 to k=n/2 of
c
c                      2.*r(2*k-2)*cos((k-1)*(i-1)*2*pi/n)
c
c                     -2.*r(2*k-1)*sin((k-1)*(i-1)*2*pi/n)
c
c           for n odd and for i = 1,...,n
c
c                r(i) = r(1) plus the sum from k=2 to k=(n+1)/2 of
c
c                     2.*r(2*k-2)*cos((k-1)*(i-1)*2*pi/n)
c
c                    -2.*r(2*k-1)*sin((k-1)*(i-1)*2*pi/n)
c
c   note:  this transform is unnormalized since a call of rfftf
c          followed by a call of rfftb will multiply the input
c          sequence by n.
c
c   wsave  contains results which must not be destroyed between
c          calls of rfftb or rfftf.
c
c***references  p. n. swarztrauber, vectorizing the ffts, in parallel
c                 computations (g. rodrigue, ed.), academic press,
c                 1982, pp. 51-83.
c***routines called  rfftb1
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
c***end prologue  rfftb
      dimension r(*), wsave(*)
c***first executable statement  rfftb
      if (n .eq. 1) return
      call rfftb1 (n,r,wsave,wsave(n+1),wsave(2*n+1))
      return
      end
