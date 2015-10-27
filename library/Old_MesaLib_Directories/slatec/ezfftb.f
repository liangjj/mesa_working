*deck ezfftb
      subroutine ezfftb (n, r, azero, a, b, wsave)
c***begin prologue  ezfftb
c***purpose  a simplified real, periodic, backward fast fourier
c            transform.
c***library   slatec (fftpack)
c***category  j1a1
c***type      single precision (ezfftb-s)
c***keywords  fftpack, fourier transform
c***author  swarztrauber, p. n., (ncar)
c***description
c
c  subroutine ezfftb computes a real periodic sequence from its
c  fourier coefficients (fourier synthesis).  the transform is
c  defined below at output parameter r.  ezfftb is a simplified
c  but slower version of rfftb.
c
c  input parameters
c
c  n       the length of the output array r.  the method is most
c          efficient when n is the product of small primes.
c
c  azero   the constant fourier coefficient
c
c  a,b     arrays which contain the remaining fourier coefficients.
c          these arrays are not destroyed.
c
c          the length of these arrays depends on whether n is even or
c          odd.
c
c          if n is even, n/2    locations are required.
c          if n is odd, (n-1)/2 locations are required
c
c  wsave   a work array which must be dimensioned at least 3*n+15
c          in the program that calls ezfftb.  the wsave array must be
c          initialized by calling subroutine ezffti(n,wsave), and a
c          different wsave array must be used for each different
c          value of n.  this initialization does not have to be
c          repeated so long as n remains unchanged.  thus subsequent
c          transforms can be obtained faster than the first.
c          the same wsave array can be used by ezfftf and ezfftb.
c
c  output parameters
c
c  r       if n is even, define kmax=n/2
c          if n is odd,  define kmax=(n-1)/2
c
c          then for i=1,...,n
c
c               r(i)=azero plus the sum from k=1 to k=kmax of
c
c               a(k)*cos(k*(i-1)*2*pi/n)+b(k)*sin(k*(i-1)*2*pi/n)
c
c  ********************* complex notation **************************
c
c          for j=1,...,n
c
c          r(j) equals the sum from k=-kmax to k=kmax of
c
c               c(k)*exp(i*k*(j-1)*2*pi/n)
c
c          where
c
c               c(k) = .5*cmplx(a(k),-b(k))   for k=1,...,kmax
c
c               c(-k) = conjg(c(k))
c
c               c(0) = azero
c
c                    and i=sqrt(-1)
c
c  *************** amplitude - phase notation ***********************
c
c          for i=1,...,n
c
c          r(i) equals azero plus the sum from k=1 to k=kmax of
c
c               alpha(k)*cos(k*(i-1)*2*pi/n+beta(k))
c
c          where
c
c               alpha(k) = sqrt(a(k)*a(k)+b(k)*b(k))
c
c               cos(beta(k))=a(k)/alpha(k)
c
c               sin(beta(k))=-b(k)/alpha(k)
c
c***references  p. n. swarztrauber, vectorizing the ffts, in parallel
c                 computations (g. rodrigue, ed.), academic press,
c                 1982, pp. 51-83.
c***routines called  rfftb
c***revision history  (yymmdd)
c   790601  date written
c   830401  modified to use slatec library source file format.
c   860115  modified by ron boisvert to adhere to fortran 77 by
c           changing dummy array size declarations (1) to (*)
c   861211  revision date from version 3.2
c   881128  modified by dick valent to meet prologue standards.
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  ezfftb
      dimension r(*), a(*), b(*), wsave(*)
c***first executable statement  ezfftb
      if (n-2) 101,102,103
  101 r(1) = azero
      return
  102 r(1) = azero+a(1)
      r(2) = azero-a(1)
      return
  103 ns2 = (n-1)/2
      do 104 i=1,ns2
         r(2*i) = .5*a(i)
         r(2*i+1) = -.5*b(i)
  104 continue
      r(1) = azero
      if (mod(n,2) .eq. 0) r(n) = a(ns2+1)
      call rfftb (n,r,wsave(n+1))
      return
      end
