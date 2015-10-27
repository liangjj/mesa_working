*deck ezfftf
      subroutine ezfftf (n, r, azero, a, b, wsave)
c***begin prologue  ezfftf
c***purpose  compute a simplified real, periodic, fast fourier forward
c            transform.
c***library   slatec (fftpack)
c***category  j1a1
c***type      single precision (ezfftf-s)
c***keywords  fftpack, fourier transform
c***author  swarztrauber, p. n., (ncar)
c***description
c
c  subroutine ezfftf computes the fourier coefficients of a real
c  periodic sequence (fourier analysis).  the transform is defined
c  below at output parameters azero, a and b.  ezfftf is a simplified
c  but slower version of rfftf.
c
c  input parameters
c
c  n       the length of the array r to be transformed.  the method
c          is most efficient when n is the product of small primes.
c
c  r       a real array of length n which contains the sequence
c          to be transformed.  r is not destroyed.
c
c
c  wsave   a work array which must be dimensioned at least 3*n+15
c          in the program that calls ezfftf.  the wsave array must be
c          initialized by calling subroutine ezffti(n,wsave), and a
c          different wsave array must be used for each different
c          value of n.  this initialization does not have to be
c          repeated so long as n remains unchanged.  thus subsequent
c          transforms can be obtained faster than the first.
c          the same wsave array can be used by ezfftf and ezfftb.
c
c  output parameters
c
c  azero   the sum from i=1 to i=n of r(i)/n
c
c  a,b     for n even b(n/2)=0. and a(n/2) is the sum from i=1 to
c          i=n of (-1)**(i-1)*r(i)/n
c
c          for n even define kmax=n/2-1
c          for n odd  define kmax=(n-1)/2
c
c          then for  k=1,...,kmax
c
c               a(k) equals the sum from i=1 to i=n of
c
c                    2./n*r(i)*cos(k*(i-1)*2*pi/n)
c
c               b(k) equals the sum from i=1 to i=n of
c
c                    2./n*r(i)*sin(k*(i-1)*2*pi/n)
c
c***references  p. n. swarztrauber, vectorizing the ffts, in parallel
c                 computations (g. rodrigue, ed.), academic press,
c                 1982, pp. 51-83.
c***routines called  rfftf
c***revision history  (yymmdd)
c   790601  date written
c   830401  modified to use slatec library source file format.
c   860115  modified by ron boisvert to adhere to fortran 77 by
c           (a) changing dummy array size declarations (1) to (*),
c           (b) changing references to intrinsic function float
c               to real.
c   881128  modified by dick valent to meet prologue standards.
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  ezfftf
      dimension r(*), a(*), b(*), wsave(*)
c***first executable statement  ezfftf
      if (n-2) 101,102,103
  101 azero = r(1)
      return
  102 azero = .5*(r(1)+r(2))
      a(1) = .5*(r(1)-r(2))
      return
  103 do 104 i=1,n
         wsave(i) = r(i)
  104 continue
      call rfftf (n,wsave,wsave(n+1))
      cf = 2./n
      cfm = -cf
      azero = .5*cf*wsave(1)
      ns2 = (n+1)/2
      ns2m = ns2-1
      do 105 i=1,ns2m
         a(i) = cf*wsave(2*i)
         b(i) = cfm*wsave(2*i+1)
  105 continue
      if (mod(n,2) .eq. 0) a(ns2) = .5*cf*wsave(n)
      return
      end
