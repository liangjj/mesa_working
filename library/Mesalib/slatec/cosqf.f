*deck cosqf
      subroutine cosqf (n, x, wsave)
c***begin prologue  cosqf
c***purpose  compute the forward cosine transform with odd wave numbers.
c***library   slatec (fftpack)
c***category  j1a3
c***type      single precision (cosqf-s)
c***keywords  cosine fourier transform, fftpack
c***author  swarztrauber, p. n., (ncar)
c***description
c
c  subroutine cosqf computes the fast fourier transform of quarter
c  wave data. that is, cosqf computes the coefficients in a cosine
c  series representation with only odd wave numbers.  the transform
c  is defined below at output parameter x
c
c  cosqf is the unnormalized inverse of cosqb since a call of cosqf
c  followed by a call of cosqb will multiply the input sequence x
c  by 4*n.
c
c  the array wsave which is used by subroutine cosqf must be
c  initialized by calling subroutine cosqi(n,wsave).
c
c
c  input parameters
c
c  n       the length of the array x to be transformed.  the method
c          is most efficient when n is a product of small primes.
c
c  x       an array which contains the sequence to be transformed
c
c  wsave   a work array which must be dimensioned at least 3*n+15
c          in the program that calls cosqf.  the wsave array must be
c          initialized by calling subroutine cosqi(n,wsave), and a
c          different wsave array must be used for each different
c          value of n.  this initialization does not have to be
c          repeated so long as n remains unchanged.  thus subsequent
c          transforms can be obtained faster than the first.
c
c  output parameters
c
c  x       for i=1,...,n
c
c               x(i) = x(1) plus the sum from k=2 to k=n of
c
c                  2*x(k)*cos((2*i-1)*(k-1)*pi/(2*n))
c
c               a call of cosqf followed by a call of
c               cosqb will multiply the sequence x by 4*n.
c               therefore cosqb is the unnormalized inverse
c               of cosqf.
c
c  wsave   contains initialization calculations which must not
c          be destroyed between calls of cosqf or cosqb.
c
c***references  p. n. swarztrauber, vectorizing the ffts, in parallel
c                 computations (g. rodrigue, ed.), academic press,
c                 1982, pp. 51-83.
c***routines called  cosqf1
c***revision history  (yymmdd)
c   790601  date written
c   830401  modified to use slatec library source file format.
c   860115  modified by ron boisvert to adhere to fortran 77 by
c           (a) changing dummy array size declarations (1) to (*),
c           (b) changing definition of variable sqrt2 by using
c               fortran intrinsic function sqrt instead of a data
c               statement.
c   861211  revision date from version 3.2
c   881128  modified by dick valent to meet prologue standards.
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  cosqf
      dimension x(*), wsave(*)
c***first executable statement  cosqf
      sqrt2 = sqrt(2.)
      if (n-2) 102,101,103
  101 tsqx = sqrt2*x(2)
      x(2) = x(1)-tsqx
      x(1) = x(1)+tsqx
  102 return
  103 call cosqf1 (n,x,wsave,wsave(n+1))
      return
      end
