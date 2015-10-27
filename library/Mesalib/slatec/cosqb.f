*deck cosqb
      subroutine cosqb (n, x, wsave)
c***begin prologue  cosqb
c***purpose  compute the unnormalized inverse cosine transform.
c***library   slatec (fftpack)
c***category  j1a3
c***type      single precision (cosqb-s)
c***keywords  fftpack, inverse cosine fourier transform
c***author  swarztrauber, p. n., (ncar)
c***description
c
c  subroutine cosqb computes the fast fourier transform of quarter
c  wave data. that is, cosqb computes a sequence from its
c  representation in terms of a cosine series with odd wave numbers.
c  the transform is defined below at output parameter x.
c
c  cosqb is the unnormalized inverse of cosqf since a call of cosqb
c  followed by a call of cosqf will multiply the input sequence x
c  by 4*n.
c
c  the array wsave which is used by subroutine cosqb must be
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
c          in the program that calls cosqb.  the wsave array must be
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
c               x(i)= the sum from k=1 to k=n of
c
c                  2*x(k)*cos((2*k-1)*(i-1)*pi/(2*n))
c
c               a call of cosqb followed by a call of
c               cosqf will multiply the sequence x by 4*n.
c               therefore cosqf is the unnormalized inverse
c               of cosqb.
c
c  wsave   contains initialization calculations which must not
c          be destroyed between calls of cosqb or cosqf.
c
c***references  p. n. swarztrauber, vectorizing the ffts, in parallel
c                 computations (g. rodrigue, ed.), academic press,
c                 1982, pp. 51-83.
c***routines called  cosqb1
c***revision history  (yymmdd)
c   790601  date written
c   830401  modified to use slatec library source file format.
c   860115  modified by ron boisvert to adhere to fortran 77 by
c           (a) changing dummy array size declarations (1) to (*),
c           (b) changing definition of variable tsqrt2 by using
c               fortran intrinsic function sqrt instead of a data
c               statement.
c   861211  revision date from version 3.2
c   881128  modified by dick valent to meet prologue standards.
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  cosqb
      dimension x(*), wsave(*)
c***first executable statement  cosqb
      tsqrt2 = 2.*sqrt(2.)
      if (n-2) 101,102,103
  101 x(1) = 4.*x(1)
      return
  102 x1 = 4.*(x(1)+x(2))
      x(2) = tsqrt2*(x(1)-x(2))
      x(1) = x1
      return
  103 call cosqb1 (n,x,wsave,wsave(n+1))
      return
      end
