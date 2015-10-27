*deck sinqf
      subroutine sinqf (n, x, wsave)
c***begin prologue  sinqf
c***purpose  compute the forward sine transform with odd wave numbers.
c***library   slatec (fftpack)
c***category  j1a3
c***type      single precision (sinqf-s)
c***keywords  fftpack, fourier transform
c***author  swarztrauber, p. n., (ncar)
c***description
c
c  subroutine sinqf computes the fast fourier transform of quarter
c  wave data.  that is, sinqf computes the coefficients in a sine
c  series representation with only odd wave numbers.  the transform
c  is defined below at output parameter x.
c
c  sinqb is the unnormalized inverse of sinqf since a call of sinqf
c  followed by a call of sinqb will multiply the input sequence x
c  by 4*n.
c
c  the array wsave which is used by subroutine sinqf must be
c  initialized by calling subroutine sinqi(n,wsave).
c
c  input parameters
c
c  n       the length of the array x to be transformed.  the method
c          is most efficient when n is a product of small primes.
c
c  x       an array which contains the sequence to be transformed
c
c  wsave   a work array which must be dimensioned at least 3*n+15
c          in the program that calls sinqf.  the wsave array must be
c          initialized by calling subroutine sinqi(n,wsave), and a
c          different wsave array must be used for each different
c          value of n.  this initialization does not have to be
c          repeated so long as n remains unchanged.  thus subsequent
c          transforms can be obtained faster than the first.
c
c  output parameters
c
c  x       for i=1,...,n
c
c               x(i) = (-1)**(i-1)*x(n)
c
c                  + the sum from k=1 to k=n-1 of
c
c                  2*x(k)*sin((2*i-1)*k*pi/(2*n))
c
c               a call of sinqf followed by a call of
c               sinqb will multiply the sequence x by 4*n.
c               therefore sinqb is the unnormalized inverse
c               of sinqf.
c
c  wsave   contains initialization calculations which must not
c          be destroyed between calls of sinqf or sinqb.
c
c***references  p. n. swarztrauber, vectorizing the ffts, in parallel
c                 computations (g. rodrigue, ed.), academic press,
c                 1982, pp. 51-83.
c***routines called  cosqf
c***revision history  (yymmdd)
c   790601  date written
c   830401  modified to use slatec library source file format.
c   860115  modified by ron boisvert to adhere to fortran 77 by
c           changing dummy array size declarations (1) to (*)
c   861211  revision date from version 3.2
c   881128  modified by dick valent to meet prologue standards.
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  sinqf
      dimension x(*), wsave(*)
c***first executable statement  sinqf
      if (n .eq. 1) return
      ns2 = n/2
      do 101 k=1,ns2
         kc = n-k
         xhold = x(k)
         x(k) = x(kc+1)
         x(kc+1) = xhold
  101 continue
      call cosqf (n,x,wsave)
      do 102 k=2,n,2
         x(k) = -x(k)
  102 continue
      return
      end
