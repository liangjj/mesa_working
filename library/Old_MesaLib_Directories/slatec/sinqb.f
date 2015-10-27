*deck sinqb
      subroutine sinqb (n, x, wsave)
c***begin prologue  sinqb
c***purpose  compute the unnormalized inverse of sinqf.
c***library   slatec (fftpack)
c***category  j1a3
c***type      single precision (sinqb-s)
c***keywords  fftpack, fourier transform
c***author  swarztrauber, p. n., (ncar)
c***description
c
c  subroutine sinqb computes the fast fourier transform of quarter
c  wave data.  that is, sinqb computes a sequence from its
c  representation in terms of a sine series with odd wave numbers.
c  the transform is defined below at output parameter x.
c
c  sinqf is the unnormalized inverse of sinqb since a call of sinqb
c  followed by a call of sinqf will multiply the input sequence x
c  by 4*n.
c
c  the array wsave which is used by subroutine sinqb must be
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
c          in the program that calls sinqb.  the wsave array must be
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
c               x(i)= the sum from k=1 to k=n of
c
c                 4*x(k)*sin((2*k-1)*i*pi/(2*n))
c
c               a call of sinqb followed by a call of
c               sinqf will multiply the sequence x by 4*n.
c               therefore sinqf is the unnormalized inverse
c               of sinqb.
c
c  wsave   contains initialization calculations which must not
c          be destroyed between calls of sinqb or sinqf.
c
c***references  p. n. swarztrauber, vectorizing the ffts, in parallel
c                 computations (g. rodrigue, ed.), academic press,
c                 1982, pp. 51-83.
c***routines called  cosqb
c***revision history  (yymmdd)
c   790601  date written
c   830401  modified to use slatec library source file format.
c   860115  modified by ron boisvert to adhere to fortran 77 by
c           changing dummy array size declarations (1) to (*).
c   861211  revision date from version 3.2
c   881128  modified by dick valent to meet prologue standards.
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  sinqb
      dimension x(*), wsave(*)
c***first executable statement  sinqb
      if (n .gt. 1) go to 101
      x(1) = 4.*x(1)
      return
  101 ns2 = n/2
      do 102 k=2,n,2
         x(k) = -x(k)
  102 continue
      call cosqb (n,x,wsave)
      do 103 k=1,ns2
         kc = n-k
         xhold = x(k)
         x(k) = x(kc+1)
         x(kc+1) = xhold
  103 continue
      return
      end
