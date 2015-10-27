*deck sint
      subroutine sint (n, x, wsave)
c***begin prologue  sint
c***purpose  compute the sine transform of a real, odd sequence.
c***library   slatec (fftpack)
c***category  j1a3
c***type      single precision (sint-s)
c***keywords  fftpack, fourier transform
c***author  swarztrauber, p. n., (ncar)
c***description
c
c  subroutine sint computes the discrete fourier sine transform
c  of an odd sequence x(i).  the transform is defined below at
c  output parameter x.
c
c  sint is the unnormalized inverse of itself since a call of sint
c  followed by another call of sint will multiply the input sequence
c  x by 2*(n+1).
c
c  the array wsave which is used by subroutine sint must be
c  initialized by calling subroutine sinti(n,wsave).
c
c  input parameters
c
c  n       the length of the sequence to be transformed.  the method
c          is most efficient when n+1 is the product of small primes.
c
c  x       an array which contains the sequence to be transformed
c
c
c  wsave   a work array with dimension at least int(3.5*n+16)
c          in the program that calls sint.  the wsave array must be
c          initialized by calling subroutine sinti(n,wsave), and a
c          different wsave array must be used for each different
c          value of n.  this initialization does not have to be
c          repeated so long as n remains unchanged.  thus subsequent
c          transforms can be obtained faster than the first.
c
c  output parameters
c
c  x       for i=1,...,n
c
c               x(i)= the sum from k=1 to k=n
c
c                    2*x(k)*sin(k*i*pi/(n+1))
c
c               a call of sint followed by another call of
c               sint will multiply the sequence x by 2*(n+1).
c               hence sint is the unnormalized inverse
c               of itself.
c
c  wsave   contains initialization calculations which must not be
c          destroyed between calls of sint.
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
c           (b) changing definition of variable sqrt3 by using
c               fortran intrinsic function sqrt instead of a data
c               statement.
c   881128  modified by dick valent to meet prologue standards.
c   891009  removed unreferenced statement label.  (wrb)
c   891009  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  sint
      dimension x(*), wsave(*)
c***first executable statement  sint
      sqrt3 = sqrt(3.)
      if (n-2) 101,102,103
  101 x(1) = x(1)+x(1)
      return
  102 xh = sqrt3*(x(1)+x(2))
      x(2) = sqrt3*(x(1)-x(2))
      x(1) = xh
      return
  103 np1 = n+1
      ns2 = n/2
      wsave(1) = 0.
      kw = np1
      do 104 k=1,ns2
         kw = kw+1
         kc = np1-k
         t1 = x(k)-x(kc)
         t2 = wsave(kw)*(x(k)+x(kc))
         wsave(k+1) = t1+t2
         wsave(kc+1) = t2-t1
  104 continue
      modn = mod(n,2)
      if (modn .ne. 0) wsave(ns2+2) = 4.*x(ns2+1)
      nf = np1+ns2+1
      call rfftf (np1,wsave,wsave(nf))
      x(1) = .5*wsave(1)
      do 105 i=3,n,2
         x(i-1) = -wsave(i)
         x(i) = x(i-2)+wsave(i-1)
  105 continue
      if (modn .ne. 0) return
      x(n) = -wsave(n+1)
      return
      end
