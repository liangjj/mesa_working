*deck cost
      subroutine cost (n, x, wsave)
c***begin prologue  cost
c***purpose  compute the cosine transform of a real, even sequence.
c***library   slatec (fftpack)
c***category  j1a3
c***type      single precision (cost-s)
c***keywords  cosine fourier transform, fftpack
c***author  swarztrauber, p. n., (ncar)
c***description
c
c  subroutine cost computes the discrete fourier cosine transform
c  of an even sequence x(i).  the transform is defined below at output
c  parameter x.
c
c  cost is the unnormalized inverse of itself since a call of cost
c  followed by another call of cost will multiply the input sequence
c  x by 2*(n-1).  the transform is defined below at output parameter x.
c
c  the array wsave which is used by subroutine cost must be
c  initialized by calling subroutine costi(n,wsave).
c
c  input parameters
c
c  n       the length of the sequence x.  n must be greater than 1.
c          the method is most efficient when n-1 is a product of
c          small primes.
c
c  x       an array which contains the sequence to be transformed
c
c  wsave   a work array which must be dimensioned at least 3*n+15
c          in the program that calls cost.  the wsave array must be
c          initialized by calling subroutine costi(n,wsave), and a
c          different wsave array must be used for each different
c          value of n.  this initialization does not have to be
c          repeated so long as n remains unchanged.  thus subsequent
c          transforms can be obtained faster than the first.
c
c  output parameters
c
c  x       for i=1,...,n
c
c             x(i) = x(1)+(-1)**(i-1)*x(n)
c
c               + the sum from k=2 to k=n-1
c
c                 2*x(k)*cos((k-1)*(i-1)*pi/(n-1))
c
c               a call of cost followed by another call of
c               cost will multiply the sequence x by 2*(n-1).
c               hence cost is the unnormalized inverse
c               of itself.
c
c  wsave   contains initialization calculations which must not be
c          destroyed between calls of cost.
c
c***references  p. n. swarztrauber, vectorizing the ffts, in parallel
c                 computations (g. rodrigue, ed.), academic press,
c                 1982, pp. 51-83.
c***routines called  rfftf
c***revision history  (yymmdd)
c   790601  date written
c   830401  modified to use slatec library source file format.
c   860115  modified by ron boisvert to adhere to fortran 77 by
c           changing dummy array size declarations (1) to (*)
c   861211  revision date from version 3.2
c   881128  modified by dick valent to meet prologue standards.
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  cost
      dimension x(*), wsave(*)
c***first executable statement  cost
      nm1 = n-1
      np1 = n+1
      ns2 = n/2
      if (n-2) 106,101,102
  101 x1h = x(1)+x(2)
      x(2) = x(1)-x(2)
      x(1) = x1h
      return
  102 if (n .gt. 3) go to 103
      x1p3 = x(1)+x(3)
      tx2 = x(2)+x(2)
      x(2) = x(1)-x(3)
      x(1) = x1p3+tx2
      x(3) = x1p3-tx2
      return
  103 c1 = x(1)-x(n)
      x(1) = x(1)+x(n)
      do 104 k=2,ns2
         kc = np1-k
         t1 = x(k)+x(kc)
         t2 = x(k)-x(kc)
         c1 = c1+wsave(kc)*t2
         t2 = wsave(k)*t2
         x(k) = t1-t2
         x(kc) = t1+t2
  104 continue
      modn = mod(n,2)
      if (modn .ne. 0) x(ns2+1) = x(ns2+1)+x(ns2+1)
      call rfftf (nm1,x,wsave(n+1))
      xim2 = x(2)
      x(2) = c1
      do 105 i=4,n,2
         xi = x(i)
         x(i) = x(i-2)-x(i-1)
         x(i-1) = xim2
         xim2 = xi
  105 continue
      if (modn .ne. 0) x(n) = xim2
  106 return
      end
