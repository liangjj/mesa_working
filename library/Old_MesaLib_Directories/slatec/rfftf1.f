*deck rfftf1
      subroutine rfftf1 (n, c, ch, wa, ifac)
c***begin prologue  rfftf1
c***purpose  compute the forward transform of a real, periodic sequence.
c***library   slatec (fftpack)
c***category  j1a1
c***type      single precision (rfftf1-s, cfftf1-c)
c***keywords  fftpack, fourier transform
c***author  swarztrauber, p. n., (ncar)
c***description
c
c   subroutine rfftf1 computes the fourier coefficients of a real
c   periodic sequence (fourier analysis).  the transform is defined
c   below at output parameter c.
c
c   the arrays wa and ifac which are used by subroutine rfftb1 must be
c   initialized by calling subroutine rffti1.
c
c   input arguments
c
c   n       the length of the array r to be transformed.  the method
c           is most efficient when n is a product of small primes.
c           n may change so long as different work arrays are provided.
c
c   c       a real array of length n which contains the sequence
c           to be transformed.
c
c   ch      a real work array of length at least n.
c
c   wa      a real work array which must be dimensioned at least n.
c
c   ifac    an integer work array which must be dimensioned at least 15.
c
c           the wa and ifac arrays must be initialized by calling
c           subroutine rffti1, and different wa and ifac arrays must be
c           used for each different value of n.  this initialization
c           does not have to be repeated so long as n remains unchanged.
c           thus subsequent transforms can be obtained faster than the
c           first.  the same wa and ifac arrays can be used by rfftf1
c           and rfftb1.
c
c   output argument
c
c   c       c(1) = the sum from i=1 to i=n of r(i)
c
c           if n is even set l = n/2; if n is odd set l = (n+1)/2
c
c             then for k = 2,...,l
c
c                c(2*k-2) = the sum from i = 1 to i = n of
c
c                     c(i)*cos((k-1)*(i-1)*2*pi/n)
c
c                c(2*k-1) = the sum from i = 1 to i = n of
c
c                    -c(i)*sin((k-1)*(i-1)*2*pi/n)
c
c           if n is even
c
c                c(n) = the sum from i = 1 to i = n of
c
c                     (-1)**(i-1)*c(i)
c
c   notes:  this transform is unnormalized since a call of rfftf1
c           followed by a call of rfftb1 will multiply the input
c           sequence by n.
c
c           wa and ifac contain initialization calculations which must
c           not be destroyed between calls of subroutine rfftf1 or
c           rfftb1.
c
c***references  p. n. swarztrauber, vectorizing the ffts, in parallel
c                 computations (g. rodrigue, ed.), academic press,
c                 1982, pp. 51-83.
c***routines called  radf2, radf3, radf4, radf5, radfg
c***revision history  (yymmdd)
c   790601  date written
c   830401  modified to use slatec library source file format.
c   860115  modified by ron boisvert to adhere to fortran 77 by
c           changing dummy array size declarations (1) to (*).
c   881128  modified by dick valent to meet prologue standards.
c   891214  prologue converted to version 4.0 format.  (bab)
c   900131  routine changed from subsidiary to user-callable.  (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  rfftf1
      dimension ch(*), c(*), wa(*), ifac(*)
c***first executable statement  rfftf1
      nf = ifac(2)
      na = 1
      l2 = n
      iw = n
      do 111 k1=1,nf
         kh = nf-k1
         ip = ifac(kh+3)
         l1 = l2/ip
         ido = n/l2
         idl1 = ido*l1
         iw = iw-(ip-1)*ido
         na = 1-na
         if (ip .ne. 4) go to 102
         ix2 = iw+ido
         ix3 = ix2+ido
         if (na .ne. 0) go to 101
         call radf4 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
         go to 110
  101    call radf4 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
         go to 110
  102    if (ip .ne. 2) go to 104
         if (na .ne. 0) go to 103
         call radf2 (ido,l1,c,ch,wa(iw))
         go to 110
  103    call radf2 (ido,l1,ch,c,wa(iw))
         go to 110
  104    if (ip .ne. 3) go to 106
         ix2 = iw+ido
         if (na .ne. 0) go to 105
         call radf3 (ido,l1,c,ch,wa(iw),wa(ix2))
         go to 110
  105    call radf3 (ido,l1,ch,c,wa(iw),wa(ix2))
         go to 110
  106    if (ip .ne. 5) go to 108
         ix2 = iw+ido
         ix3 = ix2+ido
         ix4 = ix3+ido
         if (na .ne. 0) go to 107
         call radf5 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 110
  107    call radf5 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 110
  108    if (ido .eq. 1) na = 1-na
         if (na .ne. 0) go to 109
         call radfg (ido,ip,l1,idl1,c,c,c,ch,ch,wa(iw))
         na = 1
         go to 110
  109    call radfg (ido,ip,l1,idl1,ch,ch,ch,c,c,wa(iw))
         na = 0
  110    l2 = l1
  111 continue
      if (na .eq. 1) return
      do 112 i=1,n
         c(i) = ch(i)
  112 continue
      return
      end
