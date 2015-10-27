*deck rfftb1
      subroutine rfftb1 (n, c, ch, wa, ifac)
c***begin prologue  rfftb1
c***purpose  compute the backward fast fourier transform of a real
c            coefficient array.
c***library   slatec (fftpack)
c***category  j1a1
c***type      single precision (rfftb1-s, cfftb1-c)
c***keywords  fftpack, fourier transform
c***author  swarztrauber, p. n., (ncar)
c***description
c
c   subroutine rfftb1 computes the real periodic sequence from its
c   fourier coefficients (fourier synthesis).  the transform is defined
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
c   c       for n even and for i = 1,...,n
c
c                c(i) = c(1)+(-1)**(i-1)*c(n)
c
c                     plus the sum from k=2 to k=n/2 of
c
c                      2.*c(2*k-2)*cos((k-1)*(i-1)*2*pi/n)
c
c                     -2.*c(2*k-1)*sin((k-1)*(i-1)*2*pi/n)
c
c           for n odd and for i = 1,...,n
c
c                c(i) = c(1) plus the sum from k=2 to k=(n+1)/2 of
c
c                     2.*c(2*k-2)*cos((k-1)*(i-1)*2*pi/n)
c
c                    -2.*c(2*k-1)*sin((k-1)*(i-1)*2*pi/n)
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
c***routines called  radb2, radb3, radb4, radb5, radbg
c***revision history  (yymmdd)
c   790601  date written
c   830401  modified to use slatec library source file format.
c   860115  modified by ron boisvert to adhere to fortran 77 by
c           changing dummy array size declarations (1) to (*).
c   881128  modified by dick valent to meet prologue standards.
c   891214  prologue converted to version 4.0 format.  (bab)
c   900131  routine changed from subsidiary to user-callable.  (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  rfftb1
      dimension ch(*), c(*), wa(*), ifac(*)
c***first executable statement  rfftb1
      nf = ifac(2)
      na = 0
      l1 = 1
      iw = 1
      do 116 k1=1,nf
         ip = ifac(k1+2)
         l2 = ip*l1
         ido = n/l2
         idl1 = ido*l1
         if (ip .ne. 4) go to 103
         ix2 = iw+ido
         ix3 = ix2+ido
         if (na .ne. 0) go to 101
         call radb4 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
         go to 102
  101    call radb4 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
  102    na = 1-na
         go to 115
  103    if (ip .ne. 2) go to 106
         if (na .ne. 0) go to 104
         call radb2 (ido,l1,c,ch,wa(iw))
         go to 105
  104    call radb2 (ido,l1,ch,c,wa(iw))
  105    na = 1-na
         go to 115
  106    if (ip .ne. 3) go to 109
         ix2 = iw+ido
         if (na .ne. 0) go to 107
         call radb3 (ido,l1,c,ch,wa(iw),wa(ix2))
         go to 108
  107    call radb3 (ido,l1,ch,c,wa(iw),wa(ix2))
  108    na = 1-na
         go to 115
  109    if (ip .ne. 5) go to 112
         ix2 = iw+ido
         ix3 = ix2+ido
         ix4 = ix3+ido
         if (na .ne. 0) go to 110
         call radb5 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 111
  110    call radb5 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3),wa(ix4))
  111    na = 1-na
         go to 115
  112    if (na .ne. 0) go to 113
         call radbg (ido,ip,l1,idl1,c,c,c,ch,ch,wa(iw))
         go to 114
  113    call radbg (ido,ip,l1,idl1,ch,ch,ch,c,c,wa(iw))
  114    if (ido .eq. 1) na = 1-na
  115    l1 = l2
         iw = iw+(ip-1)*ido
  116 continue
      if (na .eq. 0) return
      do 117 i=1,n
         c(i) = ch(i)
  117 continue
      return
      end
