*deck cfftb1
      subroutine cfftb1 (n, c, ch, wa, ifac)
c***begin prologue  cfftb1
c***purpose  compute the unnormalized inverse of cfftf1.
c***library   slatec (fftpack)
c***category  j1a2
c***type      complex (rfftb1-s, cfftb1-c)
c***keywords  fftpack, fourier transform
c***author  swarztrauber, p. n., (ncar)
c***description
c
c  subroutine cfftb1 computes the backward complex discrete fourier
c  transform (the fourier synthesis).  equivalently, cfftb1 computes
c  a complex periodic sequence from its fourier coefficients.
c  the transform is defined below at output parameter c.
c
c  a call of cfftf1 followed by a call of cfftb1 will multiply the
c  sequence by n.
c
c  the arrays wa and ifac which are used by subroutine cfftb1 must be
c  initialized by calling subroutine cffti1 (n, wa, ifac).
c
c  input parameters
c
c  n       the length of the complex sequence c.  the method is
c          more efficient when n is the product of small primes.
c
c  c       a complex array of length n which contains the sequence
c
c  ch      a real work array of length at least 2*n
c
c  wa      a real work array which must be dimensioned at least 2*n.
c
c  ifac    an integer work array which must be dimensioned at least 15.
c
c          the wa and ifac arrays must be initialized by calling
c          subroutine cffti1 (n, wa, ifac), and different wa and ifac
c          arrays must be used for each different value of n.  this
c          initialization does not have to be repeated so long as n
c          remains unchanged.  thus subsequent transforms can be
c          obtained faster than the first.  the same wa and ifac arrays
c          can be used by cfftf1 and cfftb1.
c
c  output parameters
c
c  c       for j=1,...,n
c
c              c(j)=the sum from k=1,...,n of
c
c                 c(k)*exp(i*(j-1)*(k-1)*2*pi/n)
c
c                         where i=sqrt(-1)
c
c  note:   wa and ifac contain initialization calculations which must
c          not be destroyed between calls of subroutine cfftf1 or cfftb1
c
c***references  p. n. swarztrauber, vectorizing the ffts, in parallel
c                 computations (g. rodrigue, ed.), academic press,
c                 1982, pp. 51-83.
c***routines called  passb, passb2, passb3, passb4, passb5
c***revision history  (yymmdd)
c   790601  date written
c   830401  modified to use slatec library source file format.
c   860115  modified by ron boisvert to adhere to fortran 77 by
c           changing dummy array size declarations (1) to (*).
c   881128  modified by dick valent to meet prologue standards.
c   891214  prologue converted to version 4.0 format.  (bab)
c   900131  routine changed from subsidiary to user-callable.  (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  cfftb1
      dimension ch(*), c(*), wa(*), ifac(*)
c***first executable statement  cfftb1
      nf = ifac(2)
      na = 0
      l1 = 1
      iw = 1
      do 116 k1=1,nf
         ip = ifac(k1+2)
         l2 = ip*l1
         ido = n/l2
         idot = ido+ido
         idl1 = idot*l1
         if (ip .ne. 4) go to 103
         ix2 = iw+idot
         ix3 = ix2+idot
         if (na .ne. 0) go to 101
         call passb4 (idot,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
         go to 102
  101    call passb4 (idot,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
  102    na = 1-na
         go to 115
  103    if (ip .ne. 2) go to 106
         if (na .ne. 0) go to 104
         call passb2 (idot,l1,c,ch,wa(iw))
         go to 105
  104    call passb2 (idot,l1,ch,c,wa(iw))
  105    na = 1-na
         go to 115
  106    if (ip .ne. 3) go to 109
         ix2 = iw+idot
         if (na .ne. 0) go to 107
         call passb3 (idot,l1,c,ch,wa(iw),wa(ix2))
         go to 108
  107    call passb3 (idot,l1,ch,c,wa(iw),wa(ix2))
  108    na = 1-na
         go to 115
  109    if (ip .ne. 5) go to 112
         ix2 = iw+idot
         ix3 = ix2+idot
         ix4 = ix3+idot
         if (na .ne. 0) go to 110
         call passb5 (idot,l1,c,ch,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 111
  110    call passb5 (idot,l1,ch,c,wa(iw),wa(ix2),wa(ix3),wa(ix4))
  111    na = 1-na
         go to 115
  112    if (na .ne. 0) go to 113
         call passb (nac,idot,ip,l1,idl1,c,c,c,ch,ch,wa(iw))
         go to 114
  113    call passb (nac,idot,ip,l1,idl1,ch,ch,ch,c,c,wa(iw))
  114    if (nac .ne. 0) na = 1-na
  115    l1 = l2
         iw = iw+(ip-1)*idot
  116 continue
      if (na .eq. 0) return
      n2 = n+n
      do 117 i=1,n2
         c(i) = ch(i)
  117 continue
      return
      end
