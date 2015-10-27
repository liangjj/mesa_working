*deck cffti1
      subroutine cffti1 (n, wa, ifac)
c***begin prologue  cffti1
c***purpose  initialize a real and an integer work array for cfftf1 and
c            cfftb1.
c***library   slatec (fftpack)
c***category  j1a2
c***type      complex (rffti1-s, cffti1-c)
c***keywords  fftpack, fourier transform
c***author  swarztrauber, p. n., (ncar)
c***description
c
c  subroutine cffti1 initializes the work arrays wa and ifac which are
c  used in both cfftf1 and cfftb1.  the prime factorization of n and a
c  tabulation of the trigonometric functions are computed and stored in
c  ifac and wa, respectively.
c
c  input parameter
c
c  n       the length of the sequence to be transformed
c
c  output parameters
c
c  wa      a real work array which must be dimensioned at least 2*n.
c
c  ifac    an integer work array which must be dimensioned at least 15.
c
c          the same work arrays can be used for both cfftf1 and cfftb1
c          as long as n remains unchanged.  different wa and ifac arrays
c          are required for different values of n.  the contents of
c          wa and ifac must not be changed between calls of cfftf1 or
c          cfftb1.
c
c***references  p. n. swarztrauber, vectorizing the ffts, in parallel
c                 computations (g. rodrigue, ed.), academic press,
c                 1982, pp. 51-83.
c***routines called  (none)
c***revision history  (yymmdd)
c   790601  date written
c   830401  modified to use slatec library source file format.
c   860115  modified by ron boisvert to adhere to fortran 77 by
c           (a) changing dummy array size declarations (1) to (*),
c           (b) changing references to intrinsic function float
c               to real, and
c           (c) changing definition of variable tpi by using
c               fortran intrinsic function atan instead of a data
c               statement.
c   881128  modified by dick valent to meet prologue standards.
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900131  routine changed from subsidiary to user-callable.  (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  cffti1
      dimension wa(*), ifac(*), ntryh(4)
      save ntryh
      data ntryh(1),ntryh(2),ntryh(3),ntryh(4)/3,4,2,5/
c***first executable statement  cffti1
      nl = n
      nf = 0
      j = 0
  101 j = j+1
      if (j-4) 102,102,103
  102 ntry = ntryh(j)
      go to 104
  103 ntry = ntry+2
  104 nq = nl/ntry
      nr = nl-ntry*nq
      if (nr) 101,105,101
  105 nf = nf+1
      ifac(nf+2) = ntry
      nl = nq
      if (ntry .ne. 2) go to 107
      if (nf .eq. 1) go to 107
      do 106 i=2,nf
         ib = nf-i+2
         ifac(ib+2) = ifac(ib+1)
  106 continue
      ifac(3) = 2
  107 if (nl .ne. 1) go to 104
      ifac(1) = n
      ifac(2) = nf
      tpi = 8.*atan(1.)
      argh = tpi/n
      i = 2
      l1 = 1
      do 110 k1=1,nf
         ip = ifac(k1+2)
         ld = 0
         l2 = l1*ip
         ido = n/l2
         idot = ido+ido+2
         ipm = ip-1
         do 109 j=1,ipm
            i1 = i
            wa(i-1) = 1.
            wa(i) = 0.
            ld = ld+l1
            fi = 0.
            argld = ld*argh
            do 108 ii=4,idot,2
               i = i+2
               fi = fi+1.
               arg = fi*argld
               wa(i-1) = cos(arg)
               wa(i) = sin(arg)
  108       continue
            if (ip .le. 5) go to 109
            wa(i1-1) = wa(i-1)
            wa(i1) = wa(i)
  109    continue
         l1 = l2
  110 continue
      return
      end
