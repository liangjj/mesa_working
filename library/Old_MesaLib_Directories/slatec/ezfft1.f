*deck ezfft1
      subroutine ezfft1 (n, wa, ifac)
c***begin prologue  ezfft1
c***subsidiary
c***purpose  ezffti calls ezfft1 with appropriate work array
c            partitioning.
c***library   slatec (fftpack)
c***type      single precision (ezfft1-s)
c***author  swarztrauber, p. n., (ncar)
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
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c***end prologue  ezfft1
      dimension wa(*), ifac(*), ntryh(4)
      save ntryh
      data ntryh(1),ntryh(2),ntryh(3),ntryh(4)/4,2,3,5/
c***first executable statement  ezfft1
      tpi = 8.*atan(1.)
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
      argh = tpi/n
      is = 0
      nfm1 = nf-1
      l1 = 1
      if (nfm1 .eq. 0) return
      do 111 k1=1,nfm1
         ip = ifac(k1+2)
         l2 = l1*ip
         ido = n/l2
         ipm = ip-1
         arg1 = l1*argh
         ch1 = 1.
         sh1 = 0.
         dch1 = cos(arg1)
         dsh1 = sin(arg1)
         do 110 j=1,ipm
            ch1h = dch1*ch1-dsh1*sh1
            sh1 = dch1*sh1+dsh1*ch1
            ch1 = ch1h
            i = is+2
            wa(i-1) = ch1
            wa(i) = sh1
            if (ido .lt. 5) go to 109
            do 108 ii=5,ido,2
               i = i+2
               wa(i-1) = ch1*wa(i-3)-sh1*wa(i-2)
               wa(i) = ch1*wa(i-2)+sh1*wa(i-3)
  108       continue
  109       is = is+ido
  110    continue
         l1 = l2
  111 continue
      return
      end
