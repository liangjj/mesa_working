*deck cmpcsg
      subroutine cmpcsg (n, ijump, fnum, fden, a)
c***begin prologue  cmpcsg
c***subsidiary
c***purpose  subsidiary to cmgnbn
c***library   slatec
c***type      complex (cosgen-s, cmpcsg-c)
c***author  (unknown)
c***description
c
c     this subroutine computes required cosine values in ascending
c     order.  when ijump .gt. 1 the routine computes values
c
c        2*cos(j*pi/l) , j=1,2,...,l and j .ne. 0(mod n/ijump+1)
c
c     where l = ijump*(n/ijump+1).
c
c
c     when ijump = 1 it computes
c
c            2*cos((j-fnum)*pi/(n+fden)) ,  j=1, 2, ... ,n
c
c     where
c        fnum = 0.5, fden = 0.0,  for regular reduction values.
c        fnum = 0.0, fden = 1.0, for b-r and c-r when istag = 1
c        fnum = 0.0, fden = 0.5, for b-r and c-r when istag = 2
c        fnum = 0.5, fden = 0.5, for b-r and c-r when istag = 2
c                                in cmposn only.
c
c***see also  cmgnbn
c***routines called  pimach
c***revision history  (yymmdd)
c   801001  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c***end prologue  cmpcsg
      complex         a
      dimension       a(*)
c
c
c***first executable statement  cmpcsg
      pi = pimach(dum)
      if (n .eq. 0) go to 105
      if (ijump .eq. 1) go to 103
      k3 = n/ijump+1
      k4 = k3-1
      pibyn = pi/(n+ijump)
      do 102 k=1,ijump
         k1 = (k-1)*k3
         k5 = (k-1)*k4
         do 101 i=1,k4
            x = k1+i
            k2 = k5+i
            a(k2) = cmplx(-2.*cos(x*pibyn),0.)
  101    continue
  102 continue
      go to 105
  103 continue
      np1 = n+1
      y = pi/(n+fden)
      do 104 i=1,n
         x = np1-i-fnum
         a(i) = cmplx(2.*cos(x*y),0.)
  104 continue
  105 continue
      return
      end
