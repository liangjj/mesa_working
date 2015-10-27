*deck cosqb1
      subroutine cosqb1 (n, x, w, xh)
c***begin prologue  cosqb1
c***subsidiary
c***purpose  compute the unnormalized inverse of cosqf1.
c***library   slatec (fftpack)
c***category  j1a3
c***type      single precision (cosqb1-s)
c***keywords  fftpack, fourier transform
c***author  swarztrauber, p. n., (ncar)
c***description
c
c  subroutine cosqb1 computes the fast fourier transform of quarter
c  wave data. that is, cosqb1 computes a sequence from its
c  representation in terms of a cosine series with odd wave numbers.
c  the transform is defined below at output parameter x.
c
c***references  p. n. swarztrauber, vectorizing the ffts, in parallel
c                 computations (g. rodrigue, ed.), academic press,
c                 1982, pp. 51-83.
c***routines called  rfftb
c***revision history  (yymmdd)
c   790601  date written
c   830401  modified to use slatec library source file format.
c   860115  modified by ron boisvert to adhere to fortran 77 by
c           changing dummy array size declarations (1) to (*).
c   881128  modified by dick valent to meet prologue standards.
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  cosqb1
      dimension x(*), w(*), xh(*)
c***first executable statement  cosqb1
      ns2 = (n+1)/2
      np2 = n+2
      do 101 i=3,n,2
         xim1 = x(i-1)+x(i)
         x(i) = x(i)-x(i-1)
         x(i-1) = xim1
  101 continue
      x(1) = x(1)+x(1)
      modn = mod(n,2)
      if (modn .eq. 0) x(n) = x(n)+x(n)
      call rfftb (n,x,xh)
      do 102 k=2,ns2
         kc = np2-k
         xh(k) = w(k-1)*x(kc)+w(kc-1)*x(k)
         xh(kc) = w(k-1)*x(k)-w(kc-1)*x(kc)
  102 continue
      if (modn .eq. 0) x(ns2+1) = w(ns2)*(x(ns2+1)+x(ns2+1))
      do 103 k=2,ns2
         kc = np2-k
         x(k) = xh(k)+xh(kc)
         x(kc) = xh(k)-xh(kc)
  103 continue
      x(1) = x(1)+x(1)
      return
      end
