*deck cosqf1
      subroutine cosqf1 (n, x, w, xh)
c***begin prologue  cosqf1
c***subsidiary
c***purpose  compute the forward cosine transform with odd wave numbers.
c***library   slatec (fftpack)
c***category  j1a3
c***type      single precision (cosqf1-s)
c***keywords  fftpack, fourier transform
c***author  swarztrauber, p. n., (ncar)
c***description
c
c  subroutine cosqf1 computes the fast fourier transform of quarter
c  wave data. that is, cosqf1 computes the coefficients in a cosine
c  series representation with only odd wave numbers.  the transform
c  is defined below at output parameter x
c
c***references  p. n. swarztrauber, vectorizing the ffts, in parallel
c                 computations (g. rodrigue, ed.), academic press,
c                 1982, pp. 51-83.
c***routines called  rfftf
c***revision history  (yymmdd)
c   790601  date written
c   830401  modified to use slatec library source file format.
c   860115  modified by ron boisvert to adhere to fortran 77 by
c           changing dummy array size declarations (1) to (*).
c   881128  modified by dick valent to meet prologue standards.
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  cosqf1
      dimension x(*), w(*), xh(*)
c***first executable statement  cosqf1
      ns2 = (n+1)/2
      np2 = n+2
      do 101 k=2,ns2
         kc = np2-k
         xh(k) = x(k)+x(kc)
         xh(kc) = x(k)-x(kc)
  101 continue
      modn = mod(n,2)
      if (modn .eq. 0) xh(ns2+1) = x(ns2+1)+x(ns2+1)
      do 102 k=2,ns2
         kc = np2-k
         x(k) = w(k-1)*xh(kc)+w(kc-1)*xh(k)
         x(kc) = w(k-1)*xh(k)-w(kc-1)*xh(kc)
  102 continue
      if (modn .eq. 0) x(ns2+1) = w(ns2)*xh(ns2+1)
      call rfftf (n,x,xh)
      do 103 i=3,n,2
         xim1 = x(i-1)-x(i)
         x(i) = x(i-1)+x(i)
         x(i-1) = xim1
  103 continue
      return
      end
