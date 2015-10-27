*deck cfftb
      subroutine cfftb (n, c, wsave)
c***begin prologue  cfftb
c***subsidiary
c***purpose  compute the unnormalized inverse of cfftf.
c***library   slatec (fftpack)
c***category  j1a2
c***type      complex (rfftb-s, cfftb-c)
c***keywords  fftpack, fourier transform
c***author  swarztrauber, p. n., (ncar)
c***description
c
c  ********************************************************************
c  *   notice   notice   notice   notice   notice   notice   notice   *
c  ********************************************************************
c  *                                                                  *
c  *   this routine uses non-standard fortran 77 constructs and will  *
c  *   be removed from the library at a future date.  you are         *
c  *   requested to use cfftb1.                                       *
c  *                                                                  *
c  ********************************************************************
c
c  subroutine cfftb computes the backward complex discrete fourier
c  transform (the fourier synthesis).  equivalently, cfftb computes
c  a complex periodic sequence from its fourier coefficients.
c  the transform is defined below at output parameter c.
c
c  a call of cfftf followed by a call of cfftb will multiply the
c  sequence by n.
c
c  the array wsave which is used by subroutine cfftb must be
c  initialized by calling subroutine cffti(n,wsave).
c
c  input parameters
c
c  n       the length of the complex sequence c.  the method is
c          more efficient when n is the product of small primes.
c
c  c       a complex array of length n which contains the sequence
c
c  wsave   a real work array which must be dimensioned at least 4*n+15
c          in the program that calls cfftb.  the wsave array must be
c          initialized by calling subroutine cffti(n,wsave), and a
c          different wsave array must be used for each different
c          value of n.  this initialization does not have to be
c          repeated so long as n remains unchanged.  thus subsequent
c          transforms can be obtained faster than the first.
c          the same wsave array can be used by cfftf and cfftb.
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
c  wsave   contains initialization calculations which must not be
c          destroyed between calls of subroutine cfftf or cfftb
c
c***references  p. n. swarztrauber, vectorizing the ffts, in parallel
c                 computations (g. rodrigue, ed.), academic press,
c                 1982, pp. 51-83.
c***routines called  cfftb1
c***revision history  (yymmdd)
c   790601  date written
c   830401  modified to use slatec library source file format.
c   860115  modified by ron boisvert to adhere to fortran 77 by
c           changing dummy array size declarations (1) to (*).
c   861211  revision date from version 3.2
c   881128  modified by dick valent to meet prologue standards.
c   891214  prologue converted to version 4.0 format.  (bab)
c   900131  routine changed from user-callable to subsidiary
c           because of non-standard fortran 77 arguments in the
c           call to cfftb1.  (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  cfftb
      complex c
      dimension c(*), wsave(*)
c***first executable statement  cfftb
      if (n .eq. 1) return
      iw1 = n+n+1
      iw2 = iw1+n+n
      call cfftb1 (n,c,wsave,wsave(iw1),wsave(iw2))
      return
      end
