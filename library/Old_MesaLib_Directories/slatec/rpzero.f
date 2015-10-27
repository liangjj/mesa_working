*deck rpzero
      subroutine rpzero (n, a, r, t, iflg, s)
c***begin prologue  rpzero
c***purpose  find the zeros of a polynomial with real coefficients.
c***library   slatec
c***category  f1a1a
c***type      single precision (rpzero-s, cpzero-c)
c***keywords  polynomial roots, polynomial zeros, real roots
c***author  kahaner, d. k., (nbs)
c***description
c
c      find the zeros of the real polynomial
c         p(x)= a(1)*x**n + a(2)*x**(n-1) +...+ a(n+1)
c
c    input...
c       n = degree of p(x)
c       a = real vector containing coefficients of p(x),
c            a(i) = coefficient of x**(n+1-i)
c       r = n word complex vector containing initial estimates for zeros
c            if these are known.
c       t = 6(n+1) word array used for temporary storage
c       iflg = flag to indicate if initial estimates of
c              zeros are input.
c            if iflg .eq. 0, no estimates are input.
c            if iflg .ne. 0, the vector r contains estimates of
c               the zeros
c       ** warning ****** if estimates are input, they must
c                         be separated; that is, distinct or
c                         not repeated.
c       s = an n word array
c
c    output...
c       r(i) = ith zero,
c       s(i) = bound for r(i) .
c       iflg = error diagnostic
c    error diagnostics...
c       if iflg .eq. 0 on return, all is well.
c       if iflg .eq. 1 on return, a(1)=0.0 or n=0 on input.
c       if iflg .eq. 2 on return, the program failed to converge
c                after 25*n iterations.  best current estimates of the
c                zeros are in r(i).  error bounds are not calculated.
c
c***references  (none)
c***routines called  cpzero
c***revision history  (yymmdd)
c   810223  date written
c   890206  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c***end prologue  rpzero
c
      complex r(*), t(*)
      real a(*), s(*)
c***first executable statement  rpzero
      n1=n+1
      do 1 i=1,n1
      t(i)= cmplx(a(i),0.0)
    1 continue
      call cpzero(n,t,r,t(n+2),iflg,s)
      return
      end
