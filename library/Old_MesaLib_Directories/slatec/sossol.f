*deck sossol
      subroutine sossol (k, n, l, x, c, b, m)
c***begin prologue  sossol
c***subsidiary
c***purpose  subsidiary to sos
c***library   slatec
c***type      single precision (sossol-s, dsossl-d)
c***author  (unknown)
c***description
c
c     sossol solves an upper triangular type of linear system by back
c     substitution.
c
c     the matrix c is upper trapezoidal and stored as a linear array by
c     rows. the equations have been normalized so that the diagonal
c     entries of c are understood to be unity. the off diagonal entries
c     and the elements of the constant right hand side vector b have
c     already been stored as the negatives of the corresponding equation
c     values.
c     with each call to sossol a (k-1) by (k-1) triangular system is
c     resolved. for l greater than k, column l of c is included in the
c     right hand side vector.
c
c***see also  sos
c***routines called  (none)
c***revision history  (yymmdd)
c   801001  date written
c   890831  modified array declarations.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c***end prologue  sossol
c
c
      dimension x(*), c(*), b(*)
c
c***first executable statement  sossol
      np1 = n + 1
      km1 = k - 1
      lk = km1
      if (l .eq. k) lk = k
      kn = m
c
c
      do 40 kj=1,km1
        kmm1 = k - kj
        km = kmm1 + 1
        xmax = 0.
        kn = kn - np1 + kmm1
        if (km .gt. lk) go to 20
        jkm = kn
c
        do 10 j=km,lk
          jkm = jkm + 1
          xmax = xmax + c(jkm)*x(j)
   10   continue
c
   20   if (l .le. k) go to 30
        jkm = kn + l - kmm1
        xmax = xmax + c(jkm)*x(l)
   30   x(kmm1) = xmax + b(kmm1)
   40 continue
c
      return
      end
