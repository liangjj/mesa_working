*deck dcdot
      subroutine dcdot (n, fm, cx, incx, cy, incy, dcr, dci)
c***begin prologue  dcdot
c***purpose  compute the inner product of two vectors with extended
c            precision accumulation and result.
c***library   slatec (blas)
c***category  d1a4
c***type      complex (dsdot-d, dcdot-c)
c***keywords  blas, complex vectors, dot product, inner product,
c             linear algebra, vector
c***author  (unknown)
c***description
c
c    compute the dot product of 2 complex vectors, cx and cy, e.g.
c    cx dot cy, or, cxconjugate dot cy.  the real and imaginary
c    parts of cx and cy are converted to double precision, the dot
c    product accumulation is done in double precision and the output
c    is given as 2 double precision numbers, corresponding to the real
c    and imaginary part of the result.
c     input
c      n:  number of complex components of cx and cy.
c      fm: =+1.0   compute cx dot cy.
c          =-1.0   compute cxconjugate dot cy.
c      cx(n):
c      cy(n):  complex arrays of length n.
c      incx:(integer)   spacing of elements of cx to use
c      incy:(integer)   spacing of elements of cy to use.
c     output
c      dcr:(double precision) real part of dot product.
c      dci:(double precision) imaginary part of dot product.
c
c***references  (none)
c***routines called  (none)
c***revision history  (yymmdd)
c   790101  date written
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c***end prologue  dcdot
      integer i, incx, incy, kx, ky, n
      complex cx(*), cy(*)
      double precision dcr, dci, dt1, dt2, dt3, dt4, fm
c***first executable statement  dcdot
      dcr = 0.0d0
      dci = 0.0d0
      if (n .le. 0) go to 20
c
      kx = 1
      ky = 1
      if (incx .lt. 0) kx = 1+(1-n)*incx
      if (incy .lt. 0) ky = 1+(1-n)*incy
      do 10 i = 1,n
        dt1 = dble(real(cx(kx)))
        dt2 = dble(real(cy(ky)))
        dt3 = dble(aimag(cx(kx)))
        dt4 = dble(aimag(cy(ky)))
        dcr = dcr+(dt1*dt2)-fm*(dt3*dt4)
        dci = dci+(dt1*dt4)+fm*(dt3*dt2)
        kx = kx+incx
        ky = ky+incy
   10 continue
   20 return
      end
