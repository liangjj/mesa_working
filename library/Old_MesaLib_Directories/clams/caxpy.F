*deck caxpy
      subroutine caxpy(n,ca,cx,incx,cy,incy)
c***begin prologue  caxpy
c***date written   791001   (yymmdd)
c***revision date  820801   (yymmdd)
c***category no.  d1a7
c***keywords  blas,complex,linear algebra,triad,vector
c***author  lawson, c. l., (jpl)
c           hanson, r. j., (snla)
c           kincaid, d. r., (u. of texas)
c           krogh, f. t., (jpl)
c***purpose  complex computation y = a*x + y
c***description
c
c                b l a s  subprogram
c    description of parameters
c
c     --input--
c        n  number of elements in input vector(s)
c       ca  complex scalar multiplier
c       cx  complex vector with n elements
c     incx  storage spacing between elements of cx
c       cy  complex vector with n elements
c     incy  storage spacing between elements of cy
c
c     --output--
c       cy  complex result (unchanged if n .le. 0)
c
c     overwrite complex cy with complex  ca*cx + cy.
c     for i = 0 to n-1, replace  cy(ly+i*incy) with ca*cx(lx+i*incx) +
c       cy(ly+i*incy), where lx = 1 if incx .ge. 0, else lx = (-incx)*n
c       and ly is defined in a similar way using incy.
c***references  lawson c.l., hanson r.j., kincaid d.r., krogh f.t.,
c                 *basic linear algebra subprograms for fortran usage*,
c                 algorithm no. 539, transactions on mathematical
c                 software, volume 5, number 3, september 1979, 308-323
c***routines called  (none)
c***end prologue  caxpy
c
      implicit real *8 (a-h,o-z)
      complex*16 cx(*),cy(*),ca
c***first executable statement  caxpy
#ifdef LOCALBLAS
      if(n.le.0)return
      if (dcabs1(ca) .eq. 0.0d0) return
      if (incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        cy(iy) = cy(iy) + ca*cx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
   20 do 30 i = 1,n
        cy(i) = cy(i) + ca*cx(i)
   30 continue
      return
      end
#endif LOCALBLAS
#ifdef MACHINEBLAS
      call zaxpy(n,ca,cx,incx,cy,incy)
      return
      end
#endif MACHINEBLAS
