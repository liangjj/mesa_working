*deck cdotu
      function cdotu(n,cx,incx,cy,incy)
c***begin prologue  cdotu
c***date written   791001   (yymmdd)
c***revision date  820801   (yymmdd)
c***category no.  d1a4
c***keywords  blas,complex,inner product,linear algebra,vector
c***author  lawson, c. l., (jpl)
c           hanson, r. j., (snla)
c           kincaid, d. r., (u. of texas)
c           krogh, f. t., (jpl)
c***purpose  inner product of complex vectors
c***description
c
c                b l a s  subprogram
c    description of parameters
c
c     --input--
c        n  number of elements in input vector(s)
c       cx  complex vector with n elements
c     incx  storage spacing between elements of cx
c       cy  complex vector with n elements
c     incy  storage spacing between elements of cy
c
c     --output--
c    cdotu  complex result (zero if n .le. 0)
c
c     returns the dot product for complex cx and cy, no conjugation
c     cdotu = sum for i = 0 to n-1 of  cx(lx+i*incx) * cy(ly+i*incy),
c     where lx = 1 if incx .ge. 0, else lx = (-incx)*n, and ly is
c     defined in a similar way using incy.
c***references  lawson c.l., hanson r.j., kincaid d.r., krogh f.t.,
c                 *basic linear algebra subprograms for fortran usage*,
c                 algorithm no. 539, transactions on mathematical
c                 software, volume 5, number 3, september 1979, 308-323
c***routines called  (none)
c***end prologue  cdotu
c
      implicit integer (a-z)
      complex*16 cdotu, cx(*), cy(*)
c***first executable statement  cdotu
#ifdef LOCALBLAS
      complex*16 temp
      temp = (0.0d0,0.0d0)
      cdotu = (0.0d0,0.0d0)
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        temp = temp + cx(ix)*cy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      cdotu = temp
      return
c
c        code for both increments equal to 1
c
   20 do 30 i = 1,n
        temp = temp + cx(i)*cy(i)
   30 continue
      cdotu = temp
      return
      end
#endif LOCALBLAS
#ifdef MACHINEBLAS
      complex*16 zdotu
      cdotu = zdotu(n,cx,incx,cy,incy)
      return
      end
#endif MACHINEBLAS
