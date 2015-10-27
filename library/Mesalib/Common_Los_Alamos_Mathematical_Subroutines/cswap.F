*deck cswap
      subroutine cswap(n,cx,incx,cy,incy)
c***begin prologue  cswap
c***date written   791001   (yymmdd)
c***revision date  820801   (yymmdd)
c***category no.  d1a5
c***keywords  blas,complex,interchange,linear algebra,vector
c***author  lawson, c. l., (jpl)
c           hanson, r. j., (snla)
c           kincaid, d. r., (u. of texas)
c           krogh, f. t., (jpl)
c***purpose  interchange complex vectors
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
c       cx  input vector cy (unchanged if n .le. 0)
c       cy  input vector cx (unchanged if n .le. 0)
c
c     interchange complex cx and complex cy
c     for i = 0 to n-1, interchange  cx(lx+i*incx) and cy(ly+i*incy),
c     where lx = 1 if incx .gt. 0, else lx = (-incx)*n, and ly is
c     defined in a similar way using incy.
c***references  lawson c.l., hanson r.j., kincaid d.r., krogh f.t.,
c                 *basic linear algebra subprograms for fortran usage*,
c                 algorithm no. 539, transactions on mathematical
c                 software, volume 5, number 3, september 1979, 308-323
c***routines called  (none)
c***end prologue  cswap
c
      complex*16 cx(*),cy(*)
c***first executable statement  cswap
#ifdef LOCALBLAS
      complex*16 temp
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        temp = cx(ix)
        cx(ix) = cy(iy)
        cy(iy) = temp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1
   20 do 30 i = 1,n
        temp = cx(i)
        cx(i) = cy(i)
        cy(i) = temp
   30 continue
      return
      end
#endif LOCALBLAS
#ifdef MACHINEBLAS
      call zswap(n,cx,incx,cy,incy)
      return
      end
#endif MACHINEBLAS
