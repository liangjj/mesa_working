      subroutine cscal(n,ca,cx,incx)
c***begin prologue  cscal
c***date written   791001   (yymmdd)
c***revision date  820801   (yymmdd)
c***category no.  d1a6
c***keywords  blas,complex,linear algebra,scale,vector
c***author  lawson, c. l., (jpl)
c           hanson, r. j., (snla)
c           kincaid, d. r., (u. of texas)
c           krogh, f. t., (jpl)
c***purpose  complex vector scale x = a*x
c***description
c
c                b l a s  subprogram
c    description of parameters
c
c     --input--
c        n  number of elements in input vector(s)
c       ca  complex scale factor
c       cx  complex vector with n elements
c     incx  storage spacing between elements of cx
c
c     --output--
c    cscal  complex result (unchanged if n .le. 0)
c
c     replace complex cx by complex ca*cx.
c     for i = 0 to n-1, replace cx(1+i*incx) with  ca * cx(1+i*incx)
c***references  lawson c.l., hanson r.j., kincaid d.r., krogh f.t.,
c                 *basic linear algebra subprograms for fortran usage*,
c                 algorithm no. 539, transactions on mathematical
c                 software, volume 5, number 3, september 1979, 308-323
c***routines called  (none)
c***end prologue  cscal
c
      implicit real *8 (a-h,o-z)
      complex*16 ca,cx(*)
c***first executable statement  cscal
#ifdef LOCALBLAS
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      do 10 i = 1,n
        cx(ix) = ca*cx(ix)
        ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 do 30 i = 1,n
        cx(i) = ca*cx(i)
   30 continue
      return
      end
#endif LOCALBLAS
#ifdef MACHINEBLAS
      call zscal(n,ca,cx,incx)
      return
      end
#endif MACHINEBLAS
