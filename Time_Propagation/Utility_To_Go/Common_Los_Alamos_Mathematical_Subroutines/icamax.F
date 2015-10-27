*deck icamax
      function icamax(n,cx,incx)
c***begin prologue  icamax
c***date written   791001   (yymmdd)
c***revision date  860401   (yymmdd)
c***category no.  d1a2
c***keywords  blas,complex,linear algebra,maximum component,vector
c***author  lawson, c. l., (jpl)
c           hanson, r. j., (snla)
c           kincaid, d. r., (u. of texas)
c           krogh, f. t., (jpl)
c***purpose  find the smallest index of an element of a complex vector
c            having the maximum magnitude of the sums of the real and
c            imaginary components.
c***description
c
c                b l a s  subprogram
c    description of parameters
c
c     --input--
c        n  number of elements in input vector(s)
c       cx  complex vector with n elements
c     incx  storage spacing between elements of cx
c
c     --output--
c   icamax  smallest index (zero if n .le. 0)
c
c      returns the index of the component of cx having the
c      largest sum of magnitudes of real and imaginary parts.
c     icamax = first i, i = 1 to n, to minimize
c        abs(real(cx(1-incx+i*incx))) + abs(imag(cx(1-incx+i*incx)))
c
c***references  lawson c.l., hanson r.j., kincaid d.r., krogh f.t.,
c                 *basic linear algebra subprograms for fortran usage*,
c                 algorithm no. 539, transactions on mathematical
c                 software, volume 5, number 3, september 1979, 308-323
c***routines called  (none)
c***end prologue  icamax
c
      implicit real *8 (a-h,o-z)
      integer icamax
      complex*16 cx(*)
c***first executable statement  icamax
#ifdef LOCALBLAS
      icamax = 0
      if( n.lt.1 .or. incx.le.0 )return
      icamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      smax = dcabs1(cx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dcabs1(cx(ix)).le.smax) go to 5
         icamax = i
         smax = dcabs1(cx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 smax = dcabs1(cx(1))
      do 30 i = 2,n
         if(dcabs1(cx(i)).le.smax) go to 30
         icamax = i
         smax = dcabs1(cx(i))
   30 continue
      return
      end
#endif LOCALBLAS
#ifdef MACHINEBLAS
      integer izamax
      icamax=izamax(n,cx,incx)
      return
      end
#endif MACHINEBLAS
