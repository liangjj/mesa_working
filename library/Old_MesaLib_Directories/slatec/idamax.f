*deck idamax
      integer function idamax (n, dx, incx)
c***begin prologue  idamax
c***purpose  find the smallest index of that component of a vector
c            having the maximum magnitude.
c***library   slatec (blas)
c***category  d1a2
c***type      double precision (isamax-s, idamax-d, icamax-c)
c***keywords  blas, linear algebra, maximum component, vector
c***author  lawson, c. l., (jpl)
c           hanson, r. j., (snla)
c           kincaid, d. r., (u. of texas)
c           krogh, f. t., (jpl)
c***description
c
c                b l a s  subprogram
c    description of parameters
c
c     --input--
c        n  number of elements in input vector(s)
c       dx  double precision vector with n elements
c     incx  storage spacing between elements of dx
c
c     --output--
c   idamax  smallest index (zero if n .le. 0)
c
c     find smallest index of maximum magnitude of double precision dx.
c     idamax = first i, i = 1 to n, to maximize abs(dx(ix+(i-1)*incx)),
c     where ix = 1 if incx .ge. 0, else ix = 1+(1-n)*incx.
c
c***references  c. l. lawson, r. j. hanson, d. r. kincaid and f. t.
c                 krogh, basic linear algebra subprograms for fortran
c                 usage, algorithm no. 539, transactions on mathematical
c                 software 5, 3 (september 1979), pp. 308-323.
c***routines called  (none)
c***revision history  (yymmdd)
c   791001  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900821  modified to correct problem with a negative increment.
c           (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  idamax
      double precision dx(*), dmax, xmag
      integer i, incx, ix, n
c***first executable statement  idamax
      idamax = 0
      if (n .le. 0) return
      idamax = 1
      if (n .eq. 1) return
c
      if (incx .eq. 1) goto 20
c
c     code for increments not equal to 1.
c
      ix = 1
      if (incx .lt. 0) ix = (-n+1)*incx + 1
      dmax = abs(dx(ix))
      ix = ix + incx
      do 10 i = 2,n
        xmag = abs(dx(ix))
        if (xmag .gt. dmax) then
          idamax = i
          dmax = xmag
        endif
        ix = ix + incx
   10 continue
      return
c
c     code for increments equal to 1.
c
   20 dmax = abs(dx(1))
      do 30 i = 2,n
        xmag = abs(dx(i))
        if (xmag .gt. dmax) then
          idamax = i
          dmax = xmag
        endif
   30 continue
      return
      end
