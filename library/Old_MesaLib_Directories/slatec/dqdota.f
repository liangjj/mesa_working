*deck dqdota
      double precision function dqdota (n, db, qc, dx, incx, dy, incy)
c***begin prologue  dqdota
c***purpose  compute the inner product of two vectors with extended
c            precision accumulation and result.
c***library   slatec
c***category  d1a4
c***type      double precision (dqdota-d)
c***keywords  dot product, inner product
c***author  lawson, c. l., (jpl)
c           hanson, r. j., (snla)
c           kincaid, d. r., (u. of texas)
c           krogh, f. t., (jpl)
c***description
c
c               b l a s  subprogram
c   description of parameters
c
c    --input--
c       n  number of elements in input vector(s)
c      db  double precision scalar to be added to inner product
c      qc  extended precision scalar to be added to inner product
c      dx  double precision vector with n elements
c    incx  storage spacing between elements of dx
c      dy  double precision vector with n elements
c    incy  storage spacing between elements of dy
c
c    --output--
c  dqdota  double precision result
c      qc  extended precision result
c
c    d.p. dot product with extended precision accumulation (and result)
c    qc and dqdota are set = db + qc + sum for i = 0 to n-1 of
c      dx(lx+i*incx) * dy(ly+i*incy),  where qc is an extended
c      precision result previously computed by dqdoti or dqdota
c      and lx = 1 if incx .ge. 0, else lx = (-incx)*n, and ly is
c      defined in a similar way using incy.  the mp package by
c      richard p. brent is used for the extended precision arithmetic.
c
c    fred t. krogh,  jpl,  1977,  june 1
c
c    the common block for the mp package is name mpcom.  if local
c    variable i1 is zero, dqdota calls mpblas to initialize
c    the mp package and reset i1 to 1.
c
c    the argument qc(*) and the local variables qx and qy are integer
c    arrays of size 30.  see the comments in the routine mpblas for the
c    reason for this choice.
c
c***references  c. l. lawson, r. j. hanson, d. r. kincaid and f. t.
c                 krogh, basic linear algebra subprograms for fortran
c                 usage, algorithm no. 539, transactions on mathematical
c                 software 5, 3 (september 1979), pp. 308-323.
c***routines called  mpadd, mpblas, mpcdm, mpcmd, mpmul
c***common blocks    mpcom
c***revision history  (yymmdd)
c   791001  date written
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c   930124  increased array sizes for sun -r8.  (rwc)
c***end prologue  dqdota
      double precision dx(*), dy(*), db
      integer  qc(30), qx(30), qy(30)
      common /mpcom/  mpb, mpt, mpm, mplun, mpmxr, mpr(30)
      save i1
      data  i1 / 0 /
c***first executable statement  dqdota
      if (i1 .eq. 0) call mpblas(i1)
      if (db .eq. 0.d0) go to 20
      call mpcdm(db, qx)
      call mpadd(qc, qx, qc)
   20 if (n .eq. 0) go to 40
      ix = 1
      iy = 1
      if (incx .lt. 0) ix = (-n + 1) * incx + 1
      if (incy .lt. 0) iy = (-n + 1) * incy + 1
      do  30  i = 1,n
         call mpcdm(dx(ix), qx)
         call mpcdm(dy(iy), qy)
         call mpmul(qx, qy, qx)
         call mpadd(qc, qx, qc)
         ix = ix + incx
         iy = iy + incy
   30 continue
   40 call mpcmd(qc, dqdota)
      return
      end
