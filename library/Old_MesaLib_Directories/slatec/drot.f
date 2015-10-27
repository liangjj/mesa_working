*deck drot
      subroutine drot (n, dx, incx, dy, incy, dc, ds)
c***begin prologue  drot
c***purpose  apply a plane givens rotation.
c***library   slatec (blas)
c***category  d1a8
c***type      double precision (srot-s, drot-d, csrot-c)
c***keywords  blas, givens rotation, givens transformation,
c             linear algebra, plane rotation, vector
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
c       dy  double precision vector with n elements
c     incy  storage spacing between elements of dy
c       dc  d.p. element of rotation matrix
c       ds  d.p. element of rotation matrix
c
c     --output--
c       dx  rotated vector dx (unchanged if n .le. 0)
c       dy  rotated vector dy (unchanged if n .le. 0)
c
c     multiply the 2 x 2 matrix  ( dc ds) times the 2 x n matrix (dx**t)
c                                (-ds dc)                        (dy**t)
c     where **t indicates transpose.  the elements of dx are in
c     dx(lx+i*incx), i = 0 to n-1, where lx = 1 if incx .ge. 0, else
c     lx = 1+(1-n)*incx, and similarly for dy using ly and incy.
c
c***references  c. l. lawson, r. j. hanson, d. r. kincaid and f. t.
c                 krogh, basic linear algebra subprograms for fortran
c                 usage, algorithm no. 539, transactions on mathematical
c                 software 5, 3 (september 1979), pp. 308-323.
c***routines called  (none)
c***revision history  (yymmdd)
c   791001  date written
c   861211  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920310  corrected definition of lx in description.  (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  drot
      double precision dx, dy, dc, ds, zero, one, w, z
      dimension dx(*), dy(*)
      save zero, one
      data zero, one /0.0d0, 1.0d0/
c***first executable statement  drot
      if (n .le. 0 .or. (ds .eq. zero .and. dc .eq. one)) go to 40
      if (.not. (incx .eq. incy .and. incx .gt. 0)) go to 20
c
c          code for equal and positive increments.
c
           nsteps=incx*n
           do 10 i = 1,nsteps,incx
                w=dx(i)
                z=dy(i)
                dx(i)=dc*w+ds*z
                dy(i)=-ds*w+dc*z
   10           continue
           go to 40
c
c     code for unequal or nonpositive increments.
c
   20 continue
           kx=1
           ky=1
c
           if (incx .lt. 0) kx = 1-(n-1)*incx
           if (incy .lt. 0) ky = 1-(n-1)*incy
c
           do 30 i = 1,n
                w=dx(kx)
                z=dy(ky)
                dx(kx)=dc*w+ds*z
                dy(ky)=-ds*w+dc*z
                kx=kx+incx
                ky=ky+incy
   30           continue
   40 continue
c
      return
      end
