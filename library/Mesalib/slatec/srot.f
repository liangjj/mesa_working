*deck srot
      subroutine srot (n, sx, incx, sy, incy, sc, ss)
c***begin prologue  srot
c***purpose  apply a plane givens rotation.
c***library   slatec (blas)
c***category  d1a8
c***type      single precision (srot-s, drot-d, csrot-c)
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
c       sx  single precision vector with n elements
c     incx  storage spacing between elements of sx
c       sy  single precision vector with n elements
c     incy  storage spacing between elements of sy
c       sc  element of rotation matrix
c       ss  element of rotation matrix
c
c     --output--
c       sx  rotated vector sx (unchanged if n .le. 0)
c       sy  rotated vector sy (unchanged if n .le. 0)
c
c     multiply the 2 x 2 matrix  ( sc ss) times the 2 x n matrix (sx**t)
c                                (-ss sc)                        (sy**t)
c     where **t indicates transpose.  the elements of sx are in
c     sx(lx+i*incx), i = 0 to n-1, where lx = 1 if incx .ge. 0, else
c     lx = 1+(1-n)*incx, and similarly for sy using ly and incy.
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
c***end prologue  srot
      real sx, sy, sc, ss, zero, one, w, z
      dimension sx(*), sy(*)
      save zero, one
      data zero, one /0.0e0, 1.0e0/
c***first executable statement  srot
      if (n .le. 0 .or. (ss .eq. zero .and. sc .eq. one)) go to 40
      if (.not. (incx .eq. incy .and. incx .gt. 0)) go to 20
c
c          code for equal and positive increments.
c
           nsteps=incx*n
           do 10 i = 1,nsteps,incx
                w=sx(i)
                z=sy(i)
                sx(i)=sc*w+ss*z
                sy(i)=-ss*w+sc*z
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
                w=sx(kx)
                z=sy(ky)
                sx(kx)=sc*w+ss*z
                sy(ky)=-ss*w+sc*z
                kx=kx+incx
                ky=ky+incy
   30           continue
   40 continue
c
      return
      end
