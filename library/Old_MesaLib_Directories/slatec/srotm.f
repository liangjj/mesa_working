*deck srotm
      subroutine srotm (n, sx, incx, sy, incy, sparam)
c***begin prologue  srotm
c***purpose  apply a modified givens transformation.
c***library   slatec (blas)
c***category  d1a8
c***type      single precision (srotm-s, drotm-d)
c***keywords  blas, linear algebra, modified givens rotation, vector
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
c   sparam  5-element vector. sparam(1) is sflag described below.
c           locations 2-5 of sparam contain elements of the
c           transformation matrix h described below.
c
c     --output--
c       sx  rotated vector (unchanged if n .le. 0)
c       sy  rotated vector (unchanged if n .le. 0)
c
c     apply the modified givens transformation, h, to the 2 by n matrix
c     (sx**t)
c     (sy**t) , where **t indicates transpose.  the elements of sx are
c     in sx(lx+i*incx), i = 0 to n-1, where lx = 1 if incx .ge. 0, else
c     lx = 1+(1-n)*incx, and similarly for sy using ly and incy.
c
c     with sparam(1)=sflag, h has one of the following forms:
c
c     sflag=-1.e0     sflag=0.e0        sflag=1.e0     sflag=-2.e0
c
c       (sh11  sh12)    (1.e0  sh12)    (sh11  1.e0)    (1.e0  0.e0)
c     h=(          )    (          )    (          )    (          )
c       (sh21  sh22),   (sh21  1.e0),   (-1.e0 sh22),   (0.e0  1.e0).
c
c     see srotmg for a description of data storage in sparam.
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
c***end prologue  srotm
      dimension sx(*), sy(*), sparam(5)
      save zero, two
      data zero, two /0.0e0, 2.0e0/
c***first executable statement  srotm
      sflag=sparam(1)
      if (n.le.0 .or. (sflag+two.eq.zero)) go to 140
          if (.not.(incx.eq.incy.and. incx .gt.0)) go to 70
c
               nsteps=n*incx
               if (sflag) 50,10,30
   10          continue
               sh12=sparam(4)
               sh21=sparam(3)
                    do 20 i = 1,nsteps,incx
                    w=sx(i)
                    z=sy(i)
                    sx(i)=w+z*sh12
                    sy(i)=w*sh21+z
   20               continue
               go to 140
   30          continue
               sh11=sparam(2)
               sh22=sparam(5)
                    do 40 i = 1,nsteps,incx
                    w=sx(i)
                    z=sy(i)
                    sx(i)=w*sh11+z
                    sy(i)=-w+sh22*z
   40               continue
               go to 140
   50          continue
               sh11=sparam(2)
               sh12=sparam(4)
               sh21=sparam(3)
               sh22=sparam(5)
                    do 60 i = 1,nsteps,incx
                    w=sx(i)
                    z=sy(i)
                    sx(i)=w*sh11+z*sh12
                    sy(i)=w*sh21+z*sh22
   60               continue
               go to 140
   70     continue
          kx=1
          ky=1
          if (incx .lt. 0) kx = 1+(1-n)*incx
          if (incy .lt. 0) ky = 1+(1-n)*incy
c
          if (sflag) 120,80,100
   80     continue
          sh12=sparam(4)
          sh21=sparam(3)
               do 90 i = 1,n
               w=sx(kx)
               z=sy(ky)
               sx(kx)=w+z*sh12
               sy(ky)=w*sh21+z
               kx=kx+incx
               ky=ky+incy
   90          continue
          go to 140
  100     continue
          sh11=sparam(2)
          sh22=sparam(5)
               do 110 i = 1,n
               w=sx(kx)
               z=sy(ky)
               sx(kx)=w*sh11+z
               sy(ky)=-w+sh22*z
               kx=kx+incx
               ky=ky+incy
  110          continue
          go to 140
  120     continue
          sh11=sparam(2)
          sh12=sparam(4)
          sh21=sparam(3)
          sh22=sparam(5)
               do 130 i = 1,n
               w=sx(kx)
               z=sy(ky)
               sx(kx)=w*sh11+z*sh12
               sy(ky)=w*sh21+z*sh22
               kx=kx+incx
               ky=ky+incy
  130          continue
  140     continue
          return
      end
