*deck drotm
      subroutine drotm (n, dx, incx, dy, incy, dparam)
c***begin prologue  drotm
c***purpose  apply a modified givens transformation.
c***library   slatec (blas)
c***category  d1a8
c***type      double precision (srotm-s, drotm-d)
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
c       dx  double precision vector with n elements
c     incx  storage spacing between elements of dx
c       dy  double precision vector with n elements
c     incy  storage spacing between elements of dy
c   dparam  5-element d.p. vector.  dparam(1) is dflag described below.
c           locations 2-5 of sparam contain elements of the
c           transformation matrix h described below.
c
c     --output--
c       dx  rotated vector (unchanged if n .le. 0)
c       dy  rotated vector (unchanged if n .le. 0)
c
c     apply the modified givens transformation, h, to the 2 by n matrix
c     (dx**t)
c     (dy**t) , where **t indicates transpose.  the elements of dx are
c     in dx(lx+i*incx), i = 0 to n-1, where lx = 1 if incx .ge. 0, else
c     lx = 1+(1-n)*incx, and similarly for dy using ly and incy.
c
c     with dparam(1)=dflag, h has one of the following forms:
c
c     dflag=-1.d0     dflag=0.d0        dflag=1.d0     dflag=-2.d0
c
c       (dh11  dh12)    (1.d0  dh12)    (dh11  1.d0)    (1.d0  0.d0)
c     h=(          )    (          )    (          )    (          )
c       (dh21  dh22),   (dh21  1.d0),   (-1.d0 dh22),   (0.d0  1.d0).
c
c     see drotmg for a description of data storage in dparam.
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
c***end prologue  drotm
      double precision dflag, dh12, dh22, dx, two, z, dh11, dh21,
     1                 dparam, dy, w, zero
      dimension dx(*), dy(*), dparam(5)
      save zero, two
      data zero, two /0.0d0, 2.0d0/
c***first executable statement  drotm
      dflag=dparam(1)
      if (n.le.0 .or. (dflag+two.eq.zero)) go to 140
          if (.not.(incx.eq.incy.and. incx .gt.0)) go to 70
c
               nsteps=n*incx
               if (dflag) 50,10,30
   10          continue
               dh12=dparam(4)
               dh21=dparam(3)
                    do 20 i = 1,nsteps,incx
                    w=dx(i)
                    z=dy(i)
                    dx(i)=w+z*dh12
                    dy(i)=w*dh21+z
   20               continue
               go to 140
   30          continue
               dh11=dparam(2)
               dh22=dparam(5)
                    do 40 i = 1,nsteps,incx
                    w=dx(i)
                    z=dy(i)
                    dx(i)=w*dh11+z
                    dy(i)=-w+dh22*z
   40               continue
               go to 140
   50          continue
               dh11=dparam(2)
               dh12=dparam(4)
               dh21=dparam(3)
               dh22=dparam(5)
                    do 60 i = 1,nsteps,incx
                    w=dx(i)
                    z=dy(i)
                    dx(i)=w*dh11+z*dh12
                    dy(i)=w*dh21+z*dh22
   60               continue
               go to 140
   70     continue
          kx=1
          ky=1
          if (incx .lt. 0) kx = 1+(1-n)*incx
          if (incy .lt. 0) ky = 1+(1-n)*incy
c
          if (dflag) 120,80,100
   80     continue
          dh12=dparam(4)
          dh21=dparam(3)
               do 90 i = 1,n
               w=dx(kx)
               z=dy(ky)
               dx(kx)=w+z*dh12
               dy(ky)=w*dh21+z
               kx=kx+incx
               ky=ky+incy
   90          continue
          go to 140
  100     continue
          dh11=dparam(2)
          dh22=dparam(5)
               do 110 i = 1,n
               w=dx(kx)
               z=dy(ky)
               dx(kx)=w*dh11+z
               dy(ky)=-w+dh22*z
               kx=kx+incx
               ky=ky+incy
  110          continue
          go to 140
  120     continue
          dh11=dparam(2)
          dh12=dparam(4)
          dh21=dparam(3)
          dh22=dparam(5)
               do 130 i = 1,n
               w=dx(kx)
               z=dy(ky)
               dx(kx)=w*dh11+z*dh12
               dy(ky)=w*dh21+z*dh22
               kx=kx+incx
               ky=ky+incy
  130          continue
  140     continue
          return
      end
