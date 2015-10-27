*deck srotg
      subroutine srotg (sa, sb, sc, ss)
c***begin prologue  srotg
c***purpose  construct a plane givens rotation.
c***library   slatec (blas)
c***category  d1b10
c***type      single precision (srotg-s, drotg-d, crotg-c)
c***keywords  blas, givens rotation, givens transformation,
c             linear algebra, vector
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
c       sa  single precision scalar
c       sb  single precision scalar
c
c     --output--
c       sa  single precision result r
c       sb  single precision result z
c       sc  single precision result
c       ss  single precision result
c
c     construct the givens transformation
c
c         ( sc  ss )
c     g = (        ) ,    sc**2 + ss**2 = 1 ,
c         (-ss  sc )
c
c     which zeros the second entry of the 2-vector  (sa,sb)**t.
c
c     the quantity r = (+/-)sqrt(sa**2 + sb**2) overwrites sa in
c     storage.  the value of sb is overwritten by a value z which
c     allows sc and ss to be recovered by the following algorithm:
c
c           if z=1  set  sc=0.0  and  ss=1.0
c           if abs(z) .lt. 1  set  sc=sqrt(1-z**2)  and  ss=z
c           if abs(z) .gt. 1  set  sc=1/z  and  ss=sqrt(1-sc**2)
c
c     normally, the subprogram srot(n,sx,incx,sy,incy,sc,ss) will
c     next be called to apply the transformation to a 2 by n matrix.
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
c   920501  reformatted the references section.  (wrb)
c***end prologue  srotg
c***first executable statement  srotg
      if (abs(sa) .le. abs(sb)) go to 10
c
c *** here abs(sa) .gt. abs(sb) ***
c
      u = sa + sa
      v = sb / u
c
c     note that u and r have the sign of sa
c
      r = sqrt(0.25e0 + v**2) * u
c
c     note that sc is positive
c
      sc = sa / r
      ss = v * (sc + sc)
      sb = ss
      sa = r
      return
c
c *** here abs(sa) .le. abs(sb) ***
c
   10 if (sb .eq. 0.0e0) go to 20
      u = sb + sb
      v = sa / u
c
c     note that u and r have the sign of sb
c     (r is immediately stored in sa)
c
      sa = sqrt(0.25e0 + v**2) * u
c
c     note that ss is positive
c
      ss = sb / sa
      sc = v * (ss + ss)
      if (sc .eq. 0.0e0) go to 15
      sb = 1.0e0 / sc
      return
   15 sb = 1.0e0
      return
c
c *** here sa = sb = 0.0 ***
c
   20 sc = 1.0e0
      ss = 0.0e0
      return
c
      end
