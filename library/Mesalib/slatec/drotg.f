*deck drotg
      subroutine drotg (da, db, dc, ds)
c***begin prologue  drotg
c***purpose  construct a plane givens rotation.
c***library   slatec (blas)
c***category  d1b10
c***type      double precision (srotg-s, drotg-d, crotg-c)
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
c       da  double precision scalar
c       db  double precision scalar
c
c     --output--
c       da  double precision result r
c       db  double precision result z
c       dc  double precision result
c       ds  double precision result
c
c     construct the givens transformation
c
c         ( dc  ds )
c     g = (        ) ,    dc**2 + ds**2 = 1 ,
c         (-ds  dc )
c
c     which zeros the second entry of the 2-vector  (da,db)**t .
c
c     the quantity r = (+/-)sqrt(da**2 + db**2) overwrites da in
c     storage.  the value of db is overwritten by a value z which
c     allows dc and ds to be recovered by the following algorithm.
c
c           if z=1  set  dc=0.0  and  ds=1.0
c           if abs(z) .lt. 1  set  dc=sqrt(1-z**2)  and  ds=z
c           if abs(z) .gt. 1  set  dc=1/z  and  ds=sqrt(1-dc**2)
c
c     normally, the subprogram drot(n,dx,incx,dy,incy,dc,ds) will
c     next be called to apply the transformation to a 2 by n matrix.
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
c   920501  reformatted the references section.  (wrb)
c***end prologue  drotg
      real*8  da, db, dc, ds, u, v, r
c***first executable statement  drotg
      if (abs(da) .le. abs(db)) go to 10
c
c *** here abs(da) .gt. abs(db) ***
c
      u = da + da
      v = db / u
c
c     note that u and r have the sign of da
c
      r = sqrt(0.25d0 + v**2) * u
c
c     note that dc is positive
c
      dc = da / r
      ds = v * (dc + dc)
      db = ds
      da = r
      return
c
c *** here abs(da) .le. abs(db) ***
c
   10 if (db .eq. 0.0d0) go to 20
      u = db + db
      v = da / u
c
c     note that u and r have the sign of db
c     (r is immediately stored in da)
c
      da = sqrt(0.25d0 + v**2) * u
c
c     note that ds is positive
c
      ds = db / da
      dc = v * (ds + ds)
      if (dc .eq. 0.0d0) go to 15
      db = 1.0d0 / dc
      return
   15 db = 1.0d0
      return
c
c *** here da = db = 0.0 ***
c
   20 dc = 1.0d0
      ds = 0.0d0
      return
c
      end
