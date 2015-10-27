*deck elmhes
      subroutine elmhes (nm, n, low, igh, a, int)
c***begin prologue  elmhes
c***purpose  reduce a real general matrix to upper hessenberg form
c            using stabilized elementary similarity transformations.
c***library   slatec (eispack)
c***category  d4c1b2
c***type      single precision (elmhes-s, comhes-c)
c***keywords  eigenvalues, eigenvectors, eispack
c***author  smith, b. t., et al.
c***description
c
c     this subroutine is a translation of the algol procedure elmhes,
c     num. math. 12, 349-368(1968) by martin and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
c
c     given a real general matrix, this subroutine
c     reduces a submatrix situated in rows and columns
c     low through igh to upper hessenberg form by
c     stabilized elementary similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of the two-dimensional
c          array parameter, a, as declared in the calling program
c          dimension statement.  nm is an integer variable.
c
c        n is the order of the matrix, a.  n is an integer variable.
c          n must be less than or equal to nm.
c
c        low and igh are two integer variables determined by the
c          balancing subroutine  balanc.  if  balanc  has not been
c          used, set low=1 and igh equal to the order of the matrix, n.
c
c        a contains the input matrix.  a is a two-dimensional real
c          array, dimensioned a(nm,n).
c
c     on output
c
c        a contains the upper hessenberg matrix.  the multipliers which
c          were used in the reduction are stored in the remaining
c          triangle under the hessenberg matrix.
c
c        int contains information on the rows and columns interchanged
c          in the reduction.  only elements low through igh are used.
c          int is a one-dimensional integer array, dimensioned int(igh).
c
c     questions and comments should be directed to b. s. garbow,
c     applied mathematics division, argonne national laboratory
c     ------------------------------------------------------------------
c
c***references  b. t. smith, j. m. boyle, j. j. dongarra, b. s. garbow,
c                 y. ikebe, v. c. klema and c. b. moler, matrix eigen-
c                 system routines - eispack guide, springer-verlag,
c                 1976.
c***routines called  (none)
c***revision history  (yymmdd)
c   760101  date written
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  elmhes
c
      integer i,j,m,n,la,nm,igh,kp1,low,mm1,mp1
      real a(nm,*)
      real x,y
      integer int(*)
c
c***first executable statement  elmhes
      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200
c
      do 180 m = kp1, la
         mm1 = m - 1
         x = 0.0e0
         i = m
c
         do 100 j = m, igh
            if (abs(a(j,mm1)) .le. abs(x)) go to 100
            x = a(j,mm1)
            i = j
  100    continue
c
         int(m) = i
         if (i .eq. m) go to 130
c    .......... interchange rows and columns of a ..........
         do 110 j = mm1, n
            y = a(i,j)
            a(i,j) = a(m,j)
            a(m,j) = y
  110    continue
c
         do 120 j = 1, igh
            y = a(j,i)
            a(j,i) = a(j,m)
            a(j,m) = y
  120    continue
c    .......... end interchange ..........
  130    if (x .eq. 0.0e0) go to 180
         mp1 = m + 1
c
         do 160 i = mp1, igh
            y = a(i,mm1)
            if (y .eq. 0.0e0) go to 160
            y = y / x
            a(i,mm1) = y
c
            do 140 j = m, n
  140       a(i,j) = a(i,j) - y * a(m,j)
c
            do 150 j = 1, igh
  150       a(j,m) = a(j,m) + y * a(j,i)
c
  160    continue
c
  180 continue
c
  200 return
      end
