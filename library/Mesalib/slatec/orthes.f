*deck orthes
      subroutine orthes (nm, n, low, igh, a, ort)
c***begin prologue  orthes
c***purpose  reduce a real general matrix to upper hessenberg form
c            using orthogonal similarity transformations.
c***library   slatec (eispack)
c***category  d4c1b2
c***type      single precision (orthes-s, corth-c)
c***keywords  eigenvalues, eigenvectors, eispack
c***author  smith, b. t., et al.
c***description
c
c     this subroutine is a translation of the algol procedure orthes,
c     num. math. 12, 349-368(1968) by martin and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
c
c     given a real general matrix, this subroutine
c     reduces a submatrix situated in rows and columns
c     low through igh to upper hessenberg form by
c     orthogonal similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of the two-dimensional
c          array parameter, a, as declared in the calling program
c          dimension statement.  nm is an integer variable.
c
c        n is the order of the matrix a.  n is an integer variable.
c          n must be less than or equal to nm.
c
c        low and igh are two integer variables determined by the
c          balancing subroutine  balanc.  if  balanc  has not been
c          used, set low=1 and igh equal to the order of the matrix, n.
c
c        a contains the general matrix to be reduced to upper
c          hessenberg form.  a is a two-dimensional real array,
c          dimensioned a(nm,n).
c
c     on output
c
c        a contains the upper hessenberg matrix.  some information about
c          the orthogonal transformations used in the reduction
c          is stored in the remaining triangle under the hessenberg
c          matrix.
c
c        ort contains further information about the orthogonal trans-
c          formations used in the reduction.  only elements low+1
c          through igh are used.  ort is a one-dimensional real array,
c          dimensioned ort(igh).
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
c***end prologue  orthes
c
      integer i,j,m,n,ii,jj,la,mp,nm,igh,kp1,low
      real a(nm,*),ort(*)
      real f,g,h,scale
c
c***first executable statement  orthes
      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200
c
      do 180 m = kp1, la
         h = 0.0e0
         ort(m) = 0.0e0
         scale = 0.0e0
c     .......... scale column (algol tol then not needed) ..........
         do 90 i = m, igh
   90    scale = scale + abs(a(i,m-1))
c
         if (scale .eq. 0.0e0) go to 180
         mp = m + igh
c     .......... for i=igh step -1 until m do -- ..........
         do 100 ii = m, igh
            i = mp - ii
            ort(i) = a(i,m-1) / scale
            h = h + ort(i) * ort(i)
  100    continue
c
         g = -sign(sqrt(h),ort(m))
         h = h - ort(m) * g
         ort(m) = ort(m) - g
c     .......... form (i-(u*ut)/h) * a ..........
         do 130 j = m, n
            f = 0.0e0
c     .......... for i=igh step -1 until m do -- ..........
            do 110 ii = m, igh
               i = mp - ii
               f = f + ort(i) * a(i,j)
  110       continue
c
            f = f / h
c
            do 120 i = m, igh
  120       a(i,j) = a(i,j) - f * ort(i)
c
  130    continue
c     .......... form (i-(u*ut)/h)*a*(i-(u*ut)/h) ..........
         do 160 i = 1, igh
            f = 0.0e0
c     .......... for j=igh step -1 until m do -- ..........
            do 140 jj = m, igh
               j = mp - jj
               f = f + ort(j) * a(i,j)
  140       continue
c
            f = f / h
c
            do 150 j = m, igh
  150       a(i,j) = a(i,j) - f * ort(j)
c
  160    continue
c
         ort(m) = scale * ort(m)
         a(m,m-1) = scale * g
  180 continue
c
  200 return
      end
