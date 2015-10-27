*deck ortbak
      subroutine ortbak (nm, low, igh, a, ort, m, z)
c***begin prologue  ortbak
c***purpose  form the eigenvectors of a general real matrix from the
c            eigenvectors of the upper hessenberg matrix output from
c            orthes.
c***library   slatec (eispack)
c***category  d4c4
c***type      single precision (ortbak-s, cortb-c)
c***keywords  eigenvalues, eigenvectors, eispack
c***author  smith, b. t., et al.
c***description
c
c     this subroutine is a translation of the algol procedure ortbak,
c     num. math. 12, 349-368(1968) by martin and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
c
c     this subroutine forms the eigenvectors of a real general
c     matrix by back transforming those of the corresponding
c     upper hessenberg matrix determined by  orthes.
c
c     on input
c
c        nm must be set to the row dimension of the two-dimensional
c          array parameters, a and z, as declared in the calling
c          program dimension statement.  nm is an integer variable.
c
c        low and igh are two integer variables determined by the
c          balancing subroutine  balanc.  if  balanc  has not been
c          used, set low=1 and igh equal to the order of the matrix.
c
c        a contains some information about the orthogonal trans-
c          formations used in the reduction to hessenberg form by
c          orthes  in its strict lower triangle.  a is a two-dimensional
c          real array, dimensioned a(nm,igh).
c
c        ort contains further information about the orthogonal trans-
c          formations used in the reduction by  orthes.  only elements
c          low through igh are used.  ort is a one-dimensional real
c          array, dimensioned ort(igh).
c
c        m is the number of columns of z to be back transformed.
c          m is an integer variable.
c
c        z contains the real and imaginary parts of the eigenvectors to
c          be back transformed in its first m columns.  z is a two-
c          dimensional real array, dimensioned z(nm,m).
c
c     on output
c
c        z contains the real and imaginary parts of the transformed
c          eigenvectors in its first m columns.
c
c        ort has been used for temporary storage as is not restored.
c
c     note that ortbak preserves vector euclidean norms.
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
c***end prologue  ortbak
c
      integer i,j,m,la,mm,mp,nm,igh,kp1,low,mp1
      real a(nm,*),ort(*),z(nm,*)
      real g
c
c***first executable statement  ortbak
      if (m .eq. 0) go to 200
      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200
c     .......... for mp=igh-1 step -1 until low+1 do -- ..........
      do 140 mm = kp1, la
         mp = low + igh - mm
         if (a(mp,mp-1) .eq. 0.0e0) go to 140
         mp1 = mp + 1
c
         do 100 i = mp1, igh
  100    ort(i) = a(i,mp-1)
c
         do 130 j = 1, m
            g = 0.0e0
c
            do 110 i = mp, igh
  110       g = g + ort(i) * z(i,j)
c     .......... divisor below is negative of h formed in orthes.
c                double division avoids possible underflow ..........
            g = (g / ort(mp)) / a(mp,mp-1)
c
            do 120 i = mp, igh
  120       z(i,j) = z(i,j) + g * ort(i)
c
  130    continue
c
  140 continue
c
  200 return
      end
