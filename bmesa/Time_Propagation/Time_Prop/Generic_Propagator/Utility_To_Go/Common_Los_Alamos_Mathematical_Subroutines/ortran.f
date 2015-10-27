*deck ortran
      subroutine ortran (nm, n, low, igh, a, ort, z)
c***begin prologue  ortran
c***purpose  accumulate orthogonal similarity transformations in the
c            reduction of real general matrix by orthes.
c***library   slatec (eispack)
c***category  d4c4
c***type      single precision (ortran-s)
c***keywords  eigenvalues, eigenvectors, eispack
c***author  smith, b. t., et al.
c***description
c
c     this subroutine is a translation of the algol procedure ortrans,
c     num. math. 16, 181-204(1970) by peters and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
c
c     this subroutine accumulates the orthogonal similarity
c     transformations used in the reduction of a real general
c     matrix to upper hessenberg form by  orthes.
c
c     on input
c
c        nm must be set to the row dimension of the two-dimensional
c          array parameters, a and z, as declared in the calling
c          program dimension statement.  nm is an integer variable.
c
c        n is the order of the matrix a.  n is an integer variable.
c          n must be less than or equal to nm.
c
c        low and igh are two integer variables determined by the
c          balancing subroutine  balanc.  if  balanc  has not been
c          used, set low=1 and igh equal to the order of the matrix, n.
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
c     on output
c
c        z contains the transformation matrix produced in the reduction
c          by  orthes  to the upper hessenberg form.  z is a two-
c          dimensional real array, dimensioned z(nm,n).
c
c        ort has been used for temporary storage as is not restored.
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
c***end prologue  ortran
c
      integer i,j,n,kl,mm,mp,nm,igh,low,mp1
      real*8 a(nm,*),ort(*),z(nm,*)
      real*8 g
c
c     .......... initialize z to identity matrix ..........
c***first executable statement  ortran
      do 80 i = 1, n
c
         do 60 j = 1, n
   60    z(i,j) = 0.0d0
c
         z(i,i) = 1.0d0
   80 continue
c
      kl = igh - low - 1
      if (kl .lt. 1) go to 200
c     .......... for mp=igh-1 step -1 until low+1 do -- ..........
      do 140 mm = 1, kl
         mp = igh - mm
         if (a(mp,mp-1) .eq. 0.0d0) go to 140
         mp1 = mp + 1
c
         do 100 i = mp1, igh
  100    ort(i) = a(i,mp-1)
c
         do 130 j = mp, igh
            g = 0.0d0
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
