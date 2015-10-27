*deck eltran
      subroutine eltran (nm, n, low, igh, a, int, z)
c***begin prologue  eltran
c***purpose  accumulates the stabilized elementary similarity
c            transformations used in the reduction of a real general
c            matrix to upper hessenberg form by elmhes.
c***library   slatec (eispack)
c***category  d4c4
c***type      single precision (eltran-s)
c***keywords  eigenvalues, eigenvectors, eispack
c***author  smith, b. t., et al.
c***description
c
c     this subroutine is a translation of the algol procedure elmtrans,
c     num. math. 16, 181-204(1970) by peters and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
c
c     this subroutine accumulates the stabilized elementary
c     similarity transformations used in the reduction of a
c     real general matrix to upper hessenberg form by  elmhes.
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
c        a contains the multipliers which were used in the reduction
c          by  elmhes  in its lower triangle below the subdiagonal.
c          a is a two-dimensional real array, dimensioned a(nm,igh).
c
c        int contains information on the rows and columns interchanged
c          in the reduction by  elmhes.  only elements low through igh
c          are used.  int is a one-dimensional integer array,
c          dimensioned int(igh).
c
c     on output
c
c        z contains the transformation matrix produced in the reduction
c          by  elmhes.  z is a two-dimensional real array, dimensioned
c          z(nm,n).
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
c***end prologue  eltran
c
      integer i,j,n,kl,mm,mp,nm,igh,low,mp1
      real a(nm,*),z(nm,*)
      integer int(*)
c
c***first executable statement  eltran
      do 80 i = 1, n
c
         do 60 j = 1, n
   60    z(i,j) = 0.0e0
c
         z(i,i) = 1.0e0
   80 continue
c
      kl = igh - low - 1
      if (kl .lt. 1) go to 200
c     .......... for mp=igh-1 step -1 until low+1 do -- ..........
      do 140 mm = 1, kl
         mp = igh - mm
         mp1 = mp + 1
c
         do 100 i = mp1, igh
  100    z(i,mp) = a(i,mp-1)
c
         i = int(mp)
         if (i .eq. mp) go to 140
c
         do 130 j = mp, igh
            z(mp,j) = z(i,j)
            z(i,j) = 0.0e0
  130    continue
c
         z(i,mp) = 1.0e0
  140 continue
c
  200 return
      end
