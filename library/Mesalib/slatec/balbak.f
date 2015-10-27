*deck balbak
      subroutine balbak (nm, n, low, igh, scale, m, z)
c***begin prologue  balbak
c***purpose  form the eigenvectors of a real general matrix from the
c            eigenvectors of matrix output from balanc.
c***library   slatec (eispack)
c***category  d4c4
c***type      single precision (balbak-s, cbabk2-c)
c***keywords  eigenvectors, eispack
c***author  smith, b. t., et al.
c***description
c
c     this subroutine is a translation of the algol procedure balbak,
c     num. math. 13, 293-304(1969) by parlett and reinsch.
c     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
c
c     this subroutine forms the eigenvectors of a real general
c     matrix by back transforming those of the corresponding
c     balanced matrix determined by  balanc.
c
c     on input
c
c        nm must be set to the row dimension of the two-dimensional
c          array parameter, z, as declared in the calling program
c          dimension statement.  nm is an integer variable.
c
c        n is the number of components of the vectors in matrix z.
c          n is an integer variable.  n must be less than or equal
c          to nm.
c
c        low and igh are integer variables determined by  balanc.
c
c        scale contains information determining the permutations and
c          scaling factors used by  balanc.  scale is a one-dimensional
c          real array, dimensioned scale(n).
c
c        m is the number of columns of z to be back transformed.
c          m is an integer variable.
c
c        z contains the real and imaginary parts of the eigen-
c          vectors to be back transformed in its first m columns.
c          z is a two-dimensional real array, dimensioned z(nm,m).
c
c     on output
c
c        z contains the real and imaginary parts of the
c          transformed eigenvectors in its first m columns.
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
c***end prologue  balbak
c
      integer i,j,k,m,n,ii,nm,igh,low
      real scale(*),z(nm,*)
      real s
c
c***first executable statement  balbak
      if (m .eq. 0) go to 200
      if (igh .eq. low) go to 120
c
      do 110 i = low, igh
         s = scale(i)
c     .......... left hand eigenvectors are back transformed
c                if the foregoing statement is replaced by
c                s=1.0e0/scale(i). ..........
         do 100 j = 1, m
  100    z(i,j) = z(i,j) * s
c
  110 continue
c     ......... for i=low-1 step -1 until 1,
c               igh+1 step 1 until n do -- ..........
  120 do 140 ii = 1, n
         i = ii
         if (i .ge. low .and. i .le. igh) go to 140
         if (i .lt. low) i = low - ii
         k = scale(i)
         if (k .eq. i) go to 140
c
         do 130 j = 1, m
            s = z(i,j)
            z(i,j) = z(k,j)
            z(k,j) = s
  130    continue
c
  140 continue
c
  200 return
      end
