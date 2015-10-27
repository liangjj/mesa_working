*deck cbabk2
      subroutine cbabk2 (nm, n, low, igh, scale, m, zr, zi)
c***begin prologue  cbabk2
c***purpose  form the eigenvectors of a complex general matrix from the
c            eigenvectors of matrix output from cbal.
c***library   slatec (eispack)
c***category  d4c4
c***type      complex (balbak-s, cbabk2-c)
c***keywords  eigenvectors, eispack
c***author  smith, b. t., et al.
c***description
c
c     this subroutine is a translation of the algol procedure
c     cbabk2, which is a complex version of balbak,
c     num. math. 13, 293-304(1969) by parlett and reinsch.
c     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
c
c     this subroutine forms the eigenvectors of a complex general
c     matrix by back transforming those of the corresponding
c     balanced matrix determined by  cbal.
c
c     on input
c
c        nm must be set to the row dimension of the two-dimensional
c          array parameters, zr and zi, as declared in the calling
c          program dimension statement.  nm is an integer variable.
c
c        n is the order of the matrix z=(zr,zi).  n is an integer
c          variable.  n must be less than or equal to nm.
c
c        low and igh are integer variables determined by  cbal.
c
c        scale contains information determining the permutations and
c          scaling factors used by  cbal.  scale is a one-dimensional
c          real array, dimensioned scale(n).
c
c        m is the number of eigenvectors to be back transformed.
c          m is an integer variable.
c
c        zr and zi contain the real and imaginary parts, respectively,
c          of the eigenvectors to be back transformed in their first
c          m columns.  zr and zi are two-dimensional real arrays,
c          dimensioned zr(nm,m) and zi(nm,m).
c
c     on output
c
c        zr and zi contain the real and imaginary parts,
c          respectively, of the transformed eigenvectors
c          in their first m columns.
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
c***end prologue  cbabk2
c
      integer i,j,k,m,n,ii,nm,igh,low
      real scale(*),zr(nm,*),zi(nm,*)
      real s
c
c***first executable statement  cbabk2
      if (m .eq. 0) go to 200
      if (igh .eq. low) go to 120
c
      do 110 i = low, igh
         s = scale(i)
c     .......... left hand eigenvectors are back transformed
c                if the foregoing statement is replaced by
c                s=1.0e0/scale(i). ..........
         do 100 j = 1, m
            zr(i,j) = zr(i,j) * s
            zi(i,j) = zi(i,j) * s
  100    continue
c
  110 continue
c     .......... for i=low-1 step -1 until 1,
c                igh+1 step 1 until n do -- ..........
  120 do 140 ii = 1, n
         i = ii
         if (i .ge. low .and. i .le. igh) go to 140
         if (i .lt. low) i = low - ii
         k = scale(i)
         if (k .eq. i) go to 140
c
         do 130 j = 1, m
            s = zr(i,j)
            zr(i,j) = zr(k,j)
            zr(k,j) = s
            s = zi(i,j)
            zi(i,j) = zi(k,j)
            zi(k,j) = s
  130    continue
c
  140 continue
c
  200 return
      end
