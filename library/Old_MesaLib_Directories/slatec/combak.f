*deck combak
      subroutine combak (nm, low, igh, ar, ai, int, m, zr, zi)
c***begin prologue  combak
c***purpose  form the eigenvectors of a complex general matrix from the
c            eigenvectors of a upper hessenberg matrix output from
c            comhes.
c***library   slatec (eispack)
c***category  d4c4
c***type      complex (elmbak-s, combak-c)
c***keywords  eigenvalues, eigenvectors, eispack
c***author  smith, b. t., et al.
c***description
c
c     this subroutine is a translation of the algol procedure combak,
c     num. math. 12, 349-368(1968) by martin and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
c
c     this subroutine forms the eigenvectors of a complex general
c     matrix by back transforming those of the corresponding
c     upper hessenberg matrix determined by  comhes.
c
c     on input
c
c        nm must be set to the row dimension of the two-dimensional
c          array parameters, ar, ai, zr and zi, as declared in the
c          calling program dimension statement.  nm is an integer
c          variable.
c
c        low and igh are two integer variables determined by the
c          balancing subroutine  cbal.  if  cbal  has not been used,
c          set low=1 and igh equal to the order of the matrix.
c
c        ar and ai contain the multipliers which were used in the
c           reduction by  comhes  in their lower triangles below
c           the subdiagonal.  ar and ai are two-dimensional real
c           arrays, dimensioned ar(nm,igh) and ai(nm,igh).
c
c        int contains information on the rows and columns
c          interchanged in the reduction by  comhes.  only
c          elements low through igh are used.  int is a
c          one-dimensional integer array, dimensioned int(igh).
c
c        m is the number of eigenvectors to be back transformed.
c          m is an integer variable.
c
c        zr and zi contain the real and imaginary parts, respectively,
c          of the eigenvectors to be back transformed in their first m
c          columns.  zr and zi are two-dimensional real arrays,
c          dimensioned zr(nm,m) and zi(nm,m).
c
c     on output
c
c        zr and zi contain the real and imaginary parts, respectively,
c          of the transformed eigenvectors in their first m columns.
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
c***end prologue  combak
c
      integer i,j,m,la,mm,mp,nm,igh,kp1,low,mp1
      real ar(nm,*),ai(nm,*),zr(nm,*),zi(nm,*)
      real xr,xi
      integer int(*)
c
c***first executable statement  combak
      if (m .eq. 0) go to 200
      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200
c     .......... for mp=igh-1 step -1 until low+1 do -- ..........
      do 140 mm = kp1, la
         mp = low + igh - mm
         mp1 = mp + 1
c
         do 110 i = mp1, igh
            xr = ar(i,mp-1)
            xi = ai(i,mp-1)
            if (xr .eq. 0.0e0 .and. xi .eq. 0.0e0) go to 110
c
            do 100 j = 1, m
               zr(i,j) = zr(i,j) + xr * zr(mp,j) - xi * zi(mp,j)
               zi(i,j) = zi(i,j) + xr * zi(mp,j) + xi * zr(mp,j)
  100       continue
c
  110    continue
c
         i = int(mp)
         if (i .eq. mp) go to 140
c
         do 130 j = 1, m
            xr = zr(i,j)
            zr(i,j) = zr(mp,j)
            zr(mp,j) = xr
            xi = zi(i,j)
            zi(i,j) = zi(mp,j)
            zi(mp,j) = xi
  130    continue
c
  140 continue
c
  200 return
      end
