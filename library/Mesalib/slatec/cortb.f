*deck cortb
      subroutine cortb (nm, low, igh, ar, ai, ortr, orti, m, zr, zi)
c***begin prologue  cortb
c***purpose  form the eigenvectors of a complex general matrix from
c            eigenvectors of upper hessenberg matrix output from
c            corth.
c***library   slatec (eispack)
c***category  d4c4
c***type      complex (ortbak-s, cortb-c)
c***keywords  eigenvalues, eigenvectors, eispack
c***author  smith, b. t., et al.
c***description
c
c     this subroutine is a translation of a complex analogue of
c     the algol procedure ortbak, num. math. 12, 349-368(1968)
c     by martin and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
c
c     this subroutine forms the eigenvectors of a complex general
c     matrix by back transforming those of the corresponding
c     upper hessenberg matrix determined by  corth.
c
c     on input
c
c        nm must be set to the row dimension of the two-dimensional
c          array parameters, ar, ai, zr, and zi, as declared in the
c          calling program dimension statement.  nm is an integer
c          variable.
c
c        low and igh are two integer variables determined by the
c          balancing subroutine  cbal.  if  cbal  has not been used,
c          set low=1 and igh equal to the order of the matrix.
c
c        ar and ai contain information about the unitary trans-
c          formations used in the reduction by  corth  in their
c          strict lower triangles.  ar and ai are two-dimensional
c          real arrays, dimensioned ar(nm,igh) and ai(nm,igh).
c
c        ortr and orti contain further information about the unitary
c          transformations used in the reduction by  corth.  only
c          elements low through igh are used.  ortr and orti are
c          one-dimensional real arrays, dimensioned ortr(igh) and
c          orti(igh).
c
c        m is the number of columns of z=(zr,zi) to be back transformed.
c          m is an integer variable.
c
c        zr and zi contain the real and imaginary parts, respectively,
c          of the eigenvectors to be back transformed in their first
c          m columns.  zr and zi are two-dimensional real arrays,
c          dimensioned zr(nm,m) and zi(nm,m).
c
c     on output
c
c        zr and zi contain the real and imaginary parts, respectively,
c          of the transformed eigenvectors in their first m columns.
c
c        ortr and orti have been altered.
c
c     note that cortb preserves vector euclidean norms.
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
c***end prologue  cortb
c
      integer i,j,m,la,mm,mp,nm,igh,kp1,low,mp1
      real ar(nm,*),ai(nm,*),ortr(*),orti(*)
      real zr(nm,*),zi(nm,*)
      real h,gi,gr
c
c***first executable statement  cortb
      if (m .eq. 0) go to 200
      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200
c     .......... for mp=igh-1 step -1 until low+1 do -- ..........
      do 140 mm = kp1, la
         mp = low + igh - mm
         if (ar(mp,mp-1) .eq. 0.0e0 .and. ai(mp,mp-1) .eq. 0.0e0)
     1      go to 140
c     .......... h below is negative of h formed in corth ..........
         h = ar(mp,mp-1) * ortr(mp) + ai(mp,mp-1) * orti(mp)
         mp1 = mp + 1
c
         do 100 i = mp1, igh
            ortr(i) = ar(i,mp-1)
            orti(i) = ai(i,mp-1)
  100    continue
c
         do 130 j = 1, m
            gr = 0.0e0
            gi = 0.0e0
c
            do 110 i = mp, igh
               gr = gr + ortr(i) * zr(i,j) + orti(i) * zi(i,j)
               gi = gi + ortr(i) * zi(i,j) - orti(i) * zr(i,j)
  110       continue
c
            gr = gr / h
            gi = gi / h
c
            do 120 i = mp, igh
               zr(i,j) = zr(i,j) + gr * ortr(i) - gi * orti(i)
               zi(i,j) = zi(i,j) + gr * orti(i) + gi * ortr(i)
  120       continue
c
  130    continue
c
  140 continue
c
  200 return
      end
