*deck corth
      subroutine corth (nm, n, low, igh, ar, ai, ortr, orti)
c***begin prologue  corth
c***purpose  reduce a complex general matrix to complex upper hessenberg
c            form using unitary similarity transformations.
c***library   slatec (eispack)
c***category  d4c1b2
c***type      complex (orthes-s, corth-c)
c***keywords  eigenvalues, eigenvectors, eispack
c***author  smith, b. t., et al.
c***description
c
c     this subroutine is a translation of a complex analogue of
c     the algol procedure orthes, num. math. 12, 349-368(1968)
c     by martin and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
c
c     given a complex general matrix, this subroutine
c     reduces a submatrix situated in rows and columns
c     low through igh to upper hessenberg form by
c     unitary similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of the two-dimensional
c          array parameters, ar and ai, as declared in the calling
c          program dimension statement.  nm is an integer variable.
c
c        n is the order of the matrix a=(ar,ai).  n is an integer
c          variable.  n must be less than or equal to nm.
c
c        low and igh are two integer variables determined by the
c          balancing subroutine  cbal.  if  cbal  has not been used,
c          set low=1 and igh equal to the order of the matrix, n.
c
c        ar and ai contain the real and imaginary parts, respectively,
c          of the complex input matrix.  ar and ai are two-dimensional
c          real arrays, dimensioned ar(nm,n) and ai(nm,n).
c
c     on output
c
c        ar and ai contain the real and imaginary parts, respectively,
c          of the hessenberg matrix.  information about the unitary
c          transformations used in the reduction is stored in the
c          remaining triangles under the hessenberg matrix.
c
c        ortr and orti contain further information about the unitary
c          transformations.  only elements low through igh are used.
c          ortr and orti are one-dimensional real arrays, dimensioned
c          ortr(igh) and orti(igh).
c
c     calls pythag(a,b) for sqrt(a**2 + b**2).
c
c     questions and comments should be directed to b. s. garbow,
c     applied mathematics division, argonne national laboratory
c     ------------------------------------------------------------------
c
c***references  b. t. smith, j. m. boyle, j. j. dongarra, b. s. garbow,
c                 y. ikebe, v. c. klema and c. b. moler, matrix eigen-
c                 system routines - eispack guide, springer-verlag,
c                 1976.
c***routines called  pythag
c***revision history  (yymmdd)
c   760101  date written
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  corth
c
      integer i,j,m,n,ii,jj,la,mp,nm,igh,kp1,low
      real ar(nm,*),ai(nm,*),ortr(*),orti(*)
      real f,g,h,fi,fr,scale
      real pythag
c
c***first executable statement  corth
      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200
c
      do 180 m = kp1, la
         h = 0.0e0
         ortr(m) = 0.0e0
         orti(m) = 0.0e0
         scale = 0.0e0
c     .......... scale column (algol tol then not needed) ..........
         do 90 i = m, igh
   90    scale = scale + abs(ar(i,m-1)) + abs(ai(i,m-1))
c
         if (scale .eq. 0.0e0) go to 180
         mp = m + igh
c     .......... for i=igh step -1 until m do -- ..........
         do 100 ii = m, igh
            i = mp - ii
            ortr(i) = ar(i,m-1) / scale
            orti(i) = ai(i,m-1) / scale
            h = h + ortr(i) * ortr(i) + orti(i) * orti(i)
  100    continue
c
         g = sqrt(h)
         f = pythag(ortr(m),orti(m))
         if (f .eq. 0.0e0) go to 103
         h = h + f * g
         g = g / f
         ortr(m) = (1.0e0 + g) * ortr(m)
         orti(m) = (1.0e0 + g) * orti(m)
         go to 105
c
  103    ortr(m) = g
         ar(m,m-1) = scale
c     .......... form (i-(u*ut)/h) * a ..........
  105    do 130 j = m, n
            fr = 0.0e0
            fi = 0.0e0
c     .......... for i=igh step -1 until m do -- ..........
            do 110 ii = m, igh
               i = mp - ii
               fr = fr + ortr(i) * ar(i,j) + orti(i) * ai(i,j)
               fi = fi + ortr(i) * ai(i,j) - orti(i) * ar(i,j)
  110       continue
c
            fr = fr / h
            fi = fi / h
c
            do 120 i = m, igh
               ar(i,j) = ar(i,j) - fr * ortr(i) + fi * orti(i)
               ai(i,j) = ai(i,j) - fr * orti(i) - fi * ortr(i)
  120       continue
c
  130    continue
c     .......... form (i-(u*ut)/h)*a*(i-(u*ut)/h) ..........
         do 160 i = 1, igh
            fr = 0.0e0
            fi = 0.0e0
c     .......... for j=igh step -1 until m do -- ..........
            do 140 jj = m, igh
               j = mp - jj
               fr = fr + ortr(j) * ar(i,j) - orti(j) * ai(i,j)
               fi = fi + ortr(j) * ai(i,j) + orti(j) * ar(i,j)
  140       continue
c
            fr = fr / h
            fi = fi / h
c
            do 150 j = m, igh
               ar(i,j) = ar(i,j) - fr * ortr(j) - fi * orti(j)
               ai(i,j) = ai(i,j) + fr * orti(j) - fi * ortr(j)
  150       continue
c
  160    continue
c
         ortr(m) = scale * ortr(m)
         orti(m) = scale * orti(m)
         ar(m,m-1) = -g * ar(m,m-1)
         ai(m,m-1) = -g * ai(m,m-1)
  180 continue
c
  200 return
      end
