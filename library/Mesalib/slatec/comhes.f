*deck comhes
      subroutine comhes (nm, n, low, igh, ar, ai, int)
c***begin prologue  comhes
c***purpose  reduce a complex general matrix to complex upper hessenberg
c            form using stabilized elementary similarity
c            transformations.
c***library   slatec (eispack)
c***category  d4c1b2
c***type      complex (elmhes-s, comhes-c)
c***keywords  eigenvalues, eigenvectors, eispack
c***author  smith, b. t., et al.
c***description
c
c     this subroutine is a translation of the algol procedure comhes,
c     num. math. 12, 349-368(1968) by martin and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
c
c     given a complex general matrix, this subroutine
c     reduces a submatrix situated in rows and columns
c     low through igh to upper hessenberg form by
c     stabilized elementary similarity transformations.
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
c          of the upper hessenberg matrix.  the multipliers which
c          were used in the reduction are stored in the remaining
c          triangles under the hessenberg matrix.
c
c        int contains information on the rows and columns
c          interchanged in the reduction.  only elements low through
c          igh are used.  int is a one-dimensional integer array,
c          dimensioned int(igh).
c
c     calls cdiv for complex division.
c
c     questions and comments should be directed to b. s. garbow,
c     applied mathematics division, argonne national laboratory
c     ------------------------------------------------------------------
c
c***references  b. t. smith, j. m. boyle, j. j. dongarra, b. s. garbow,
c                 y. ikebe, v. c. klema and c. b. moler, matrix eigen-
c                 system routines - eispack guide, springer-verlag,
c                 1976.
c***routines called  cdiv
c***revision history  (yymmdd)
c   760101  date written
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  comhes
c
      integer i,j,m,n,la,nm,igh,kp1,low,mm1,mp1
      real ar(nm,*),ai(nm,*)
      real xr,xi,yr,yi
      integer int(*)
c
c***first executable statement  comhes
      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200
c
      do 180 m = kp1, la
         mm1 = m - 1
         xr = 0.0e0
         xi = 0.0e0
         i = m
c
         do 100 j = m, igh
            if (abs(ar(j,mm1)) + abs(ai(j,mm1))
     1         .le. abs(xr) + abs(xi)) go to 100
            xr = ar(j,mm1)
            xi = ai(j,mm1)
            i = j
  100    continue
c
         int(m) = i
         if (i .eq. m) go to 130
c     .......... interchange rows and columns of ar and ai ..........
         do 110 j = mm1, n
            yr = ar(i,j)
            ar(i,j) = ar(m,j)
            ar(m,j) = yr
            yi = ai(i,j)
            ai(i,j) = ai(m,j)
            ai(m,j) = yi
  110    continue
c
         do 120 j = 1, igh
            yr = ar(j,i)
            ar(j,i) = ar(j,m)
            ar(j,m) = yr
            yi = ai(j,i)
            ai(j,i) = ai(j,m)
            ai(j,m) = yi
  120    continue
c     .......... end interchange ..........
  130    if (xr .eq. 0.0e0 .and. xi .eq. 0.0e0) go to 180
         mp1 = m + 1
c
         do 160 i = mp1, igh
            yr = ar(i,mm1)
            yi = ai(i,mm1)
            if (yr .eq. 0.0e0 .and. yi .eq. 0.0e0) go to 160
            call cdiv(yr,yi,xr,xi,yr,yi)
            ar(i,mm1) = yr
            ai(i,mm1) = yi
c
            do 140 j = m, n
               ar(i,j) = ar(i,j) - yr * ar(m,j) + yi * ai(m,j)
               ai(i,j) = ai(i,j) - yr * ai(m,j) - yi * ar(m,j)
  140       continue
c
            do 150 j = 1, igh
               ar(j,m) = ar(j,m) + yr * ar(j,i) - yi * ai(j,i)
               ai(j,m) = ai(j,m) + yr * ai(j,i) + yi * ar(j,i)
  150       continue
c
  160    continue
c
  180 continue
c
  200 return
      end
