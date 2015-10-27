*deck rebakb
      subroutine rebakb (nm, n, b, dl, m, z)
c***begin prologue  rebakb
c***purpose  form the eigenvectors of a generalized symmetric
c            eigensystem from the eigenvectors of derived matrix output
c            from reduc2.
c***library   slatec (eispack)
c***category  d4c4
c***type      single precision (rebakb-s)
c***keywords  eigenvalues, eigenvectors, eispack
c***author  smith, b. t., et al.
c***description
c
c     this subroutine is a translation of the algol procedure rebakb,
c     num. math. 11, 99-110(1968) by martin and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 303-314(1971).
c
c     this subroutine forms the eigenvectors of a generalized
c     symmetric eigensystem by back transforming those of the
c     derived symmetric matrix determined by  reduc2.
c
c     on input
c
c        nm must be set to the row dimension of the two-dimensional
c          array parameters, b and z, as declared in the calling
c          program dimension statement.  nm is an integer variable.
c
c        n is the order of the matrix system.  n is an integer
c          variable.  n must be less than or equal to nm.
c
c        b contains information about the similarity transformation
c          (cholesky decomposition) used in the reduction by  reduc2
c          in its strict lower triangle.  b is a two-dimensional
c          real array, dimensioned b(nm,n).
c
c        dl contains further information about the transformation.
c          dl is a one-dimensional real array, dimensioned dl(n).
c
c        m is the number of eigenvectors to be back transformed.
c          m is an integer variable.
c
c        z contains the eigenvectors to be back transformed in its
c          first m columns.  z is a two-dimensional real array
c          dimensioned z(nm,m).
c
c     on output
c
c        z contains the transformed eigenvectors in its first
c          m columns.
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
c***end prologue  rebakb
c
      integer i,j,k,m,n,i1,ii,nm
      real b(nm,*),dl(*),z(nm,*)
      real x
c
c***first executable statement  rebakb
      if (m .eq. 0) go to 200
c
      do 100 j = 1, m
c     .......... for i=n step -1 until 1 do -- ..........
         do 100 ii = 1, n
            i1 = n - ii
            i = i1 + 1
            x = dl(i) * z(i,j)
            if (i .eq. 1) go to 80
c
            do 60 k = 1, i1
   60       x = x + b(i,k) * z(k,j)
c
   80       z(i,j) = x
  100 continue
c
  200 return
      end
