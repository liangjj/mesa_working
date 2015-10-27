*deck snbco
      subroutine snbco (abe, lda, n, ml, mu, ipvt, rcond, z)
c***begin prologue  snbco
c***purpose  factor a band matrix using gaussian elimination and
c            estimate the condition number.
c***library   slatec
c***category  d2a2
c***type      single precision (snbco-s, dnbco-d, cnbco-c)
c***keywords  banded, linear equations, matrix factorization,
c             nonsymmetric
c***author  voorhees, e. a., (lanl)
c***description
c
c     snbco factors a real band matrix by gaussian
c     elimination and estimates the condition of the matrix.
c
c     if rcond is not needed, snbfa is slightly faster.
c     to solve  a*x = b , follow snbco by snbsl.
c     to compute  inverse(a)*c , follow snbco by snbsl.
c     to compute  determinant(a) , follow snbco by snbdi.
c
c     on entry
c
c        abe     real(lda, nc)
c                contains the matrix in band storage.  the rows
c                of the original matrix are stored in the rows
c                of abe and the diagonals of the original matrix
c                are stored in columns 1 through ml+mu+1 of abe.
c                nc must be .ge. 2*ml+mu+1 .
c                see the comments below for details.
c
c        lda     integer
c                the leading dimension of the array abe.
c                lda must be .ge. n .
c
c        n       integer
c                the order of the original matrix.
c
c        ml      integer
c                number of diagonals below the main diagonal.
c                0 .le. ml .lt. n .
c
c        mu      integer
c                number of diagonals above the main diagonal.
c                0 .le. mu .lt. n .
c                more efficient if ml .le. mu .
c
c     on return
c
c        abe     an upper triangular matrix in band storage
c                and the multipliers which were used to obtain it.
c                the factorization can be written  a = l*u , where
c                l is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        rcond   real
c                an estimate of the reciprocal condition of  a .
c                for the system  a*x = b , relative perturbations
c                in  a  and  b  of size  epsilon  may cause
c                relative perturbations in  x  of size  epsilon/rcond .
c                if  rcond  is so small that the logical expression
c                         1.0 + rcond .eq. 1.0
c                is true, then  a  may be singular to working
c                precision.  in particular,  rcond  is zero  if
c                exact singularity is detected or the estimate
c                underflows.
c
c        z       real(n)
c                a work vector whose contents are usually unimportant.
c                if  a  is close to a singular matrix, then  z  is
c                an approximate null vector in the sense that
c                norm(a*z) = rcond*norm(a)*norm(z) .
c
c     band storage
c
c           if  a  is a band matrix, the following program segment
c           will set up the input.
c
c                   ml = (band width below the diagonal)
c                   mu = (band width above the diagonal)
c                   do 20 i = 1, n
c                      j1 = max(1, i-ml)
c                      j2 = min(n, i+mu)
c                      do 10 j = j1, j2
c                         k = j - i + ml + 1
c                         abe(i,k) = a(i,j)
c                10    continue
c                20 continue
c
c           this uses columns  1  through  ml+mu+1  of abe .
c           furthermore,  ml  additional columns are needed in
c           abe  starting with column  ml+mu+2  for elements
c           generated during the triangularization.  the total
c           number of columns needed in  abe  is  2*ml+mu+1 .
c
c     example:  if the original matrix is
c
c           11 12 13  0  0  0
c           21 22 23 24  0  0
c            0 32 33 34 35  0
c            0  0 43 44 45 46
c            0  0  0 54 55 56
c            0  0  0  0 65 66
c
c      then  n = 6, ml = 1, mu = 2, lda .ge. 5  and abe should contain
c
c            * 11 12 13  +     , * = not used
c           21 22 23 24  +     , + = used for pivoting
c           32 33 34 35  +
c           43 44 45 46  +
c           54 55 56  *  +
c           65 66  *  *  +
c
c***references  j. j. dongarra, j. r. bunch, c. b. moler, and g. w.
c                 stewart, linpack users' guide, siam, 1979.
c***routines called  sasum, saxpy, sdot, snbfa, sscal
c***revision history  (yymmdd)
c   800723  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  snbco
      integer lda,n,ml,mu,ipvt(*)
      real abe(lda,*),z(*)
      real rcond
c
      real sdot,ek,t,wk,wkm
      real anorm,s,sasum,sm,ynorm
      integer i,info,j,ju,k,kb,kp1,l,ldb,lm,lz,m,ml1,mm,nl,nu
c***first executable statement  snbco
      ml1=ml+1
      ldb = lda - 1
      anorm = 0.0e0
      do 10 j = 1, n
        nu = min(mu,j-1)
        nl = min(ml,n-j)
        l = 1 + nu + nl
        anorm = max(anorm,sasum(l,abe(j+nl,ml1-nl),ldb))
   10 continue
c
c     factor
c
      call snbfa(abe,lda,n,ml,mu,ipvt,info)
c
c     rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) .
c     estimate = norm(z)/norm(y) where  a*z = y  and trans(a)*y = e .
c     trans(a)  is the transpose of a .  the components of  e  are
c     chosen to cause maximum local growth in the elements of  w where
c     trans(u)*w = e .  the vectors are frequently rescaled to avoid
c     overflow.
c
c     solve trans(u)*w = e
c
      ek = 1.0e0
      do 20 j = 1, n
        z(j) = 0.0e0
   20 continue
      m = ml + mu + 1
      ju = 0
      do 100 k = 1, n
        if (z(k) .ne. 0.0e0) ek = sign(ek,-z(k))
        if (abs(ek-z(k)) .le. abs(abe(k,ml1))) go to 30
          s = abs(abe(k,ml1))/abs(ek-z(k))
          call sscal(n,s,z,1)
          ek = s*ek
   30   continue
        wk = ek - z(k)
        wkm = -ek - z(k)
        s = abs(wk)
        sm = abs(wkm)
        if (abe(k,ml1) .eq. 0.0e0) go to 40
          wk = wk/abe(k,ml1)
          wkm = wkm/abe(k,ml1)
        go to 50
   40   continue
          wk = 1.0e0
          wkm = 1.0e0
   50   continue
        kp1 = k + 1
        ju = min(max(ju,mu+ipvt(k)),n)
        mm = ml1
        if (kp1 .gt. ju) go to 90
          do 60 i = kp1, ju
            mm = mm + 1
            sm = sm + abs(z(i)+wkm*abe(k,mm))
            z(i) = z(i) + wk*abe(k,mm)
            s = s + abs(z(i))
   60     continue
          if (s .ge. sm) go to 80
            t = wkm -wk
            wk = wkm
            mm = ml1
            do 70 i = kp1, ju
              mm = mm + 1
              z(i) = z(i) + t*abe(k,mm)
   70       continue
   80     continue
   90   continue
      z(k) = wk
  100 continue
      s = 1.0e0/sasum(n,z,1)
      call sscal(n,s,z,1)
c
c     solve trans(l)*y = w
c
      do 120 kb = 1, n
        k = n + 1 - kb
        nl = min(ml,n-k)
        if (k .lt. n) z(k) = z(k) + sdot(nl,abe(k+nl,ml1-nl),-ldb,z(k+1)
     1  ,1)
        if (abs(z(k)) .le. 1.0e0) go to 110
          s = 1.0e0/abs(z(k))
          call sscal(n,s,z,1)
  110   continue
        l = ipvt(k)
        t = z(l)
        z(l) = z(k)
        z(k) = t
  120 continue
      s = 1.0e0/sasum(n,z,1)
      call sscal(n,s,z,1)
c
      ynorm = 1.0e0
c
c     solve l*v = y
c
      do 140 k = 1, n
        l = ipvt(k)
        t = z(l)
        z(l) = z(k)
        z(k) = t
        nl = min(ml,n-k)
        if (k .lt. n) call saxpy(nl,t,abe(k+nl,ml1-nl),-ldb,z(k+1),1)
        if (abs(z(k)) .le. 1.0e0) go to 130
          s = 1.0e0/abs(z(k))
          call sscal(n,s,z,1)
          ynorm = s*ynorm
  130   continue
  140 continue
      s = 1.0e0/sasum(n,z,1)
      call sscal(n,s,z,1)
      ynorm = s*ynorm
c
c     solve  u*z = v
c
      do 160 kb = 1, n
        k = n + 1 - kb
        if (abs(z(k)) .le. abs(abe(k,ml1))) go to 150
          s = abs(abe(k,ml1))/abs(z(k))
          call sscal(n,s,z,1)
          ynorm = s*ynorm
  150   continue
        if (abe(k,ml1) .ne. 0.0e0) z(k) = z(k)/abe(k,ml1)
        if (abe(k,ml1) .eq. 0.0e0) z(k) = 1.0e0
        lm = min(k,m) - 1
        lz = k - lm
        t = -z(k)
        call saxpy(lm,t,abe(k-1,ml+2),-ldb,z(lz),1)
  160 continue
c     make znorm = 1.0e0
      s = 1.0e0/sasum(n,z,1)
      call sscal(n,s,z,1)
      ynorm = s*ynorm
c
      if (anorm .ne. 0.0e0) rcond = ynorm/anorm
      if (anorm .eq. 0.0e0) rcond = 0.0e0
      return
      end
