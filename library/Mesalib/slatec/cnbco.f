*deck cnbco
      subroutine cnbco (abe, lda, n, ml, mu, ipvt, rcond, z)
c***begin prologue  cnbco
c***purpose  factor a band matrix using gaussian elimination and
c            estimate the condition number.
c***library   slatec
c***category  d2c2
c***type      complex (snbco-s, dnbco-d, cnbco-c)
c***keywords  banded, linear equations, matrix factorization,
c             nonsymmetric
c***author  voorhees, e. a., (lanl)
c***description
c
c     cnbco factors a complex band matrix by gaussian
c     elimination and estimates the condition of the matrix.
c
c     if rcond is not needed, cnbfa is slightly faster.
c     to solve  a*x = b , follow cnbco by cnbsl.
c     to compute  inverse(a)*c , follow cnbco by cnbsl.
c     to compute  determinant(a) , follow cnbco by cnbdi.
c
c     on entry
c
c        abe     complex(lda, nc)
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
c                the factorization can be written  a = l*u  where
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
c        z       complex(n)
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
c***routines called  caxpy, cdotc, cnbfa, csscal, scasum
c***revision history  (yymmdd)
c   800730  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  cnbco
      integer lda,n,ml,mu,ipvt(*)
      complex abe(lda,*),z(*)
      real rcond
c
      complex cdotc,ek,t,wk,wkm
      real anorm,s,scasum,sm,ynorm
      integer i,info,j,ju,k,kb,kp1,l,ldb,lm,lz,m,ml1,mm,nl,nu
      complex zdum,zdum1,zdum2,csign1
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
      csign1(zdum1,zdum2) = cabs1(zdum1)*(zdum2/cabs1(zdum2))
c
c     compute 1-norm of a
c
c***first executable statement  cnbco
      ml1=ml+1
      ldb = lda - 1
      anorm = 0.0e0
      do 10 j = 1, n
        nu = min(mu,j-1)
        nl = min(ml,n-j)
        l = 1 + nu + nl
        anorm = max(anorm,scasum(l,abe(j+nl,ml1-nl),ldb))
   10 continue
c
c     factor
c
      call cnbfa(abe,lda,n,ml,mu,ipvt,info)
c
c     rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) .
c     estimate = norm(z)/norm(y) where  a*z = y  and ctrans(a)*y = e .
c     ctrans(a)  is the conjugate transpose of a .
c     the components of  e  are chosen to cause maximum local
c     growth in the elements of  w  where ctrans(u)*w = e .
c     the vectors are frequently rescaled to avoid overflow.
c
c     solve ctrans(u)*w = e
c
      ek = (1.0e0,0.0e0)
      do 20 j = 1, n
        z(j) = (0.0e0,0.0e0)
   20 continue
      m = ml + mu + 1
      ju = 0
      do 100 k = 1, n
        if (cabs1(z(k)) .ne. 0.0e0) ek = csign1(ek,-z(k))
        if (cabs1(ek-z(k)) .le. cabs1(abe(k,ml1))) go to 30
          s = cabs1(abe(k,ml1))/cabs1(ek-z(k))
          call csscal(n,s,z,1)
          ek = cmplx(s,0.0e0)*ek
   30   continue
        wk = ek - z(k)
        wkm = -ek - z(k)
        s = cabs1(wk)
        sm = cabs1(wkm)
        if (cabs1(abe(k,ml1)) .eq. 0.0e0) go to 40
          wk = wk/conjg(abe(k,ml1))
          wkm = wkm/conjg(abe(k,ml1))
        go to 50
   40   continue
          wk = (1.0e0,0.0e0)
          wkm = (1.0e0,0.0e0)
   50   continue
        kp1 = k + 1
        ju = min(max(ju,mu+ipvt(k)),n)
        mm = ml1
        if (kp1 .gt. ju) go to 90
          do 60 i = kp1, ju
            mm = mm + 1
            sm = sm + cabs1(z(i)+wkm*conjg(abe(k,mm)))
            z(i) = z(i) + wk*conjg(abe(k,mm))
            s = s + cabs1(z(i))
   60     continue
          if (s .ge. sm) go to 80
            t = wkm -wk
            wk = wkm
            mm = ml1
            do 70 i = kp1, ju
              mm = mm + 1
              z(i) = z(i) + t*conjg(abe(k,mm))
   70       continue
   80     continue
   90   continue
      z(k) = wk
  100 continue
      s = 1.0e0/scasum(n,z,1)
      call csscal(n,s,z,1)
c
c     solve ctrans(l)*y = w
c
      do 120 kb = 1, n
        k = n + 1 - kb
        nl = min(ml,n-k)
        if (k .lt. n) z(k) = z(k) + cdotc(nl,abe(k+nl,ml1-nl),-ldb,
     1  z(k+1),1)
        if (cabs1(z(k)) .le. 1.0e0) go to 110
          s = 1.0e0/cabs1(z(k))
          call csscal(n,s,z,1)
  110   continue
        l = ipvt(k)
        t = z(l)
        z(l) = z(k)
        z(k) = t
  120 continue
      s = 1.0e0/scasum(n,z,1)
      call csscal(n,s,z,1)
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
        if (k .lt. n) call caxpy(nl,t,abe(k+nl,ml1-nl),-ldb,z(k+1),1)
        if (cabs1(z(k)) .le. 1.0e0) go to 130
          s = 1.0e0/cabs1(z(k))
          call csscal(n,s,z,1)
          ynorm = s*ynorm
  130   continue
  140 continue
      s = 1.0e0/scasum(n,z,1)
      call csscal(n,s,z,1)
      ynorm = s*ynorm
c
c     solve  u*z = v
c
      do 160 kb = 1, n
        k = n + 1 - kb
        if (cabs1(z(k)) .le. cabs1(abe(k,ml1))) go to 150
          s = cabs1(abe(k,ml1))/cabs1(z(k))
          call csscal(n,s,z,1)
          ynorm = s*ynorm
  150   continue
        if (cabs1(abe(k,ml1)) .ne. 0.0e0) z(k) = z(k)/abe(k,ml1)
        if (cabs1(abe(k,ml1)) .eq. 0.0e0) z(k) = 1.0e0
        lm = min(k,m) - 1
        lz = k - lm
        t = -z(k)
        call caxpy(lm,t,abe(k-1,ml+2),-ldb,z(lz),1)
  160 continue
c     make znorm = 1.0e0
      s = 1.0e0/scasum(n,z,1)
      call csscal(n,s,z,1)
      ynorm = s*ynorm
c
      if (anorm .ne. 0.0e0) rcond = ynorm/anorm
      if (anorm .eq. 0.0e0) rcond = 0.0e0
      return
      end
