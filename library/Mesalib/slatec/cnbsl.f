*deck cnbsl
      subroutine cnbsl (abe, lda, n, ml, mu, ipvt, b, job)
c***begin prologue  cnbsl
c***purpose  solve a complex band system using the factors computed by
c            cnbco or cnbfa.
c***library   slatec
c***category  d2c2
c***type      complex (snbsl-s, dnbsl-d, cnbsl-c)
c***keywords  banded, linear equations, nonsymmetric, solve
c***author  voorhees, e. a., (lanl)
c***description
c
c     cnbsl solves the complex band system
c     a * x = b  or  ctrans(a) * x = b
c     using the factors computed by cnbco or cnbfa.
c
c     on entry
c
c        abe     complex(lda, nc)
c                the output from cnbco or cnbfa.
c                nc must be .ge. 2*ml+mu+1 .
c
c        lda     integer
c                the leading dimension of the array  abe .
c
c        n       integer
c                the order of the original matrix.
c
c        ml      integer
c                number of diagonals below the main diagonal.
c
c        mu      integer
c                number of diagonals above the main diagonal.
c
c        ipvt    integer(n)
c                the pivot vector from cnbco or cnbfa.
c
c        b       complex(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b .
c                = nonzero   to solve  ctrans(a)*x = b , where
c                            ctrans(a)  is the conjugate transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda.  it will not occur if the subroutines are
c        called correctly and if cnbco has set rcond .gt. 0.0
c        or cnbfa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call cnbco(abe,lda,n,ml,mu,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c             call cnbsl(abe,lda,n,ml,mu,ipvt,c(1,j),0)
c        10 continue
c
c***references  j. j. dongarra, j. r. bunch, c. b. moler, and g. w.
c                 stewart, linpack users' guide, siam, 1979.
c***routines called  caxpy, cdotc
c***revision history  (yymmdd)
c   800730  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  cnbsl
      integer lda,n,ml,mu,ipvt(*),job
      complex abe(lda,*),b(*)
c
      complex cdotc,t
      integer k,kb,l,lb,ldb,lm,m,mlm,nm1
c***first executable statement  cnbsl
      m=mu+ml+1
      nm1=n-1
      ldb=1-lda
      if(job.ne.0)go to 50
c
c       job = 0 , solve  a * x = b
c       first solve l*y = b
c
        if(ml.eq.0)go to 30
        if(nm1.lt.1)go to 30
          do 20 k=1,nm1
            lm=min(ml,n-k)
            l=ipvt(k)
            t=b(l)
            if(l.eq.k)go to 10
              b(l)=b(k)
              b(k)=t
   10       continue
            mlm=ml-(lm-1)
            call caxpy(lm,t,abe(k+lm,mlm),ldb,b(k+1),1)
   20     continue
   30   continue
c
c       now solve  u*x = y
c
        do 40 kb=1,n
          k=n+1-kb
          b(k)=b(k)/abe(k,ml+1)
          lm=min(k,m)-1
          lb=k-lm
          t=-b(k)
          call caxpy(lm,t,abe(k-1,ml+2),ldb,b(lb),1)
   40   continue
      go to 100
   50 continue
c
c       job = nonzero, solve ctrans(a) * x = b
c       first solve  ctrans(u)*y = b
c
        do 60 k = 1, n
          lm = min(k,m) - 1
          lb = k - lm
          t = cdotc(lm,abe(k-1,ml+2),ldb,b(lb),1)
          b(k) = (b(k) - t)/conjg(abe(k,ml+1))
   60   continue
c
c       now solve ctrans(l)*x = y
c
        if (ml .eq. 0) go to 90
        if (nm1 .lt. 1) go to 90
          do 80 kb = 1, nm1
            k = n - kb
            lm = min(ml,n-k)
            mlm = ml - (lm - 1)
            b(k) = b(k) + cdotc(lm,abe(k+lm,mlm),ldb,b(k+1),1)
            l = ipvt(k)
            if (l .eq. k) go to 70
              t = b(l)
              b(l) = b(k)
              b(k) = t
   70       continue
   80     continue
   90   continue
  100 continue
      return
      end