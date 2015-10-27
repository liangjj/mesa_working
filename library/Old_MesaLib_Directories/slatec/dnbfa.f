*deck dnbfa
      subroutine dnbfa (abe, lda, n, ml, mu, ipvt, info)
c***begin prologue  dnbfa
c***purpose  factor a band matrix by elimination.
c***library   slatec
c***category  d2a2
c***type      double precision (snbfa-s, dnbfa-d, cnbfa-c)
c***keywords  banded, linear equations, matrix factorization,
c             nonsymmetric
c***author  voorhees, e. a., (lanl)
c***description
c
c     dnbfa factors a double precision band matrix by elimination.
c
c     dnbfa is usually called by dnbco, but it can be called
c     directly with a saving in time if rcond is not needed.
c
c     on entry
c
c        abe     double precision(lda, nc)
c                contains the matrix in band storage.  the rows
c                of the original matrix are stored in the rows
c                of abe and the diagonals of the original matrix
c                are stored in columns 1 through ml+mu+1 of abe.
c                nc must be .ge. 2*ml+mu+1 .
c                see the comments below for details.
c
c        lda     integer
c                the leading dimension of the array abe.
c                lda must be .ge.  n .
c
c        n       integer
c                the order of the original matrix.
c
c        ml      integer
c                number of diagonals below the main diagonal.
c                0 .le. ml .lt.  n .
c
c        mu      integer
c                number of diagonals above the main diagonal.
c                0 .le. mu .lt.  n .
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
c        info    integer
c                =0  normal value
c                =k  if  u(k,k) .eq. 0.0 .  this is not an error
c                condition for this subroutine, but it does
c                indicate that dnbsl will divide by zero if
c                called.  use rcond in dnbco for a reliable
c                indication of singularity.
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
c***routines called  daxpy, dscal, dswap, idamax
c***revision history  (yymmdd)
c   800728  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  dnbfa
      integer lda,n,ml,mu,ipvt(*),info
      double precision abe(lda,*)
c
      integer ml1,mb,m,n1,ldb,i,j,k,l,lm,lm1,lm2,mp,idamax
      double precision t
c***first executable statement  dnbfa
      ml1=ml+1
      mb=ml+mu
      m=ml+mu+1
      n1=n-1
      ldb=lda-1
      info=0
c
c     set fill-in columns to zero
c
      if(n.le.1)go to 50
      if(ml.le.0)go to 7
      do 6 j=1,ml
        do 5 i=1,n
          abe(i,m+j)=0.0d0
    5   continue
    6 continue
    7 continue
c
c     gaussian elimination with partial elimination
c
      do 40 k=1,n1
        lm=min(n-k,ml)
        lm1=lm+1
        lm2=ml1-lm
c
c     search for pivot index
c
        l=-idamax(lm1,abe(lm+k,lm2),ldb)+lm1+k
        ipvt(k)=l
        mp=min(mb,n-k)
c
c     swap rows if necessary
c
        if(l.ne.k)call dswap(mp+1,abe(k,ml1),lda,abe(l,ml1+k-l),lda)
c
c     skip column reduction if pivot is zero
c
        if(abe(k,ml1).eq.0.0d0) go to 20
c
c     compute multipliers
c
        t=-1.0/abe(k,ml1)
        call dscal(lm,t,abe(lm+k,lm2),ldb)
c
c     row elimination with column indexing
c
        do 10 j=1,mp
          call daxpy (lm,abe(k,ml1+j),abe(lm+k,lm2),ldb,abe(lm+k,lm2+j),
     1                ldb)
   10   continue
        go to 30
   20   continue
        info=k
   30   continue
   40 continue
   50 continue
      ipvt(n)=n
      if(abe(n,ml1).eq.0.0d0) info=n
      return
      end
