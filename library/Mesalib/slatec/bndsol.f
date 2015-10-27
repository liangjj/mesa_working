*deck bndsol
      subroutine bndsol (mode, g, mdg, nb, ip, ir, x, n, rnorm)
c***begin prologue  bndsol
c***purpose  solve the least squares problem for a banded matrix using
c            sequential accumulation of rows of the data matrix.
c            exactly one right-hand side vector is permitted.
c***library   slatec
c***category  d9
c***type      single precision (bndsol-s, dbndsl-d)
c***keywords  banded matrix, curve fitting, least squares
c***author  lawson, c. l., (jpl)
c           hanson, r. j., (snla)
c***description
c
c     these subroutines solve the least squares problem ax = b for
c     banded matrices a using sequential accumulation of rows of the
c     data matrix.  exactly one right-hand side vector is permitted.
c
c     these subroutines are intended for the type of least squares
c     systems that arise in applications such as curve or surface
c     fitting of data.  the least squares equations are accumulated and
c     processed using only part of the data.  this requires a certain
c     user interaction during the solution of ax = b.
c
c     specifically, suppose the data matrix (a b) is row partitioned
c     into q submatrices.  let (e f) be the t-th one of these
c     submatrices where e = (0 c 0).  here the dimension of e is mt by n
c     and the dimension of c is mt by nb.  the value of nb is the
c     bandwidth of a.  the dimensions of the leading block of zeros in e
c     are mt by jt-1.
c
c     the user of the subroutine bndacc provides mt,jt,c and f for
c     t=1,...,q.  not all of this data must be supplied at once.
c
c     following the processing of the various blocks (e f), the matrix
c     (a b) has been transformed to the form (r d) where r is upper
c     triangular and banded with bandwidth nb.  the least squares
c     system rx = d is then easily solved using back substitution by
c     executing the statement call bndsol(1,...). the sequence of
c     values for jt must be nondecreasing.  this may require some
c     preliminary interchanges of rows and columns of the matrix a.
c
c     the primary reason for these subroutines is that the total
c     processing can take place in a working array of dimension mu by
c     nb+1.  an acceptable value for mu is
c
c                       mu = max(mt + n + 1),
c
c     where n is the number of unknowns.
c
c     here the maximum is taken over all values of mt for t=1,...,q.
c     notice that mt can be taken to be a small as one, showing that
c     mu can be as small as n+2.  the subprogram bndacc processes the
c     rows more efficiently if mu is large enough so that each new
c     block (c f) has a distinct value of jt.
c
c     the four principle parts of these algorithms are obtained by the
c     following call statements
c
c     call bndacc(...)  introduce new blocks of data.
c
c     call bndsol(1,...)compute solution vector and length of
c                       residual vector.
c
c     call bndsol(2,...)given any row vector h solve yr = h for the
c                       row vector y.
c
c     call bndsol(3,...)given any column vector w solve rz = w for
c                       the column vector z.
c
c     the dots in the above call statements indicate additional
c     arguments that will be specified in the following paragraphs.
c
c     the user must dimension the array appearing in the call list..
c     g(mdg,nb+1)
c
c     description of calling sequence for bndacc..
c
c     the entire set of parameters for bndacc are
c
c     input..
c
c     g(*,*)            the working array into which the user will
c                       place the mt by nb+1 block (c f) in rows ir
c                       through ir+mt-1, columns 1 through nb+1.
c                       see descriptions of ir and mt below.
c
c     mdg               the number of rows in the working array
c                       g(*,*).  the value of mdg should be .ge. mu.
c                       the value of mu is defined in the abstract
c                       of these subprograms.
c
c     nb                the bandwidth of the data matrix a.
c
c     ip                set by the user to the value 1 before the
c                       first call to bndacc.  its subsequent value
c                       is controlled by bndacc to set up for the
c                       next call to bndacc.
c
c     ir                index of the row of g(*,*) where the user is
c                       the user to the value 1 before the first call
c                       to bndacc.  its subsequent value is controlled
c                       by bndacc. a value of ir .gt. mdg is considered
c                       an error.
c
c     mt,jt             set by the user to indicate respectively the
c                       number of new rows of data in the block and
c                       the index of the first nonzero column in that
c                       set of rows (e f) = (0 c 0 f) being processed.
c     output..
c
c     g(*,*)            the working array which will contain the
c                       processed rows of that part of the data
c                       matrix which has been passed to bndacc.
c
c     ip,ir             the values of these arguments are advanced by
c                       bndacc to be ready for storing and processing
c                       a new block of data in g(*,*).
c
c     description of calling sequence for bndsol..
c
c     the user must dimension the arrays appearing in the call list..
c
c     g(mdg,nb+1), x(n)
c
c     the entire set of parameters for bndsol are
c
c     input..
c
c     mode              set by the user to one of the values 1, 2, or
c                       3.  these values respectively indicate that
c                       the solution of ax = b, yr = h or rz = w is
c                       required.
c
c     g(*,*),mdg,       these arguments all have the same meaning and
c      nb,ip,ir         contents as following the last call to bndacc.
c
c     x(*)              with mode=2 or 3 this array contains,
c                       respectively, the right-side vectors h or w of
c                       the systems yr = h or rz = w.
c
c     n                 the number of variables in the solution
c                       vector.  if any of the n diagonal terms are
c                       zero the subroutine bndsol prints an
c                       appropriate message.  this condition is
c                       considered an error.
c
c     output..
c
c     x(*)              this array contains the solution vectors x,
c                       y or z of the systems ax = b, yr = h or
c                       rz = w depending on the value of mode=1,
c                       2 or 3.
c
c     rnorm             if mode=1 rnorm is the euclidean length of the
c                       residual vector ax-b.  when mode=2 or 3 rnorm
c                       is set to zero.
c
c     remarks..
c
c     to obtain the upper triangular matrix and transformed right-hand
c     side vector d so that the super diagonals of r form the columns
c     of g(*,*), execute the following fortran statements.
c
c     nbp1=nb+1
c
c     do 10 j=1, nbp1
c
c  10 g(ir,j) = 0.e0
c
c     mt=1
c
c     jt=n+1
c
c     call bndacc(g,mdg,nb,ip,ir,mt,jt)
c
c***references  c. l. lawson and r. j. hanson, solving least squares
c                 problems, prentice-hall, inc., 1974, chapter 27.
c***routines called  xermsg
c***revision history  (yymmdd)
c   790101  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   891006  cosmetic changes to prologue.  (wrb)
c   891006  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  bndsol
      dimension g(mdg,*),x(*)
c***first executable statement  bndsol
      zero=0.
c
      rnorm=zero
      go to (10,90,50), mode
c                                   ********************* mode = 1
c                                   alg. step 26
   10      do 20 j=1,n
           x(j)=g(j,nb+1)
   20 continue
      rsq=zero
      np1=n+1
      irm1=ir-1
      if (np1.gt.irm1) go to 40
           do 30 j=np1,irm1
           rsq=rsq+g(j,nb+1)**2
   30 continue
      rnorm=sqrt(rsq)
   40 continue
c                                   ********************* mode = 3
c                                   alg. step 27
   50      do 80 ii=1,n
           i=n+1-ii
c                                   alg. step 28
           s=zero
           l=max(0,i-ip)
c                                   alg. step 29
           if (i.eq.n) go to 70
c                                   alg. step 30
           ie=min(n+1-i,nb)
                do 60 j=2,ie
                jg=j+l
                ix=i-1+j
                s=s+g(i,jg)*x(ix)
   60 continue
c                                   alg. step 31
   70      if (g(i,l+1)) 80,130,80
   80      x(i)=(x(i)-s)/g(i,l+1)
c                                   alg. step 32
      return
c                                   ********************* mode = 2
   90      do 120 j=1,n
           s=zero
           if (j.eq.1) go to 110
           i1=max(1,j-nb+1)
           i2=j-1
                do 100 i=i1,i2
                l=j-i+1+max(0,i-ip)
                s=s+x(i)*g(i,l)
  100 continue
  110      l=max(0,j-ip)
           if (g(j,l+1)) 120,130,120
  120      x(j)=(x(j)-s)/g(j,l+1)
      return
c
  130 continue
      nerr=1
      iopt=2
      call xermsg ('slatec', 'bndsol',
     +   'a zero diagonal term is in the n by n upper triangular ' //
     +   'matrix.', nerr, iopt)
      return
      end
