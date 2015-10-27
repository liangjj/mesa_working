*deck sglss
      subroutine sglss (a, mda, m, n, b, mdb, nb, rnorm, work, lw,
     +   iwork, liw, info)
c***begin prologue  sglss
c***purpose  solve a linear least squares problems by performing a qr
c            factorization of the matrix using householder
c            transformations.  emphasis is put on detecting possible
c            rank deficiency.
c***library   slatec
c***category  d9, d5
c***type      single precision (sglss-s, dglss-d)
c***keywords  linear least squares, lq factorization, qr factorization,
c             underdetermined linear systems
c***author  manteuffel, t. a., (lanl)
c***description
c
c     sglss solves both underdetermined and overdetermined
c     linear systems ax = b, where a is an m by n matrix
c     and b is an m by nb matrix of right hand sides. if
c     m.ge.n, the least squares solution is computed by
c     decomposing the matrix a into the product of an
c     orthogonal matrix q and an upper triangular matrix
c     r (qr factorization). if m.lt.n, the minimal
c     length solution is computed by factoring the
c     matrix a into the product of a lower triangular
c     matrix l and an orthogonal matrix q (lq factor-
c     ization). if the matrix a is determined to be rank
c     deficient, that is the rank of a is less than
c     min(m,n), then the minimal length least squares
c     solution is computed.
c
c     sglss assumes full machine precision in the data.
c     if more control over the uncertainty in the data
c     is desired, the codes llsia and ulsia are
c     recommended.
c
c     sglss requires mda*n + (mdb + 1)*nb + 5*min(m,n) dimensioned
c     real space and m+n dimensioned integer space.
c
c
c   ******************************************************************
c   *                                                                *
c   *         warning - all input arrays are changed on exit.        *
c   *                                                                *
c   ******************************************************************
c     subroutine sglss(a,mda,m,n,b,mdb,nb,rnorm,work,lw,iwork,liw,info)
c
c     input..
c
c     a(,)          linear coefficient matrix of ax=b, with mda the
c      mda,m,n      actual first dimension of a in the calling program.
c                   m is the row dimension (no. of equations of the
c                   problem) and n the col dimension (no. of unknowns).
c
c     b(,)          right hand side(s), with mdb the actual first
c      mdb,nb       dimension of b in the calling program. nb is the
c                   number of m by 1 right hand sides. must have
c                   mdb.ge.max(m,n). if nb = 0, b is never accessed.
c
c
c     rnorm()       vector of length at least nb.  on input the contents
c                   of rnorm are unused.
c
c     work()        a real work array dimensioned 5*min(m,n).
c
c     lw            actual dimension of work.
c
c     iwork()       integer work array dimensioned at least n+m.
c
c     liw           actual dimension of iwork.
c
c
c     info          a flag which provides for the efficient
c                   solution of subsequent problems involving the
c                   same a but different b.
c                   if info = 0 original call
c                      info = 1 subsequent calls
c                   on subsequent calls, the user must supply a, info,
c                   lw, iwork, liw, and the first 2*min(m,n) locations
c                   of work as output by the original call to sglss.
c
c
c     output..
c
c     a(,)          contains the triangular part of the reduced matrix
c                   and the transformation information. it together with
c                   the first 2*min(m,n) elements of work (see below)
c                   completely specify the factorization of a.
c
c     b(,)          contains the n by nb solution matrix x.
c
c
c     rnorm()       contains the euclidean length of the nb residual
c                   vectors  b(i)-ax(i), i=1,nb.
c
c     work()        the first 2*min(m,n) locations of work contain value
c                   necessary to reproduce the factorization of a.
c
c     iwork()       the first m+n locations contain the order in
c                   which the rows and columns of a were used.
c                   if m.ge.n columns then rows. if m.lt.n rows
c                   then columns.
c
c     info          flag to indicate status of computation on completion
c                  -1   parameter error(s)
c                   0 - full rank
c                   n.gt.0 - reduced rank  rank=min(m,n)-info
c
c***references  t. manteuffel, an interval analysis approach to rank
c                 determination in linear least squares problems,
c                 report sand80-0655, sandia laboratories, june 1980.
c***routines called  llsia, ulsia
c***revision history  (yymmdd)
c   810801  date written
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  sglss
      dimension a(mda,*),b(mdb,*),rnorm(*),work(*)
      integer iwork(*)
c
c***first executable statement  sglss
      re=0.
      ae=0.
      key=0
      mode=2
      np=0
c
c     if m.ge.n call llsia
c     if m.lt.n call ulsia
c
      if(m.lt.n) go to 10
      call llsia(a,mda,m,n,b,mdb,nb,re,ae,key,mode,np,
     1            krank,ksure,rnorm,work,lw,iwork,liw,info)
      if(info.eq.-1) return
      info=n-krank
      return
   10 call ulsia(a,mda,m,n,b,mdb,nb,re,ae,key,mode,np,
     1            krank,ksure,rnorm,work,lw,iwork,liw,info)
      if(info.eq.-1) return
      info=m-krank
      return
      end
