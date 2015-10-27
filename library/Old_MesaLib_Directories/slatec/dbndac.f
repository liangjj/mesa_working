*deck dbndac
      subroutine dbndac (g, mdg, nb, ip, ir, mt, jt)
c***begin prologue  dbndac
c***purpose  compute the lu factorization of a  banded matrices using
c            sequential accumulation of rows of the data matrix.
c            exactly one right-hand side vector is permitted.
c***library   slatec
c***category  d9
c***type      double precision (bndacc-s, dbndac-d)
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
c     the user of the subroutine dbndac provides mt,jt,c and f for
c     t=1,...,q.  not all of this data must be supplied at once.
c
c     following the processing of the various blocks (e f), the matrix
c     (a b) has been transformed to the form (r d) where r is upper
c     triangular and banded with bandwidth nb.  the least squares
c     system rx = d is then easily solved using back substitution by
c     executing the statement call dbndsl(1,...). the sequence of
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
c     mu can be as small as n+2.  the subprogram dbndac processes the
c     rows more efficiently if mu is large enough so that each new
c     block (c f) has a distinct value of jt.
c
c     the four principle parts of these algorithms are obtained by the
c     following call statements
c
c     call dbndac(...)  introduce new blocks of data.
c
c     call dbndsl(1,...)compute solution vector and length of
c                       residual vector.
c
c     call dbndsl(2,...)given any row vector h solve yr = h for the
c                       row vector y.
c
c     call dbndsl(3,...)given any column vector w solve rz = w for
c                       the column vector z.
c
c     the dots in the above call statements indicate additional
c     arguments that will be specified in the following paragraphs.
c
c     the user must dimension the array appearing in the call list..
c     g(mdg,nb+1)
c
c     description of calling sequence for dbndac..
c
c     the entire set of parameters for dbndac are
c
c     input.. all type real variables are double precision
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
c                       first call to dbndac.  its subsequent value
c                       is controlled by dbndac to set up for the
c                       next call to dbndac.
c
c     ir                index of the row of g(*,*) where the user is
c                       to place the new block of data (c f).  set by
c                       the user to the value 1 before the first call
c                       to dbndac.  its subsequent value is controlled
c                       by dbndac. a value of ir .gt. mdg is considered
c                       an error.
c
c     mt,jt             set by the user to indicate respectively the
c                       number of new rows of data in the block and
c                       the index of the first nonzero column in that
c                       set of rows (e f) = (0 c 0 f) being processed.
c
c     output.. all type real variables are double precision
c
c     g(*,*)            the working array which will contain the
c                       processed rows of that part of the data
c                       matrix which has been passed to dbndac.
c
c     ip,ir             the values of these arguments are advanced by
c                       dbndac to be ready for storing and processing
c                       a new block of data in g(*,*).
c
c     description of calling sequence for dbndsl..
c
c     the user must dimension the arrays appearing in the call list..
c
c     g(mdg,nb+1), x(n)
c
c     the entire set of parameters for dbndsl are
c
c     input.. all type real variables are double precision
c
c     mode              set by the user to one of the values 1, 2, or
c                       3.  these values respectively indicate that
c                       the solution of ax = b, yr = h or rz = w is
c                       required.
c
c     g(*,*),mdg,       these arguments all have the same meaning and
c      nb,ip,ir         contents as following the last call to dbndac.
c
c     x(*)              with mode=2 or 3 this array contains,
c                       respectively, the right-side vectors h or w of
c                       the systems yr = h or rz = w.
c
c     n                 the number of variables in the solution
c                       vector.  if any of the n diagonal terms are
c                       zero the subroutine dbndsl prints an
c                       appropriate message.  this condition is
c                       considered an error.
c
c     output.. all type real variables are double precision
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
c     call dbndac(g,mdg,nb,ip,ir,mt,jt)
c
c***references  c. l. lawson and r. j. hanson, solving least squares
c                 problems, prentice-hall, inc., 1974, chapter 27.
c***routines called  dh12, xermsg
c***revision history  (yymmdd)
c   790101  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891006  cosmetic changes to prologue.  (wrb)
c   891006  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   920501  reformatted the references section.  (wrb)
c***end prologue  dbndac
      implicit double precision (a-h,o-z)
      dimension g(mdg,*)
c***first executable statement  dbndac
      zero=0.d0
c
c              alg. steps 1-4 are performed external to this subroutine.
c
      nbp1=nb+1
      if (mt.le.0.or.nb.le.0) return
c
      if(.not.mdg.lt.ir) go to 5
      nerr=1
      iopt=2
      call xermsg ('slatec', 'dbndac', 'mdg.lt.ir, probable error.',
     +   nerr, iopt)
      return
    5 continue
c
c                                             alg. step 5
      if (jt.eq.ip) go to 70
c                                             alg. steps 6-7
      if (jt.le.ir) go to 30
c                                             alg. steps 8-9
      do 10 i=1,mt
        ig1=jt+mt-i
        ig2=ir+mt-i
        do 10 j=1,nbp1
        g(ig1,j)=g(ig2,j)
   10 continue
c                                             alg. step 10
      ie=jt-ir
      do 20 i=1,ie
        ig=ir+i-1
        do 20 j=1,nbp1
        g(ig,j)=zero
   20 continue
c                                             alg. step 11
      ir=jt
c                                             alg. step 12
   30 mu=min(nb-1,ir-ip-1)
      if (mu.eq.0) go to 60
c                                             alg. step 13
      do 50 l=1,mu
c                                             alg. step 14
        k=min(l,jt-ip)
c                                             alg. step 15
        lp1=l+1
        ig=ip+l
        do 40 i=lp1,nb
          jg=i-k
          g(ig,jg)=g(ig,i)
   40 continue
c                                             alg. step 16
        do 50 i=1,k
        jg=nbp1-i
        g(ig,jg)=zero
   50 continue
c                                             alg. step 17
   60 ip=jt
c                                             alg. steps 18-19
   70 mh=ir+mt-ip
      kh=min(nbp1,mh)
c                                             alg. step 20
      do 80 i=1,kh
        call dh12 (1,i,max(i+1,ir-ip+1),mh,g(ip,i),1,rho,
     1            g(ip,i+1),1,mdg,nbp1-i)
   80 continue
c                                             alg. step 21
      ir=ip+kh
c                                             alg. step 22
      if (kh.lt.nbp1) go to 100
c                                             alg. step 23
      do 90 i=1,nb
        g(ir-1,i)=zero
   90 continue
c                                             alg. step 24
  100 continue
c                                             alg. step 25
      return
      end
