*deck dlpdoc
      subroutine dlpdoc
c***begin prologue  dlpdoc
c***purpose  sparse linear algebra package version 2.0.2 documentation.
c            routines to solve large sparse symmetric and nonsymmetric
c            positive definite linear systems, ax = b, using precondi-
c            tioned iterative methods.
c***library   slatec (slap)
c***category  d2a4, d2b4, z
c***type      double precision (slpdoc-s, dlpdoc-d)
c***keywords  biconjugate gradient squared, documentation,
c             generalized minimum residual, iterative improvement,
c             normal equations, orthomin,
c             preconditioned conjugate gradient, slap,
c             sparse iterative methods
c***author  seager, mark. k., (llnl)
c             user systems division
c             lawrence livermore national laboratory
c             po box 808, l-60
c             livermore, ca 94550
c             (fts) 543-3141, (510) 423-3141
c             seager@llnl.gov
c***description
c                                 the
c                    sparse linear algebra package
c                      double precision routines
c
c                @@@@@@@  @            @@@    @@@@@@@@
c               @       @ @           @   @   @       @
c               @         @          @     @  @       @
c                @@@@@@@  @         @       @ @@@@@@@@
c                       @ @         @@@@@@@@@ @
c               @       @ @         @       @ @
c                @@@@@@@  @@@@@@@@@ @       @ @
c
c      @       @                            @@@@@@@        @@@@@
c      @       @                           @       @      @    @@
c      @       @  @@@@@@@  @ @@                    @     @    @  @
c      @       @ @       @ @@  @             @@@@@@      @   @   @
c       @     @  @@@@@@@@@ @                @            @  @    @
c        @   @   @         @               @         @@@  @@    @
c         @@@     @@@@@@@  @               @@@@@@@@@ @@@   @@@@@
c
c
c    =================================================================
c    ========================== introduction =========================
c    =================================================================
c      this package was  originally derived from a set of  iterative
c      routines written by anne greenbaum, as announced in "routines
c      for solving large sparse linear systems",  tentacle, lawrence
c      livermore  national  laboratory,  livermore  computing center
c      (january 1986), pp 15-21.
c
c    this document  contains the specifications for  the  slap version
c    2.0 package, a fortran 77  package  for  the  solution  of  large
c    sparse   linear systems, ax  =  b,  via  preconditioned iterative
c    methods.   included in  this  package are "core"  routines  to do
c    iterative   refinement  (jacobi's  method),  conjugate  gradient,
c    conjugate gradient on the normal equations, aa'y = b,  (where x =
c    a'y and  a' denotes the  transpose of   a), biconjugate gradient,
c    biconjugate  gradient  squared, orthomin and  generalized minimum
c    residual iteration.    these "core" routines   do  not  require a
c    "fixed"   data  structure   for storing  the   matrix  a  and the
c    preconditioning   matrix  m.   the  user  is free  to  choose any
c    structure that facilitates  efficient solution  of the problem at
c    hand.  the drawback  to this approach  is that the user must also
c    supply at least two routines  (matvec and msolve,  say).   matvec
c    must calculate, y = ax, given x and the user's data structure for
c    a.  msolve must solve,  r = mz, for z (*not*  r) given r  and the
c    user's data  structure for  m (or its  inverse).  the user should
c    choose m so that  inv(m)*a  is approximately the identity and the
c    solution step r = mz is "easy" to  solve.  for some of the "core"
c    routines (orthomin,  biconjugate gradient and  conjugate gradient
c    on the  normal equations)   the user must  also  supply  a matrix
c    transpose times   vector  routine  (mttvec,  say)  and (possibly,
c    depending    on the "core"  method)   a  routine  that solves the
c    transpose  of   the   preconditioning    step     (mtsolv,  say).
c    specifically, mttvec is a routine which calculates y = a'x, given
c    x and the user's data structure for a (a' is the transpose of a).
c    mtsolv is a routine which solves the system r = m'z for z given r
c    and the user's data structure for m.
c
c    this process of writing the matrix vector operations  can be time
c    consuming and error  prone.  to alleviate  these problems we have
c    written drivers   for  the  "core" methods  that  assume the user
c    supplies one of two specific data structures (slap triad and slap
c    column format), see  below.  utilizing these  data structures  we
c    have augmented   each  "core" method  with   two preconditioners:
c    diagonal  scaling and incomplete factorization.  diagonal scaling
c    is easy to implement, vectorizes very  well and for problems that
c    are  not too  ill-conditioned  reduces the  number  of iterations
c    enough   to warrant its use.  on   the other  hand, an incomplete
c    factorization  (incomplete  cholesky for  symmetric systems   and
c    incomplete lu for nonsymmetric  systems) may  take much longer to
c    calculate, but it reduces the iteration count (for most problems)
c    significantly.  our implementations  of ic and ilu  vectorize for
c    machines with hardware gather scatter, but the vector lengths can
c    be quite short if  the  number  of non-zeros  in a column is  not
c    large.
c
c    =================================================================
c    ==================== supplied data structures ===================
c    =================================================================
c    the following describes the data   structures supplied  with  the
c    package: slap triad and column formats.
c
c    ====================== s l a p triad format =====================
c
c    in the slap triad format only the non-zeros are stored.  they may
c    appear in *any* order.  the user supplies three  arrays of length
c    nelt, where nelt  is the   number of  non-zeros  in the   matrix:
c    (ia(nelt),  ja(nelt), a(nelt)).  if  the matrix is symmetric then
c    one need only store the lower triangle (including  the  diagonal)
c    and nelt would be the corresponding  number  of non-zeros stored.
c    for each non-zero the user puts the row and column  index of that
c    matrix  element   in the  ia  and ja  arrays.  the  value  of the
c    non-zero matrix element is placed  in  the corresponding location
c    of  the a array.   this  is an extremely  easy  data structure to
c    generate.  on the other hand, it is not very  efficient on vector
c    computers for the iterative  solution of  linear systems.  hence,
c    slap changes this input data structure to  the slap column format
c    for the iteration (but does not change it back).
c
c    here  is an example   of  the  slap  triad storage  format  for a
c    nonsymmetric 5x5 matrix.  nelt=11.   recall that the  entries may
c    appear in any order.
c
c     5x5 matrix       slap triad format for 5x5 matrix on left.
c                           1  2  3  4  5  6  7  8  9 10 11
c    |11 12  0  0 15|   a: 51 12 11 33 15 53 55 22 35 44 21
c    |21 22  0  0  0|  ia:  5  1  1  3  1  5  5  2  3  4  2
c    | 0  0 33  0 35|  ja:  1  2  1  3  5  3  5  2  5  4  1
c    | 0  0  0 44  0|
c    |51  0 53  0 55|
c
c    ====================== s l a p column format ====================
c
c    in the slap column format  the non-zeros are stored counting down
c    columns (except for the  diagonal entry,  which must appear first
c    in each "column") and are stored in the double precision array a.
c    in  other words,  for each  column  in the matrix  first put  the
c    diagonal  entry  in a.  then  put in the other  non-zero elements
c    going  down the column  (except the  diagonal)  in order.  the ia
c    array holds the row index  for each non-zero.  the ja array holds
c    the  offsets  into the  ia, a  arrays for the  beginning of  each
c    column. that is, ia(ja(icol)), a(ja(icol)) are the first elements
c    of  the  icol-th  column  in  ia  and  a,  and  ia(ja(icol+1)-1),
c    a(ja(icol+1)-1) are the last elements of the icol-th column. note
c    that we  always have  ja(n+1) = nelt+1, where  n is the number of
c    columns in the matrix  and nelt is the number of non-zeros in the
c    matrix.  if the matrix is symmetric one need only store the lower
c    triangle  (including the diagonal)  and nelt would be the  corre-
c    sponding number of non-zeros stored.
c
c    here is  an  example of the  slap   column storage format  for  a
c    nonsymmetric 5x5 matrix (in the  a and  ia arrays '|' denotes the
c    end of a column):
c
c       5x5 matrix      slap column format for 5x5 matrix on left.
c                           1  2  3    4  5    6  7    8    9 10 11
c    |11 12  0  0 15|   a: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35
c    |21 22  0  0  0|  ia:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
c    | 0  0 33  0 35|  ja:  1  4  6    8  9   12
c    | 0  0  0 44  0|
c    |51  0 53  0 55|
c
c    =================================================================
c    ====================== which method to use ======================
c    =================================================================
c
c                          background
c    in solving a large sparse linear system ax = b using an iterative
c    method, it   is  not necessary to actually   store  the matrix a.
c    rather, what is needed is a procedure  for multiplying the matrix
c    a times a given vector y to obtain the matrix-vector product, ay.
c    slap has been written to take advantage of this fact.  the higher
c    level routines in the package require storage only of the non-zero
c    elements of   a (and  their  positions), and  even this   can  be
c    avoided, if the  user  writes his own subroutine for  multiplying
c    the matrix times a vector  and   calls the lower-level  iterative
c    routines in the package.
c
c    if  the matrix a is ill-conditioned,  then most iterative methods
c    will be slow to converge (if they converge  at all!).  to improve
c    the  convergence  rate,  one  may use  a "matrix  splitting," or,
c    "preconditioning matrix," say, m.  it is then necessary to solve,
c    at each iteration, a linear system  with coefficient matrix m.  a
c    good preconditioner  m should have  two  properties: (1) m should
c    "approximate" a, in the sense that the  matrix inv(m)*a  (or some
c    variant  thereof) is better conditioned  than the original matrix
c    a; and  (2) linear  systems with coefficient  matrix m should  be
c    much easier  to solve  than  the original system with coefficient
c    matrix   a.   preconditioning routines  in the   slap package are
c    separate from the  iterative   routines,  so   that any of    the
c    preconditioners provided in the package,   or one that the   user
c    codes himself, can be used with any of the iterative routines.
c
c                        choice of preconditioner
c    if you  willing   to live with   either the slap triad or  column
c    matrix data structure  you  can then  choose one  of two types of
c    preconditioners   to   use:   diagonal  scaling    or  incomplete
c    factorization.  to  choose   between these two   methods requires
c    knowing  something  about the computer you're going  to run these
c    codes on  and how well incomplete factorization  approximates the
c    inverse of your matrix.
c
c    let us  suppose you have   a scalar  machine.   then,  unless the
c    incomplete factorization is very,  very poor this  is *generally*
c    the method to choose.  it  will reduce the  number of  iterations
c    significantly and is not all  that expensive  to compute.  so  if
c    you have just one  linear system to solve  and  "just want to get
c    the job  done" then try  incomplete factorization first.   if you
c    are thinking of integrating some slap  iterative method into your
c    favorite   "production  code" then  try incomplete  factorization
c    first,  but  also check  to see that  diagonal  scaling is indeed
c    slower for a large sample of test problems.
c
c    let us now suppose  you have  a  vector  computer  with  hardware
c    gather/scatter support (cray x-mp, y-mp, scs-40 or cyber 205, eta
c    10,  eta piper,  convex c-1,  etc.).   then  it is much harder to
c    choose  between the  two  methods.   the  versions  of incomplete
c    factorization in slap do in fact vectorize, but have short vector
c    lengths and the factorization step is relatively  more expensive.
c    hence,  for  most problems (i.e.,  unless  your  problem  is  ill
c    conditioned,  sic!)  diagonal  scaling is  faster,  with its very
c    fast    set up  time    and  vectorized  (with   long    vectors)
c    preconditioning step (even though  it  may take more iterations).
c    if you have several systems (or  right hand sides) to  solve that
c    can  utilize  the  same  preconditioner  then the   cost   of the
c    incomplete factorization can   be  amortized over these  several
c    solutions.  this situation gives more advantage to the incomplete
c    factorization methods.  if  you have  a  vector  machine  without
c    hardware  gather/scatter (cray  1,  cray  2  &  cray 3) then  the
c    advantages for incomplete factorization are even less.
c
c    if you're trying to shoehorn slap into your  favorite "production
c    code" and can not easily generate either the slap triad or column
c    format  then  you are  left  to   your  own  devices in terms  of
c    preconditioning.  also,  you may  find that the   preconditioners
c    supplied with slap are not sufficient  for your problem.  in this
c    situation we would  recommend  that you   talk  with a  numerical
c    analyst  versed in   iterative   methods   about   writing  other
c    preconditioning  subroutines (e.g.,  polynomial  preconditioning,
c    shifted incomplete factorization,  sor  or ssor  iteration).  you
c    can always "roll your own"  by using the "core" iterative methods
c    and supplying your own msolve and matvec (and possibly mtsolv and
c    mttvec) routines.
c
c                          symmetric systems
c    if your matrix is symmetric then you would want to use one of the
c    symmetric system  solvers.    if  your  system  is  also positive
c    definite,   (ax,x) (ax dot  product  with x) is  positive for all
c    non-zero  vectors x,  then use   conjugate gradient (dcg,  dsdcg,
c    dsicsg).  if you're  not sure it's spd   (symmetric and  positive
c    definite)  then try dcg anyway and  if it works, fine.  if you're
c    sure your matrix is not  positive definite  then you  may want to
c    try the iterative refinement   methods  (dir)  or the  gmres code
c    (dgmres) if dir converges too slowly.
c
c                         nonsymmetric systems
c    this   is currently  an  area  of  active research  in  numerical
c    analysis  and   there   are   new  strategies  being   developed.
c    consequently take the following advice with a grain of salt.   if
c    you matrix is positive definite, (ax,x)  (ax  dot product  with x
c    is positive for all non-zero  vectors x), then you can use any of
c    the    methods   for   nonsymmetric   systems (orthomin,   gmres,
c    biconjugate gradient, biconjugate gradient  squared and conjugate
c    gradient applied to the normal equations).  if your system is not
c    too ill conditioned then try  biconjugate gradient squared (bcgs)
c    or gmres (dgmres).  both  of  these methods converge very quickly
c    and do  not require a'  or m' ('  denotes transpose) information.
c    dgmres  does require  some  additional storage,  though.  if  the
c    system is very  ill conditioned  or   nearly positive  indefinite
c    ((ax,x) is positive,  but may be  very small),  then gmres should
c    be the first choice,  but try the  other  methods  if you have to
c    fine tune  the solution process for a  "production code".  if you
c    have a great preconditioner for the normal  equations (i.e., m is
c    an approximation to the inverse of aa' rather than  just  a) then
c    this is not a bad route to travel.  old wisdom would say that the
c    normal equations are a disaster  (since it squares the  condition
c    number of the system and dcg convergence is linked to this number
c    of    infamy), but   some     preconditioners    (like incomplete
c    factorization) can reduce the condition number back below that of
c    the original system.
c
c    =================================================================
c    ======================= naming conventions ======================
c    =================================================================
c    slap  iterative  methods,    matrix vector    and  preconditioner
c    calculation  routines   follow a naming   convention  which, when
c    understood, allows one to determine the iterative method and data
c    structure(s) used.  the  subroutine  naming convention  takes the
c    following form:
c                          p[s][m]desc
c    where
c        p  stands for the precision (or data type) of the routine and
c           is required in all names,
c        s  denotes whether or not the routine requires the slap triad
c           or column format (it does if the second letter of the name
c           is s and does not otherwise),
c        m  stands for the type of preconditioner used (only appears
c           in drivers for "core" routines), and
c     desc  is some number of letters describing the method or purpose
c           of the routine.  the following is a list of the "desc"
c           fields for iterative methods and their meaning:
c             bcg,bc:       biconjugate gradient
c             cg:           conjugate gradient
c             cgn,cn:       conjugate gradient on the normal equations
c             cgs,cs:       biconjugate gradient squared
c             gmres,gmr,gm: generalized minimum residual
c             ir,r:         iterative refinement
c             jac:          jacobi's method
c             gs:           gauss-seidel
c             omn,om:       orthomin
c
c    in the double precision version of slap, all routine names start
c    with a d. the brackets around the s and m designate that these
c    fields are optional.
c
c    here are some examples of the routines:
c    1) dbcg: double precision biconjugate gradient "core" routine.
c       one can deduce that this is a "core" routine, because the s and
c       m fields are missing and biconjugate gradient is an iterative
c       method.
c    2) dsdbcg: double precision, slap data structure bcg with diagonal
c       scaling.
c    3) dslubc: double precision, slap data structure bcg with incom-
c       plete lu factorization as the preconditioning.
c    4) dcg: double precision conjugate gradient "core" routine.
c    5) dsdcg: double precision, slap data structure conjugate gradient
c       with diagonal scaling.
c    6) dsiccg: double precision, slap data structure conjugate gra-
c       dient with incomplete cholesky factorization preconditioning.
c
c
c    =================================================================
c    ===================== user callable routines ====================
c    =================================================================
c    the following is a list of  the "user callable" slap routines and
c    their one line descriptions.  the headers denote  the  file names
c    where the routines can be found, as distributed for unix systems.
c
c    note:  each core routine, dxxx, has a corresponding stop routine,
c         isdxxx.  if the stop routine does not have the specific stop
c         test the user requires (e.g., weighted infinity norm),  then
c         the user should modify the source for isdxxx accordingly.
c
c    ============================= dir.f =============================
c    dir: preconditioned iterative refinement sparse ax = b solver.
c    dsjac: jacobi's method iterative sparse ax = b solver.
c    dsgs: gauss-seidel method iterative sparse ax = b solver.
c    dsilur: incomplete lu iterative refinement sparse ax = b solver.
c
c    ============================= dcg.f =============================
c    dcg: preconditioned conjugate gradient sparse ax=b solver.
c    dsdcg: diagonally scaled conjugate gradient sparse ax=b solver.
c    dsiccg: incomplete cholesky conjugate gradient sparse ax=b solver.
c
c    ============================= dcgn.f ============================
c    dcgn: preconditioned cg sparse ax=b solver for normal equations.
c    dsdcgn: diagonally scaled cg sparse ax=b solver for normal eqn's.
c    dslucn: incomplete lu cg sparse ax=b solver for normal equations.
c
c    ============================= dbcg.f ============================
c    dbcg: preconditioned biconjugate gradient sparse ax = b solver.
c    dsdbcg: diagonally scaled biconjugate gradient sparse ax=b solver.
c    dslubc: incomplete lu biconjugate gradient sparse ax=b solver.
c
c    ============================= dcgs.f ============================
c    dcgs: preconditioned biconjugate gradient squared ax=b solver.
c    dsdcgs: diagonally scaled cgs sparse ax=b solver.
c    dslucs: incomplete lu biconjugate gradient squared ax=b solver.
c
c    ============================= domn.f ============================
c    domn: preconditioned orthomin sparse iterative ax=b solver.
c    dsdomn: diagonally scaled orthomin sparse iterative ax=b solver.
c    dsluom: incomplete lu orthomin sparse iterative ax=b solver.
c
c    ============================ dgmres.f ===========================
c    dgmres: preconditioned gmres iterative sparse ax=b solver.
c    dsdgmr: diagonally scaled gmres iterative sparse ax=b solver.
c    dslugm: incomplete lu gmres iterative sparse ax=b solver.
c
c    ============================ dmset.f ============================
c       the following routines are used to set up preconditioners.
c
c    dsds: diagonal scaling preconditioner slap set up.
c    dsdscl: diagonally scales/unscales a slap column matrix.
c    dsd2s: diagonal scaling preconditioner slap normal eqns set up.
c    ds2lt: lower triangle preconditioner slap set up.
c    dsics: incomplete cholesky decomp. preconditioner slap set up.
c    dsilus: incomplete lu decomposition preconditioner slap set up.
c
c    ============================ dmvops.f ===========================
c       most of the incomplete  factorization  (ll' and ldu) solvers
c       in this  file require an  intermediate routine  to translate
c       from the slap msolve(n, r, z, nelt, ia,  ja, a, isym, rwork,
c       iwork) calling  convention to the calling  sequence required
c       by  the solve routine.   this generally  is  accomplished by
c       fishing out pointers to the preconditioner (stored in rwork)
c       from the  iwork  array and then making a call to the routine
c       that actually does the backsolve.
c
c    dsmv: slap column format sparse matrix vector product.
c    dsmtv: slap column format sparse matrix (transpose) vector prod.
c    dsdi: diagonal matrix vector multiply.
c    dsli: slap msolve for lower triangle matrix (set up for dsli2).
c    dsli2: lower triangle matrix backsolve.
c    dsllti: slap msolve for ldl' (ic) fact. (set up for dllti2).
c    dllti2: backsolve routine for ldl' factorization.
c    dslui: slap msolve for ldu factorization (set up for dslui2).
c    dslui2: slap backsolve for ldu factorization.
c    dsluti: slap mtsolv for ldu factorization (set up for dslui4).
c    dslui4: slap backsolve for ldu factorization.
c    dsmmti: slap msolve for ldu fact of normal eq (set up for dsmmi2).
c    dsmmi2: slap backsolve for ldu factorization of normal equations.
c
c    =========================== dlaputil.f ==========================
c       the following utility routines are useful additions to slap.
c
c    dbhin: read sparse linear system in the boeing/harwell format.
c    dchkw: slap work/iwork array bounds checker.
c    dcpplt: printer plot of slap column format matrix.
c    ds2y: slap triad to slap column format converter.
c    qs2i1d: quick sort integer array, moving integer and dp arrays.
c            (used by ds2y.)
c    dtin: read in slap triad format linear system.
c    dtout: write out slap triad format linear system.
c
c
c***references  1. mark k. seager, a slap for the masses, in
c                  g. f. carey, ed., parallel supercomputing: methods,
c                  algorithms and applications, wiley, 1989, pp.135-155.
c***routines called  (none)
c***revision history  (yymmdd)
c   890404  date written
c   890404  previous revision date
c   890915  made changes requested at july 1989 cml meeting.  (mks)
c   890921  removed tex from comments.  (fnf)
c   890922  numerous changes to prologue to make closer to slatec
c           standard.  (fnf)
c   890929  numerous changes to reduce sp/dp differences.  (fnf)
c           -----( this produced version 2.0.1. )-----
c   891003  rearranged list of user callable routines to agree with
c           order in source deck.  (fnf)
c   891004  updated reference.
c   910411  prologue converted to version 4.0 format.  (bab)
c           -----( this produced version 2.0.2. )-----
c   910506  minor improvements to prologue.  (fnf)
c   920511  added complete declaration section.  (wrb)
c   920929  corrected format of reference.  (fnf)
c   921019  improved one-line descriptions, reordering some.  (fnf)
c***end prologue  dlpdoc
c***first executable statement  dlpdoc
c
c     this is a *dummy* subroutine and should never be called.
c
      return
c------------- last line of dlpdoc follows -----------------------------
      end
