*deck sslui2
      subroutine sslui2 (n, b, x, il, jl, l, dinv, iu, ju, u)
c***begin prologue  sslui2
c***purpose  slap backsolve for ldu factorization.
c            routine to solve a system of the form  l*d*u x = b,
c            where l is a unit lower triangular matrix, d is a diagonal
c            matrix, and u is a unit upper triangular matrix.
c***library   slatec (slap)
c***category  d2e
c***type      single precision (sslui2-s, dslui2-d)
c***keywords  iterative precondition, non-symmetric linear system solve,
c             slap, sparse
c***author  greenbaum, anne, (courant institute)
c           seager, mark k., (llnl)
c             lawrence livermore national laboratory
c             po box 808, l-60
c             livermore, ca 94550 (510) 423-3141
c             seager@llnl.gov
c***description
c
c *usage:
c     integer n, il(nl), jl(nl), iu(nu), ju(nu)
c     real    b(n), x(n), l(nl), dinv(n), u(nu)
c
c     call sslui2( n, b, x, il, jl, l, dinv, iu, ju, u )
c
c *arguments:
c n      :in       integer
c         order of the matrix.
c b      :in       real b(n).
c         right hand side.
c x      :out      real x(n).
c         solution of l*d*u x = b.
c il     :in       integer il(nl).
c jl     :in       integer jl(nl).
c l      :in       real     l(nl).
c         il, jl, l contain the unit  lower triangular factor of the
c         incomplete decomposition of some matrix stored in slap row
c         format.  the diagonal of ones *is* stored.  this structure
c         can   be   set  up  by   the  ssilus  routine.   see   the
c         "description", below  for more   details about   the  slap
c         format.  (nl is the number of non-zeros in the l array.)
c dinv   :in       real dinv(n).
c         inverse of the diagonal matrix d.
c iu     :in       integer iu(nu).
c ju     :in       integer ju(nu).
c u      :in       real     u(nu).
c         iu, ju, u contain the unit upper triangular factor  of the
c         incomplete decomposition  of  some  matrix stored in  slap
c         column format.   the diagonal of ones  *is* stored.   this
c         structure can be set up  by the ssilus routine.  see   the
c         "description", below   for  more   details about  the slap
c         format.  (nu is the number of non-zeros in the u array.)
c
c *description:
c       this routine is supplied with  the slap package as a routine
c       to  perform  the  msolve operation  in   the  sir and   sbcg
c       iteration routines for  the  drivers ssilur and sslubc.   it
c       must  be called  via   the  slap  msolve  calling   sequence
c       convention interface routine sslui.
c         **** this routine itself does not conform to the ****
c               **** slap msolve calling convention ****
c
c       il, jl, l should contain the unit lower triangular factor of
c       the incomplete decomposition of the a matrix  stored in slap
c       row format.  iu, ju, u should contain  the unit upper factor
c       of the  incomplete decomposition of  the a matrix  stored in
c       slap column format this ilu factorization can be computed by
c       the ssilus routine. the diagonals (which are all one's) are
c       stored.
c
c       =================== s l a p column format ==================
c
c       this routine  requires that  the matrix a  be stored in  the
c       slap column format.  in this format the non-zeros are stored
c       counting down columns (except for  the diagonal entry, which
c       must appear first in each  "column")  and are stored  in the
c       real array a.  in other words, for each column in the matrix
c       put the diagonal entry in a.  then put in the other non-zero
c       elements going down   the  column (except  the diagonal)  in
c       order.  the ia array holds the row  index for each non-zero.
c       the ja array holds the offsets into the ia, a arrays for the
c       beginning of   each    column.    that  is,    ia(ja(icol)),
c       a(ja(icol)) points to the beginning of the icol-th column in
c       ia and  a.  ia(ja(icol+1)-1),  a(ja(icol+1)-1) points to the
c       end  of   the icol-th  column.  note   that  we  always have
c       ja(n+1) = nelt+1, where  n  is the number of columns in  the
c       matrix and  nelt   is the number of non-zeros in the matrix.
c
c       here is an example of the  slap column  storage format for a
c       5x5 matrix (in the a and ia arrays '|'  denotes the end of a
c       column):
c
c           5x5 matrix      slap column format for 5x5 matrix on left.
c                              1  2  3    4  5    6  7    8    9 10 11
c       |11 12  0  0 15|   a: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35
c       |21 22  0  0  0|  ia:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
c       | 0  0 33  0 35|  ja:  1  4  6    8  9   12
c       | 0  0  0 44  0|
c       |51  0 53  0 55|
c
c       ==================== s l a p row format ====================
c
c       this routine requires  that the matrix a  be  stored  in the
c       slap  row format.   in this format  the non-zeros are stored
c       counting across  rows (except for the diagonal  entry, which
c       must appear first in each "row") and  are stored in the real
c       array a.  in other words, for each row in the matrix put the
c       diagonal entry in  a.   then   put  in the   other  non-zero
c       elements   going  across the  row (except   the diagonal) in
c       order.   the  ja array  holds   the column   index for  each
c       non-zero.   the ia  array holds the  offsets into  the ja, a
c       arrays  for   the   beginning  of   each  row.   that    is,
c       ja(ia(irow)),  a(ia(irow)) points  to  the beginning  of the
c       irow-th row in ja and a.   ja(ia(irow+1)-1), a(ia(irow+1)-1)
c       points to the  end of the  irow-th row.  note that we always
c       have ia(n+1) =  nelt+1, where  n  is  the number of rows  in
c       the matrix  and nelt  is the  number   of  non-zeros in  the
c       matrix.
c
c       here is an example of the slap row storage format for a  5x5
c       matrix (in the a and ja arrays '|' denotes the end of a row):
c
c           5x5 matrix         slap row format for 5x5 matrix on left.
c                              1  2  3    4  5    6  7    8    9 10 11
c       |11 12  0  0 15|   a: 11 12 15 | 22 21 | 33 35 | 44 | 55 51 53
c       |21 22  0  0  0|  ja:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
c       | 0  0 33  0 35|  ia:  1  4  6    8  9   12
c       | 0  0  0 44  0|
c       |51  0 53  0 55|
c
c       with  the slap  format  the "inner  loops" of  this  routine
c       should vectorize   on machines with   hardware  support  for
c       vector gather/scatter operations.  your compiler may require
c       a  compiler directive  to  convince   it that there  are  no
c       implicit vector  dependencies.  compiler directives  for the
c       alliant fx/fortran and cri cft/cft77 compilers  are supplied
c       with the standard slap distribution.
c
c***see also  ssilus
c***references  (none)
c***routines called  (none)
c***revision history  (yymmdd)
c   871119  date written
c   881213  previous revision date
c   890915  made changes requested at july 1989 cml meeting.  (mks)
c   890922  numerous changes to prologue to make closer to slatec
c           standard.  (fnf)
c   890929  numerous changes to reduce sp/dp differences.  (fnf)
c   910411  prologue converted to version 4.0 format.  (bab)
c   920511  added complete declaration section.  (wrb)
c   921113  corrected c***category line.  (fnf)
c   930701  updated category section.  (fnf, wrb)
c***end prologue  sslui2
c     .. scalar arguments ..
      integer n
c     .. array arguments ..
      real b(n), dinv(n), l(*), u(*), x(n)
      integer il(*), iu(*), jl(*), ju(*)
c     .. local scalars ..
      integer i, icol, irow, j, jbgn, jend
c***first executable statement  sslui2
c
c         solve  l*y = b,  storing result in x, l stored by rows.
c
      do 10 i = 1, n
         x(i) = b(i)
 10   continue
      do 30 irow = 2, n
         jbgn = il(irow)
         jend = il(irow+1)-1
         if( jbgn.le.jend ) then
clll. option assert (nohazard)
cdir$ ivdep
cvd$ assoc
cvd$ nodepchk
            do 20 j = jbgn, jend
               x(irow) = x(irow) - l(j)*x(jl(j))
 20         continue
         endif
 30   continue
c
c         solve  d*z = y,  storing result in x.
      do 40 i=1,n
         x(i) = x(i)*dinv(i)
 40   continue
c
c         solve  u*x = z, u stored by columns.
      do 60 icol = n, 2, -1
         jbgn = ju(icol)
         jend = ju(icol+1)-1
         if( jbgn.le.jend ) then
clll. option assert (nohazard)
cdir$ ivdep
cvd$ nodepchk
            do 50 j = jbgn, jend
               x(iu(j)) = x(iu(j)) - u(j)*x(icol)
 50         continue
         endif
 60   continue
c
      return
c------------- last line of sslui2 follows ----------------------------
      end
