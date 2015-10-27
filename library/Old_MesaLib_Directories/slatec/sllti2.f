*deck sllti2
      subroutine sllti2 (n, b, x, nel, iel, jel, el, dinv)
c***begin prologue  sllti2
c***purpose  slap backsolve routine for ldl' factorization.
c            routine to solve a system of the form  l*d*l' x = b,
c            where l is a unit lower triangular matrix and d is a
c            diagonal matrix and ' means transpose.
c***library   slatec (slap)
c***category  d2e
c***type      single precision (sllti2-s, dllti2-d)
c***keywords  incomplete factorization, iterative precondition, slap,
c             sparse, symmetric linear system solve
c***author  greenbaum, anne, (courant institute)
c           seager, mark k., (llnl)
c             lawrence livermore national laboratory
c             po box 808, l-60
c             livermore, ca 94550 (510) 423-3141
c             seager@llnl.gov
c***description
c
c *usage:
c     integer n, nel, iel(nel), jel(nel)
c     real    b(n), x(n), el(nel), dinv(n)
c
c     call sllti2( n, b, x, nel, iel, jel, el, dinv )
c
c *arguments:
c n      :in       integer
c         order of the matrix.
c b      :in       real b(n).
c         right hand side vector.
c x      :out      real x(n).
c         solution to l*d*l' x = b.
c nel    :in       integer.
c         number of non-zeros in the el array.
c iel    :in       integer iel(nel).
c jel    :in       integer jel(nel).
c el     :in       real     el(nel).
c         iel, jel, el contain the unit lower triangular factor   of
c         the incomplete decomposition   of the a  matrix  stored in
c         slap row format.   the diagonal of ones *is* stored.  this
c         structure can be set  up  by  the ss2lt routine.  see  the
c         "description", below for more details about the  slap  row
c         format.
c dinv   :in       real dinv(n).
c         inverse of the diagonal matrix d.
c
c *description:
c       this routine is supplied with  the slap package as a routine
c       to perform the msolve operation in the scg iteration routine
c       for  the driver  routine ssiccg.   it must be called via the
c       slap  msolve calling sequence  convention  interface routine
c       sslli.
c         **** this routine itself does not conform to the ****
c               **** slap msolve calling convention ****
c
c       iel, jel, el should contain the unit lower triangular factor
c       of  the incomplete decomposition of  the a matrix  stored in
c       slap row format.   this ic factorization  can be computed by
c       the  ssics routine.  the  diagonal  (which is all one's) is
c       stored.
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
c       with  the slap  row format  the "inner loop" of this routine
c       should vectorize   on machines with   hardware  support  for
c       vector gather/scatter operations.  your compiler may require
c       a  compiler directive  to  convince   it that there  are  no
c       implicit vector  dependencies.  compiler directives  for the
c       alliant fx/fortran and cri cft/cft77 compilers  are supplied
c       with the standard slap distribution.
c
c***see also  ssiccg, ssics
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
c***end prologue  sllti2
c     .. scalar arguments ..
      integer n, nel
c     .. array arguments ..
      real b(n), dinv(n), el(nel), x(n)
      integer iel(nel), jel(nel)
c     .. local scalars ..
      integer i, ibgn, iend, irow
c***first executable statement  sllti2
c
c         solve  l*y = b,  storing result in x.
c
      do 10 i=1,n
         x(i) = b(i)
 10   continue
      do 30 irow = 1, n
         ibgn = iel(irow) + 1
         iend = iel(irow+1) - 1
         if( ibgn.le.iend ) then
clll. option assert (nohazard)
cdir$ ivdep
cvd$ noconcur
cvd$ nodepchk
            do 20 i = ibgn, iend
               x(irow) = x(irow) - el(i)*x(jel(i))
 20         continue
         endif
 30   continue
c
c         solve  d*z = y,  storing result in x.
c
      do 40 i=1,n
         x(i) = x(i)*dinv(i)
 40   continue
c
c         solve  l-trans*x = z.
c
      do 60 irow = n, 2, -1
         ibgn = iel(irow) + 1
         iend = iel(irow+1) - 1
         if( ibgn.le.iend ) then
clll. option assert (nohazard)
cdir$ ivdep
cvd$ noconcur
cvd$ nodepchk
            do 50 i = ibgn, iend
               x(jel(i)) = x(jel(i)) - el(i)*x(irow)
 50         continue
         endif
 60   continue
c
      return
c------------- last line of sllti2 follows ----------------------------
      end
