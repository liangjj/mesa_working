*deck ssli2
      subroutine ssli2 (n, b, x, nel, iel, jel, el)
c***begin prologue  ssli2
c***purpose  slap lower triangle matrix backsolve.
c            routine to solve a system of the form  lx = b , where l
c            is a lower triangular matrix.
c***library   slatec (slap)
c***category  d2a3
c***type      single precision (ssli2-s, dsli2-d)
c***keywords  iterative precondition, linear system solve, slap, sparse
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
c     real    b(n), x(n), el(nel)
c
c     call ssli2( n, b, x, nel, iel, jel, el )
c
c *arguments:
c n      :in       integer
c         order of the matrix.
c b      :in       real b(n).
c         right hand side vector.
c x      :out      real x(n).
c         solution to lx = b.
c nel    :in       integer.
c         number of non-zeros in the el array.
c iel    :in       integer iel(nel).
c jel    :in       integer jel(nel).
c el     :in       real     el(nel).
c         iel, jel, el contain the unit lower triangular factor   of
c         the incomplete decomposition   of the a  matrix  stored in
c         slap row format.  the diagonal of  ones *is* stored.  this
c         structure can be  set up by the  ss2lt  routine.  see  the
c         "description", below, for more details about the  slap row
c         format.
c
c *description:
c       this routine is supplied with the slap package  as a routine
c       to perform the msolve operation in the sir iteration routine
c       for the driver routine ssgs.  it must be called via the slap
c       msolve calling sequence convention interface routine ssli.
c         **** this routine itself does not conform to the ****
c               **** slap msolve calling convention ****
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
c***see also  ssli
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
c***end prologue  ssli2
c     .. scalar arguments ..
      integer n, nel
c     .. array arguments ..
      real b(n), el(nel), x(n)
      integer iel(nel), jel(nel)
c     .. local scalars ..
      integer i, icol, j, jbgn, jend
c***first executable statement  ssli2
c
c         initialize the solution by copying the right hands side
c         into it.
c
      do 10 i=1,n
         x(i) = b(i)
 10   continue
c
cvd$ noconcur
      do 30 icol = 1, n
         x(icol) = x(icol)/el(jel(icol))
         jbgn = jel(icol) + 1
         jend = jel(icol+1) - 1
         if( jbgn.le.jend ) then
clll. option assert (nohazard)
cdir$ ivdep
cvd$ noconcur
cvd$ nodepchk
            do 20 j = jbgn, jend
               x(iel(j)) = x(iel(j)) - el(j)*x(icol)
 20         continue
         endif
 30   continue
c
      return
c------------- last line of ssli2 follows ----------------------------
      end
