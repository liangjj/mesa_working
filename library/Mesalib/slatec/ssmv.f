*deck ssmv
      subroutine ssmv (n, x, y, nelt, ia, ja, a, isym)
c***begin prologue  ssmv
c***purpose  slap column format sparse matrix vector product.
c            routine to calculate the sparse matrix vector product:
c            y = a*x.
c***library   slatec (slap)
c***category  d1b4
c***type      single precision (ssmv-s, dsmv-d)
c***keywords  matrix vector multiply, slap, sparse
c***author  greenbaum, anne, (courant institute)
c           seager, mark k., (llnl)
c             lawrence livermore national laboratory
c             po box 808, l-60
c             livermore, ca 94550 (510) 423-3141
c             seager@llnl.gov
c***description
c
c *usage:
c     integer  n, nelt, ia(nelt), ja(nelt), isym
c     real     x(n), y(n), a(nelt)
c
c     call ssmv(n, x, y, nelt, ia, ja, a, isym )
c
c *arguments:
c n      :in       integer.
c         order of the matrix.
c x      :in       real x(n).
c         the vector that should be multiplied by the matrix.
c y      :out      real y(n).
c         the product of the matrix and the vector.
c nelt   :in       integer.
c         number of non-zeros stored in a.
c ia     :in       integer ia(nelt).
c ja     :in       integer ja(nelt).
c a      :in       real a(nelt).
c         these arrays should hold the matrix a in the slap column
c         format.  see "description", below.
c isym   :in       integer.
c         flag to indicate symmetric storage format.
c         if isym=0, all non-zero entries of the matrix are stored.
c         if isym=1, the matrix is symmetric, and only the upper
c         or lower triangle of the matrix is stored.
c
c *description
c       =================== s l a p column format ==================
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
c       with  the slap  format  the "inner  loops" of  this  routine
c       should vectorize   on machines with   hardware  support  for
c       vector gather/scatter operations.  your compiler may require
c       a  compiler directive  to  convince   it that there  are  no
c       implicit vector  dependencies.  compiler directives  for the
c       alliant fx/fortran and cri cft/cft77 compilers  are supplied
c       with the standard slap distribution.
c
c *cautions:
c     this   routine   assumes  that  the matrix a is stored in slap
c     column format.  it does not check  for  this (for  speed)  and
c     evil, ugly, ornery and nasty things  will happen if the matrix
c     data  structure  is,  in fact, not slap column.  beware of the
c     wrong data structure!!!
c
c***see also  ssmtv
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
c   930701  updated category section.  (fnf, wrb)
c***end prologue  ssmv
c     .. scalar arguments ..
      integer isym, n, nelt
c     .. array arguments ..
      real a(nelt), x(n), y(n)
      integer ia(nelt), ja(nelt)
c     .. local scalars ..
      integer i, ibgn, icol, iend, irow, j, jbgn, jend
c***first executable statement  ssmv
c
c         zero out the result vector.
c
      do 10 i = 1, n
         y(i) = 0
 10   continue
c
c         multiply by a.
c
cvd$r noconcur
      do 30 icol = 1, n
         ibgn = ja(icol)
         iend = ja(icol+1)-1
clll. option assert (nohazard)
cdir$ ivdep
cvd$ nodepchk
         do 20 i = ibgn, iend
            y(ia(i)) = y(ia(i)) + a(i)*x(icol)
 20      continue
 30   continue
c
      if( isym.eq.1 ) then
c
c         the matrix is non-symmetric.  need to get the other half in...
c         this loops assumes that the diagonal is the first entry in
c         each column.
c
         do 50 irow = 1, n
            jbgn = ja(irow)+1
            jend = ja(irow+1)-1
            if( jbgn.gt.jend ) goto 50
            do 40 j = jbgn, jend
               y(irow) = y(irow) + a(j)*x(ia(j))
 40         continue
 50      continue
      endif
      return
c------------- last line of ssmv follows ----------------------------
      end
