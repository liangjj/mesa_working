*deck ssd2s
      subroutine ssd2s (n, nelt, ia, ja, a, isym, dinv)
c***begin prologue  ssd2s
c***purpose  diagonal scaling preconditioner slap normal eqns set up.
c            routine to compute the inverse of the diagonal of the
c            matrix a*a', where a is stored in slap-column format.
c***library   slatec (slap)
c***category  d2e
c***type      single precision (ssd2s-s, dsd2s-d)
c***keywords  diagonal, slap sparse
c***author  greenbaum, anne, (courant institute)
c           seager, mark k., (llnl)
c             lawrence livermore national laboratory
c             po box 808, l-60
c             livermore, ca 94550 (510) 423-3141
c             seager@llnl.gov
c***description
c
c *usage:
c     integer n, nelt, ia(nelt), ja(nelt), isym
c     real    a(nelt), dinv(n)
c
c     call ssd2s( n, nelt, ia, ja, a, isym, dinv )
c
c *arguments:
c n      :in       integer
c         order of the matrix.
c nelt   :in       integer.
c         number of elements in arrays ia, ja, and a.
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
c dinv   :out      real dinv(n).
c         upon return this array holds 1./diag(a*a').
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
c       with the slap  format  all  of  the   "inner  loops" of this
c       routine should vectorize  on  machines with hardware support
c       for vector   gather/scatter  operations.  your compiler  may
c       require a compiler directive to  convince it that  there are
c       no  implicit  vector  dependencies.  compiler directives for
c       the alliant    fx/fortran and cri   cft/cft77 compilers  are
c       supplied with the standard slap distribution.
c
c
c *cautions:
c       this routine assumes that the diagonal of a is all  non-zero
c       and that the operation dinv = 1.0/diag(a*a') will not under-
c       flow or overflow. this is done so that the loop  vectorizes.
c       matrices  with zero or near zero or very  large entries will
c       have numerical difficulties  and  must  be fixed before this
c       routine is called.
c
c***see also  ssdcgn
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
c***end prologue  ssd2s
c     .. scalar arguments ..
      integer isym, n, nelt
c     .. array arguments ..
      real a(nelt), dinv(n)
      integer ia(nelt), ja(nelt)
c     .. local scalars ..
      integer i, k, kbgn, kend
c***first executable statement  ssd2s
      do 10 i = 1, n
         dinv(i) = 0
 10   continue
c
c         loop over each column.
cvd$r noconcur
      do 40 i = 1, n
         kbgn = ja(i)
         kend = ja(i+1) - 1
c
c         add in the contributions for each row that has a non-zero
c         in this column.
clll. option assert (nohazard)
cdir$ ivdep
cvd$ nodepchk
         do 20 k = kbgn, kend
            dinv(ia(k)) = dinv(ia(k)) + a(k)**2
 20      continue
         if( isym.eq.1 ) then
c
c         lower triangle stored by columns => upper triangle stored by
c         rows with diagonal being the first entry.  loop across the
c         rest of the row.
            kbgn = kbgn + 1
            if( kbgn.le.kend ) then
               do 30 k = kbgn, kend
                  dinv(i) = dinv(i) + a(k)**2
 30            continue
            endif
         endif
 40   continue
      do 50 i=1,n
         dinv(i) = 1.0e0/dinv(i)
 50   continue
c
      return
c------------- last line of ssd2s follows ----------------------------
      end