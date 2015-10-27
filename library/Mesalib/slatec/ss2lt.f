*deck ss2lt
      subroutine ss2lt (n, nelt, ia, ja, a, isym, nel, iel, jel, el)
c***begin prologue  ss2lt
c***purpose  lower triangle preconditioner slap set up.
c            routine to store the lower triangle of a matrix stored
c            in the slap column format.
c***library   slatec (slap)
c***category  d2e
c***type      single precision (ss2lt-s, ds2lt-d)
c***keywords  linear system, lower triangle, slap sparse
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
c     integer nel, iel(nel), jel(nel)
c     real    a(nelt), el(nel)
c
c     call ss2lt( n, nelt, ia, ja, a, isym, nel, iel, jel, el )
c
c *arguments:
c n      :in       integer
c         order of the matrix.
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
c         if isym=1, the matrix is symmetric, and only the lower
c         triangle of the matrix is stored.
c nel    :out      integer.
c         number of non-zeros in the lower triangle of a.   also
c         corresponds to the length of the iel, jel, el arrays.
c iel    :out      integer iel(nel).
c jel    :out      integer jel(nel).
c el     :out      real     el(nel).
c         iel, jel, el contain the lower triangle of the a matrix
c         stored in slap column format.  see "description", below,
c         for more details bout the slap column format.
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
c***end prologue  ss2lt
c     .. scalar arguments ..
      integer isym, n, nel, nelt
c     .. array arguments ..
      real a(nelt), el(nelt)
      integer ia(nelt), iel(nel), ja(nelt), jel(nel)
c     .. local scalars ..
      integer i, icol, j, jbgn, jend
c***first executable statement  ss2lt
      if( isym.eq.0 ) then
c
c         the matrix is stored non-symmetricly.  pick out the lower
c         triangle.
c
         nel = 0
         do 20 icol = 1, n
            jel(icol) = nel+1
            jbgn = ja(icol)
            jend = ja(icol+1)-1
cvd$ novector
            do 10 j = jbgn, jend
               if( ia(j).ge.icol ) then
                  nel = nel + 1
                  iel(nel) = ia(j)
                  el(nel)  = a(j)
               endif
 10         continue
 20      continue
         jel(n+1) = nel+1
      else
c
c         the matrix is symmetric and only the lower triangle is
c         stored.  copy it to iel, jel, el.
c
         nel = nelt
         do 30 i = 1, nelt
            iel(i) = ia(i)
            el(i) = a(i)
 30      continue
         do 40 i = 1, n+1
            jel(i) = ja(i)
 40      continue
      endif
      return
c------------- last line of ss2lt follows ----------------------------
      end
