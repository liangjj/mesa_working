*deck dsics
      subroutine dsics (n, nelt, ia, ja, a, isym, nel, iel, jel, el, d,
     +   r, iwarn)
c***begin prologue  dsics
c***purpose  incompl. cholesky decomposition preconditioner slap set up.
c            routine to generate the incomplete cholesky decomposition,
c            l*d*l-trans, of a symmetric positive definite matrix, a,
c            which is stored in slap column format.  the unit lower
c            triangular matrix l is stored by rows, and the inverse of
c            the diagonal matrix d is stored.
c***library   slatec (slap)
c***category  d2e
c***type      double precision (ssics-s, dsics-d)
c***keywords  incomplete cholesky factorization,
c             iterative precondition, linear system, slap sparse
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
c     integer nel, iel(nel), jel(nel), iwarn
c     double precision a(nelt), el(nel), d(n), r(n)
c
c     call dsics( n, nelt, ia, ja, a, isym, nel, iel, jel, el, d, r,
c    $    iwarn )
c
c *arguments:
c n      :in       integer.
c         order of the matrix.
c nelt   :in       integer.
c         number of elements in arrays ia, ja, and a.
c ia     :inout    integer ia(nelt).
c ja     :inout    integer ja(nelt).
c a      :inout    double precision a(nelt).
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
c el     :out      double precision el(nel).
c         iel, jel, el contain the unit lower triangular factor  of the
c         incomplete decomposition   of the a  matrix  stored  in  slap
c         row format.   the diagonal of   ones   *is*   stored.     see
c         "description", below for more details about the slap row fmt.
c d      :out      double precision d(n)
c         upon return this array holds d(i) = 1./diag(a).
c r      :work     double precision r(n).
c         temporary double precision workspace needed for the
c         factorization.
c iwarn  :out      integer.
c         this is a warning variable and is zero if the ic factoriza-
c         tion goes well.  it is set to the row index corresponding to
c         the last zero pivot found.  see "description", below.
c
c *description
c       =================== s l a p column format ==================
c       this routine  requires that  the matrix a  be stored in  the
c       slap column format.  in this format the non-zeros are stored
c       counting down columns (except for  the diagonal entry, which
c       must appear first in each  "column")  and are stored  in the
c       double precision array a.   in other words,  for each column
c       in the matrix put the diagonal entry in  a.  then put in the
c       other non-zero  elements going down  the column (except  the
c       diagonal) in order.   the  ia array holds the  row index for
c       each non-zero.  the ja array holds the offsets  into the ia,
c       a arrays  for  the  beginning  of each   column.   that  is,
c       ia(ja(icol)),  a(ja(icol)) points   to the beginning  of the
c       icol-th   column    in    ia and   a.      ia(ja(icol+1)-1),
c       a(ja(icol+1)-1) points to  the  end of the   icol-th column.
c       note that we always have  ja(n+1) = nelt+1,  where n is  the
c       number of columns in  the matrix and nelt  is the number  of
c       non-zeros in the matrix.
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
c       must  appear first  in each  "row")  and  are stored  in the
c       double precision  array a.  in other words, for each row  in
c       the matrix  put the diagonal  entry in a.   then put in  the
c       other  non-zero elements  going across  the row  (except the
c       diagonal) in order.  the ja array holds the column index for
c       each non-zero.  the ia array holds the offsets  into the ja,
c       a  arrays  for  the   beginning  of  each  row.    that  is,
c       ja(ia(irow)),a(ia(irow)) are the first elements of the irow-
c       th row in  ja and a,  and  ja(ia(irow+1)-1), a(ia(irow+1)-1)
c       are  the last elements  of the  irow-th row.   note  that we
c       always have  ia(n+1) = nelt+1, where n is the number of rows
c       in the matrix  and  nelt is the  number of non-zeros  in the
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
c       with the slap  format some  of  the   "inner  loops" of this
c       routine should vectorize  on  machines with hardware support
c       for vector   gather/scatter  operations.  your compiler  may
c       require a compiler directive to  convince it that  there are
c       no  implicit  vector  dependencies.  compiler directives for
c       the alliant    fx/fortran and cri   cft/cft77 compilers  are
c       supplied with the standard slap distribution.
c
c       the ic factorization does not always exist for spd matrices.
c       in the event that a zero pivot is found it is set  to be 1.0
c       and the factorization proceeds.   the integer variable iwarn
c       is set to the last row where the diagonal was fudged.  this
c       eventuality hardly ever occurs in practice.
c
c***see also  dcg, dsiccg
c***references  1. gene golub and charles van loan, matrix computations,
c                  johns hopkins university press, baltimore, maryland,
c                  1983.
c***routines called  xermsg
c***revision history  (yymmdd)
c   890404  date written
c   890404  previous revision date
c   890915  made changes requested at july 1989 cml meeting.  (mks)
c   890922  numerous changes to prologue to make closer to slatec
c           standard.  (fnf)
c   890929  numerous changes to reduce sp/dp differences.  (fnf)
c   900805  changed xerrwv calls to calls to xermsg.  (rwc)
c   910411  prologue converted to version 4.0 format.  (bab)
c   920511  added complete declaration section.  (wrb)
c   920929  corrected format of reference.  (fnf)
c   930701  updated category section.  (fnf, wrb)
c***end prologue  dsics
c     .. scalar arguments ..
      integer isym, iwarn, n, nel, nelt
c     .. array arguments ..
      double precision a(nelt), d(n), el(nel), r(n)
      integer ia(nelt), iel(nel), ja(nelt), jel(nel)
c     .. local scalars ..
      double precision eltmp
      integer i, ibgn, ic, icbgn, icend, icol, iend, ir, irbgn, irend,
     +        irow, irr, j, jbgn, jeltmp, jend
      character xern1*8
c     .. external subroutines ..
      external xermsg
c***first executable statement  dsics
c
c         set the lower triangle in iel, jel, el
c
      iwarn = 0
c
c         all matrix elements stored in ia, ja, a.  pick out the lower
c         triangle (making sure that the diagonal of el is one) and
c         store by rows.
c
      nel = 1
      iel(1) = 1
      jel(1) = 1
      el(1) = 1
      d(1) = a(1)
cvd$r noconcur
      do 30 irow = 2, n
c         put in the diagonal.
         nel = nel + 1
         iel(irow) = nel
         jel(nel) = irow
         el(nel) = 1
         d(irow) = a(ja(irow))
c
c         look in all the lower triangle columns for a matching row.
c         since the matrix is symmetric, we can look across the
c         irow-th row by looking down the irow-th column (if it is
c         stored isym=0)...
         if( isym.eq.0 ) then
            icbgn = ja(irow)
            icend = ja(irow+1)-1
         else
            icbgn = 1
            icend = irow-1
         endif
         do 20 ic = icbgn, icend
            if( isym.eq.0 ) then
               icol = ia(ic)
               if( icol.ge.irow ) goto 20
            else
               icol = ic
            endif
            jbgn = ja(icol)+1
            jend = ja(icol+1)-1
            if( jbgn.le.jend .and. ia(jend).ge.irow ) then
cvd$ novector
               do 10 j = jbgn, jend
                  if( ia(j).eq.irow ) then
                     nel = nel + 1
                     jel(nel) = icol
                     el(nel)  = a(j)
                     goto 20
                  endif
 10            continue
            endif
 20      continue
 30   continue
      iel(n+1) = nel+1
c
c         sort rows of lower triangle into descending order (count out
c         along rows out from diagonal).
c
      do 60 irow = 2, n
         ibgn = iel(irow)+1
         iend = iel(irow+1)-1
         if( ibgn.lt.iend ) then
            do 50 i = ibgn, iend-1
cvd$ novector
               do 40 j = i+1, iend
                  if( jel(i).gt.jel(j) ) then
                     jeltmp = jel(j)
                     jel(j) = jel(i)
                     jel(i) = jeltmp
                     eltmp = el(j)
                     el(j) = el(i)
                     el(i) = eltmp
                  endif
 40            continue
 50         continue
         endif
 60   continue
c
c         perform the incomplete cholesky decomposition by looping
c         over the rows.
c         scale the first column.  use the structure of a to pick out
c         the rows with something in column 1.
c
      irbgn = ja(1)+1
      irend = ja(2)-1
      do 65 irr = irbgn, irend
         ir = ia(irr)
c         find the index into el for el(1,ir).
c         hint: it's the second entry.
         i = iel(ir)+1
         el(i) = el(i)/d(1)
 65   continue
c
      do 110 irow = 2, n
c
c         update the irow-th diagonal.
c
         do 66 i = 1, irow-1
            r(i) = 0
 66      continue
         ibgn = iel(irow)+1
         iend = iel(irow+1)-1
         if( ibgn.le.iend ) then
clll. option assert (nohazard)
cdir$ ivdep
cvd$ nodepchk
            do 70 i = ibgn, iend
               r(jel(i)) = el(i)*d(jel(i))
               d(irow) = d(irow) - el(i)*r(jel(i))
 70         continue
c
c         check to see if we have a problem with the diagonal.
c
            if( d(irow).le.0.0d0 ) then
               if( iwarn.eq.0 ) iwarn = irow
               d(irow) = 1
            endif
         endif
c
c         update each el(irow+1:n,irow), if there are any.
c         use the structure of a to determine the non-zero elements
c         of the irow-th column of el.
c
         irbgn = ja(irow)
         irend = ja(irow+1)-1
         do 100 irr = irbgn, irend
            ir = ia(irr)
            if( ir.le.irow ) goto 100
c         find the index into el for el(ir,irow)
            ibgn = iel(ir)+1
            iend = iel(ir+1)-1
            if( jel(ibgn).gt.irow ) goto 100
            do 90 i = ibgn, iend
               if( jel(i).eq.irow ) then
                  icend = iend
 91               if( jel(icend).ge.irow ) then
                     icend = icend - 1
                     goto 91
                  endif
c         sum up the el(ir,1:irow-1)*r(1:irow-1) contributions.
clll. option assert (nohazard)
cdir$ ivdep
cvd$ nodepchk
                  do 80 ic = ibgn, icend
                     el(i) = el(i) - el(ic)*r(jel(ic))
 80               continue
                  el(i) = el(i)/d(irow)
                  goto 100
               endif
 90         continue
c
c         if we get here, we have real problems...
            write (xern1, '(i8)') irow
            call xermsg ('slatec', 'dsics',
     $         'a and el data structure mismatch in row '// xern1, 1, 2)
 100     continue
 110  continue
c
c         replace diagonals by their inverses.
c
cvd$ concur
      do 120 i =1, n
         d(i) = 1.0d0/d(i)
 120  continue
      return
c------------- last line of dsics follows ----------------------------
      end
