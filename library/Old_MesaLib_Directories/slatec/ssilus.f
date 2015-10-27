*deck ssilus
      subroutine ssilus (n, nelt, ia, ja, a, isym, nl, il, jl, l, dinv,
     +   nu, iu, ju, u, nrow, ncol)
c***begin prologue  ssilus
c***purpose  incomplete lu decomposition preconditioner slap set up.
c            routine to generate the incomplete ldu decomposition of a
c            matrix.  the unit lower triangular factor l is stored by
c            rows and the unit upper triangular factor u is stored by
c            columns.  the inverse of the diagonal matrix d is stored.
c            no fill in is allowed.
c***library   slatec (slap)
c***category  d2e
c***type      single precision (ssilus-s, dsilus-d)
c***keywords  incomplete lu factorization, iterative precondition,
c             non-symmetric linear system, slap, sparse
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
c     integer nl, il(nl), jl(nl), nu, iu(nu), ju(nu)
c     integer nrow(n), ncol(n)
c     real    a(nelt), l(nl), dinv(n), u(nu)
c
c     call ssilus( n, nelt, ia, ja, a, isym, nl, il, jl, l,
c    $    dinv, nu, iu, ju, u, nrow, ncol )
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
c         if isym=1, the matrix is symmetric, and only the lower
c         triangle of the matrix is stored.
c nl     :out      integer.
c         number of non-zeros in the l array.
c il     :out      integer il(nl).
c jl     :out      integer jl(nl).
c l      :out      real     l(nl).
c         il, jl, l  contain the unit lower triangular factor of  the
c         incomplete decomposition  of some  matrix stored  in   slap
c         row format.     the   diagonal  of ones  *is*  stored.  see
c         "description", below for more details about the slap format.
c nu     :out      integer.
c         number of non-zeros in the u array.
c iu     :out      integer iu(nu).
c ju     :out      integer ju(nu).
c u      :out      real     u(nu).
c         iu, ju, u contain   the unit upper triangular factor of the
c         incomplete  decomposition    of some matrix  stored in slap
c         column  format.   the diagonal of ones   *is*  stored.  see
c         "description", below  for  more  details  about  the   slap
c         format.
c nrow   :work     integer nrow(n).
c         nrow(i) is the number of non-zero elements in the i-th row
c         of l.
c ncol   :work     integer ncol(n).
c         ncol(i) is the number of non-zero elements in the i-th
c         column of u.
c
c *description
c       il, jl, l should contain the unit  lower triangular factor of
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
c***see also  silur
c***references  1. gene golub and charles van loan, matrix computations,
c                  johns hopkins university press, baltimore, maryland,
c                  1983.
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
c   920929  corrected format of reference.  (fnf)
c   930701  updated category section.  (fnf, wrb)
c***end prologue  ssilus
c     .. scalar arguments ..
      integer isym, n, nelt, nl, nu
c     .. array arguments ..
      real a(nelt), dinv(n), l(nl), u(nu)
      integer ia(nelt), il(nl), iu(nu), ja(nelt), jl(nl), ju(nu),
     +        ncol(n), nrow(n)
c     .. local scalars ..
      real temp
      integer i, ibgn, icol, iend, indx, indx1, indx2, indxc1, indxc2,
     +        indxr1, indxr2, irow, itemp, j, jbgn, jend, jtemp, k, kc,
     +        kr
c***first executable statement  ssilus
c
c         count number of elements in each row of the lower triangle.
c
      do 10 i=1,n
         nrow(i) = 0
         ncol(i) = 0
 10   continue
cvd$r noconcur
cvd$r novector
      do 30 icol = 1, n
         jbgn = ja(icol)+1
         jend = ja(icol+1)-1
         if( jbgn.le.jend ) then
            do 20 j = jbgn, jend
               if( ia(j).lt.icol ) then
                  ncol(icol) = ncol(icol) + 1
               else
                  nrow(ia(j)) = nrow(ia(j)) + 1
                  if( isym.ne.0 ) ncol(ia(j)) = ncol(ia(j)) + 1
               endif
 20         continue
         endif
 30   continue
      ju(1) = 1
      il(1) = 1
      do 40 icol = 1, n
         il(icol+1) = il(icol) + nrow(icol)
         ju(icol+1) = ju(icol) + ncol(icol)
         nrow(icol) = il(icol)
         ncol(icol) = ju(icol)
 40   continue
c
c         copy the matrix a into the l and u structures.
      do 60 icol = 1, n
         dinv(icol) = a(ja(icol))
         jbgn = ja(icol)+1
         jend = ja(icol+1)-1
         if( jbgn.le.jend ) then
            do 50 j = jbgn, jend
               irow = ia(j)
               if( irow.lt.icol ) then
c         part of the upper triangle.
                  iu(ncol(icol)) = irow
                  u(ncol(icol)) = a(j)
                  ncol(icol) = ncol(icol) + 1
               else
c         part of the lower triangle (stored by row).
                  jl(nrow(irow)) = icol
                  l(nrow(irow)) = a(j)
                  nrow(irow) = nrow(irow) + 1
                  if( isym.ne.0 ) then
c         symmetric...copy lower triangle into upper triangle as well.
                     iu(ncol(irow)) = icol
                     u(ncol(irow)) = a(j)
                     ncol(irow) = ncol(irow) + 1
                  endif
               endif
 50         continue
         endif
 60   continue
c
c         sort the rows of l and the columns of u.
      do 110 k = 2, n
         jbgn = ju(k)
         jend = ju(k+1)-1
         if( jbgn.lt.jend ) then
            do 80 j = jbgn, jend-1
               do 70 i = j+1, jend
                  if( iu(j).gt.iu(i) ) then
                     itemp = iu(j)
                     iu(j) = iu(i)
                     iu(i) = itemp
                     temp = u(j)
                     u(j) = u(i)
                     u(i) = temp
                  endif
 70            continue
 80         continue
         endif
         ibgn = il(k)
         iend = il(k+1)-1
         if( ibgn.lt.iend ) then
            do 100 i = ibgn, iend-1
               do 90 j = i+1, iend
                  if( jl(i).gt.jl(j) ) then
                     jtemp = ju(i)
                     ju(i) = ju(j)
                     ju(j) = jtemp
                     temp = l(i)
                     l(i) = l(j)
                     l(j) = temp
                  endif
 90            continue
 100        continue
         endif
 110  continue
c
c         perform the incomplete ldu decomposition.
      do 300 i=2,n
c
c           i-th row of l
         indx1 = il(i)
         indx2 = il(i+1) - 1
         if(indx1 .gt. indx2) go to 200
         do 190 indx=indx1,indx2
            if(indx .eq. indx1) go to 180
            indxr1 = indx1
            indxr2 = indx - 1
            indxc1 = ju(jl(indx))
            indxc2 = ju(jl(indx)+1) - 1
            if(indxc1 .gt. indxc2) go to 180
 160        kr = jl(indxr1)
 170        kc = iu(indxc1)
            if(kr .gt. kc) then
               indxc1 = indxc1 + 1
               if(indxc1 .le. indxc2) go to 170
            elseif(kr .lt. kc) then
               indxr1 = indxr1 + 1
               if(indxr1 .le. indxr2) go to 160
            elseif(kr .eq. kc) then
               l(indx) = l(indx) - l(indxr1)*dinv(kc)*u(indxc1)
               indxr1 = indxr1 + 1
               indxc1 = indxc1 + 1
               if(indxr1 .le. indxr2 .and. indxc1 .le. indxc2) go to 160
            endif
 180        l(indx) = l(indx)/dinv(jl(indx))
 190     continue
c
c         i-th column of u
 200     indx1 = ju(i)
         indx2 = ju(i+1) - 1
         if(indx1 .gt. indx2) go to 260
         do 250 indx=indx1,indx2
            if(indx .eq. indx1) go to 240
            indxc1 = indx1
            indxc2 = indx - 1
            indxr1 = il(iu(indx))
            indxr2 = il(iu(indx)+1) - 1
            if(indxr1 .gt. indxr2) go to 240
 210        kr = jl(indxr1)
 220        kc = iu(indxc1)
            if(kr .gt. kc) then
               indxc1 = indxc1 + 1
               if(indxc1 .le. indxc2) go to 220
            elseif(kr .lt. kc) then
               indxr1 = indxr1 + 1
               if(indxr1 .le. indxr2) go to 210
            elseif(kr .eq. kc) then
               u(indx) = u(indx) - l(indxr1)*dinv(kc)*u(indxc1)
               indxr1 = indxr1 + 1
               indxc1 = indxc1 + 1
               if(indxr1 .le. indxr2 .and. indxc1 .le. indxc2) go to 210
            endif
 240        u(indx) = u(indx)/dinv(iu(indx))
 250     continue
c
c         i-th diagonal element
 260     indxr1 = il(i)
         indxr2 = il(i+1) - 1
         if(indxr1 .gt. indxr2) go to 300
         indxc1 = ju(i)
         indxc2 = ju(i+1) - 1
         if(indxc1 .gt. indxc2) go to 300
 270     kr = jl(indxr1)
 280     kc = iu(indxc1)
         if(kr .gt. kc) then
            indxc1 = indxc1 + 1
            if(indxc1 .le. indxc2) go to 280
         elseif(kr .lt. kc) then
            indxr1 = indxr1 + 1
            if(indxr1 .le. indxr2) go to 270
         elseif(kr .eq. kc) then
            dinv(i) = dinv(i) - l(indxr1)*dinv(kc)*u(indxc1)
            indxr1 = indxr1 + 1
            indxc1 = indxc1 + 1
            if(indxr1 .le. indxr2 .and. indxc1 .le. indxc2) go to 270
         endif
c
 300  continue
c
c         replace diagonal elements by their inverses.
cvd$ vector
      do 430 i=1,n
         dinv(i) = 1.0e0/dinv(i)
 430  continue
c
      return
c------------- last line of ssilus follows ----------------------------
      end
