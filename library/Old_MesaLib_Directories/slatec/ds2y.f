*deck ds2y
      subroutine ds2y (n, nelt, ia, ja, a, isym)
c***begin prologue  ds2y
c***purpose  slap triad to slap column format converter.
c            routine to convert from the slap triad to slap column
c            format.
c***library   slatec (slap)
c***category  d1b9
c***type      double precision (ss2y-s, ds2y-d)
c***keywords  linear system, slap sparse
c***author  seager, mark k., (llnl)
c             lawrence livermore national laboratory
c             po box 808, l-60
c             livermore, ca 94550 (510) 423-3141
c             seager@llnl.gov
c***description
c
c *usage:
c     integer n, nelt, ia(nelt), ja(nelt), isym
c     double precision a(nelt)
c
c     call ds2y( n, nelt, ia, ja, a, isym )
c
c *arguments:
c n      :in       integer
c         order of the matrix.
c nelt   :in       integer.
c         number of non-zeros stored in a.
c ia     :inout    integer ia(nelt).
c ja     :inout    integer ja(nelt).
c a      :inout    double precision a(nelt).
c         these arrays should hold the matrix a in either the slap
c         triad format or the slap column format.  see "description",
c         below.  if the slap triad format is used, this format is
c         translated to the slap column format by this routine.
c isym   :in       integer.
c         flag to indicate symmetric storage format.
c         if isym=0, all non-zero entries of the matrix are stored.
c         if isym=1, the matrix is symmetric, and only the lower
c         triangle of the matrix is stored.
c
c *description:
c       the sparse linear algebra package (slap) utilizes two matrix
c       data structures: 1) the  slap triad  format or  2)  the slap
c       column format.  the user can hand this routine either of the
c       of these data structures.  if the slap triad format is give
c       as input then this routine transforms it into slap column
c       format.  the way this routine tells which format is given as
c       input is to look at ja(n+1).  if ja(n+1) = nelt+1 then we
c       have the slap column format.  if that equality does not hold
c       then it is assumed that the ia, ja, a arrays contain the
c       slap triad format.
c
c       =================== s l a p triad format ===================
c       this routine requires that the  matrix a be   stored in  the
c       slap  triad format.  in  this format only the non-zeros  are
c       stored.  they may appear in  *any* order.  the user supplies
c       three arrays of  length nelt, where  nelt is  the number  of
c       non-zeros in the matrix: (ia(nelt), ja(nelt), a(nelt)).  for
c       each non-zero the user puts the row and column index of that
c       matrix element  in the ia and  ja arrays.  the  value of the
c       non-zero   matrix  element is  placed  in  the corresponding
c       location of the a array.   this is  an  extremely  easy data
c       structure to generate.  on  the  other hand it   is  not too
c       efficient on vector computers for  the iterative solution of
c       linear systems.  hence,   slap changes   this  input    data
c       structure to the slap column format  for  the iteration (but
c       does not change it back).
c
c       here is an example of the  slap triad   storage format for a
c       5x5 matrix.  recall that the entries may appear in any order.
c
c           5x5 matrix      slap triad format for 5x5 matrix on left.
c                              1  2  3  4  5  6  7  8  9 10 11
c       |11 12  0  0 15|   a: 51 12 11 33 15 53 55 22 35 44 21
c       |21 22  0  0  0|  ia:  5  1  1  3  1  5  5  2  3  4  2
c       | 0  0 33  0 35|  ja:  1  2  1  3  5  3  5  2  5  4  1
c       | 0  0  0 44  0|
c       |51  0 53  0 55|
c
c       =================== s l a p column format ==================
c
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
c***references  (none)
c***routines called  qs2i1d
c***revision history  (yymmdd)
c   871119  date written
c   881213  previous revision date
c   890915  made changes requested at july 1989 cml meeting.  (mks)
c   890922  numerous changes to prologue to make closer to slatec
c           standard.  (fnf)
c   890929  numerous changes to reduce sp/dp differences.  (fnf)
c   910411  prologue converted to version 4.0 format.  (bab)
c   910502  corrected c***first executable statement line.  (fnf)
c   920511  added complete declaration section.  (wrb)
c   930701  updated category section.  (fnf, wrb)
c***end prologue  ds2y
c     .. scalar arguments ..
      integer isym, n, nelt
c     .. array arguments ..
      double precision a(nelt)
      integer ia(nelt), ja(nelt)
c     .. local scalars ..
      double precision temp
      integer i, ibgn, icol, iend, itemp, j
c     .. external subroutines ..
      external qs2i1d
c***first executable statement  ds2y
c
c         check to see if the (ia,ja,a) arrays are in slap column
c         format.  if it's not then transform from slap triad.
c
      if( ja(n+1).eq.nelt+1 ) return
c
c         sort into ascending order by column (on the ja array).
c         this will line up the columns.
c
      call qs2i1d( ja, ia, a, nelt, 1 )
c
c         loop over each column to see where the column indices change
c         in the column index array ja.  this marks the beginning of the
c         next column.
c
cvd$r novector
      ja(1) = 1
      do 20 icol = 1, n-1
         do 10 j = ja(icol)+1, nelt
            if( ja(j).ne.icol ) then
               ja(icol+1) = j
               goto 20
            endif
 10      continue
 20   continue
      ja(n+1) = nelt+1
c
c         mark the n+2 element so that future calls to a slap routine
c         utilizing the ysmp-column storage format will be able to tell.
c
      ja(n+2) = 0
c
c         now loop through the ia array making sure that the diagonal
c         matrix element appears first in the column.  then sort the
c         rest of the column in ascending order.
c
      do 70 icol = 1, n
         ibgn = ja(icol)
         iend = ja(icol+1)-1
         do 30 i = ibgn, iend
            if( ia(i).eq.icol ) then
c
c              swap the diagonal element with the first element in the
c              column.
c
               itemp = ia(i)
               ia(i) = ia(ibgn)
               ia(ibgn) = itemp
               temp = a(i)
               a(i) = a(ibgn)
               a(ibgn) = temp
               goto 40
            endif
 30      continue
 40      ibgn = ibgn + 1
         if( ibgn.lt.iend ) then
            do 60 i = ibgn, iend
               do 50 j = i+1, iend
                  if( ia(i).gt.ia(j) ) then
                     itemp = ia(i)
                     ia(i) = ia(j)
                     ia(j) = itemp
                     temp = a(i)
                     a(i) = a(j)
                     a(j) = temp
                  endif
 50            continue
 60         continue
         endif
 70   continue
      return
c------------- last line of ds2y follows ----------------------------
      end
