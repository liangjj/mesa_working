*deck dcpplt
      subroutine dcpplt (n, nelt, ia, ja, a, isym, iunit)
c***begin prologue  dcpplt
c***purpose  printer plot of slap column format matrix.
c            routine to print out a slap column format matrix in a
c            "printer plot" graphical representation.
c***library   slatec (slap)
c***category  n1
c***type      double precision (scpplt-s, dcpplt-d)
c***keywords  diagnostics, linear system, slap sparse
c***author  seager, mark k., (llnl)
c             lawrence livermore national laboratory
c             po box 808, l-60
c             livermore, ca 94550 (510) 423-3141
c             seager@llnl.gov
c***description
c
c *usage:
c     integer n, nelt, ia(nelt), ja(nelt), isym, iunit
c     double precision a(nelt)
c
c     call dcpplt( n, nelt, ia, ja, a, isym, iunit )
c
c *arguments:
c n      :in       integer
c         order of the matrix.
c         if n.gt.maxord, only the leading maxord x maxord
c         submatrix will be printed.  (currently maxord = 225.)
c nelt   :in       integer.
c         number of non-zeros stored in a.
c ia     :in       integer ia(nelt).
c ja     :in       integer ja(nelt).
c a      :in       double precision a(nelt).
c         these arrays should hold the matrix a in the slap
c         column format.  see "description", below.
c isym   :in       integer.
c         flag to indicate symmetric storage format.
c         if isym=0, all non-zero entries of the matrix are stored.
c         if isym=1, the matrix is symmetric, and only the lower
c         triangle of the matrix is stored.
c iunit  :in       integer.
c         fortran logical i/o device unit number to write the matrix
c         to.  this unit must be connected in a system dependent fashion
c         to a file or the console or you will get a nasty message
c         from the fortran i/o libraries.
c
c *description:
c       this routine prints out a slap  column format matrix  to the
c       fortran logical i/o unit   number  iunit.  the  numbers them
c       selves  are not printed  out, but   rather  a one  character
c       representation of the numbers.   elements of the matrix that
c       are not represented in the (ia,ja,a)  arrays are  denoted by
c       ' ' character (a blank).  elements of a that are *zero* (and
c       hence  should  really not be  stored) are  denoted  by a '0'
c       character.  elements of a that are *positive* are denoted by
c       'd' if they are diagonal elements  and '#' if  they are off
c       diagonal  elements.  elements of  a that are *negative* are
c       denoted by 'n'  if they  are diagonal  elements and  '*' if
c       they are off diagonal elements.
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
c *cautions:
c     this routine will attempt to write to the fortran logical output
c     unit iunit, if iunit .ne. 0.  thus, the user must make sure that
c     this logical unit is attached to a file or terminal before calling
c     this routine with a non-zero value for iunit.  this routine does
c     not check for the validity of a non-zero iunit unit number.
c
c *portability:
c     this routine, as distributed, can generate lines up to 229
c     characters long.  some fortran systems have more restricted
c     line lengths.  change parameter maxord and the large number
c     in format 1010 to reduce this line length.
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
c   921007  replaced hard-wired 225 with parameter maxord.  (fnf)
c   921021  corrected syntax of character declaration.  (fnf)
c   921026  corrected d to e in output format.  (fnf)
c   930701  updated category section.  (fnf, wrb)
c***end prologue  dcpplt
c     .. scalar arguments ..
      integer isym, iunit, n, nelt
c     .. array arguments ..
      double precision a(nelt)
      integer ia(nelt), ja(nelt)
c     .. parameters ..
      integer  maxord
      parameter (maxord=225)
c     .. local scalars ..
      integer i, icol, irow, j, jbgn, jend, nmax
c     .. local arrays ..
      character chmat(maxord)*(maxord)
c     .. intrinsic functions ..
      intrinsic min, mod, real
c***first executable statement  dcpplt
c
c         set up the character matrix...
c
      nmax = min( maxord, n )
      do 10 i = 1, nmax
         chmat(i)(1:nmax) = ' '
 10   continue
      do 30 icol = 1, nmax
         jbgn = ja(icol)
         jend = ja(icol+1)-1
         do 20 j = jbgn, jend
            irow = ia(j)
            if( irow.le.nmax ) then
               if( isym.ne.0 ) then
c         put in non-sym part as well...
                  if( a(j).eq.0.0d0 ) then
                     chmat(irow)(icol:icol) = '0'
                  elseif( a(j).gt.0.0d0 ) then
                     chmat(irow)(icol:icol) = '#'
                  else
                     chmat(irow)(icol:icol) = '*'
                  endif
               endif
               if( irow.eq.icol ) then
c         diagonal entry.
                  if( a(j).eq.0.0d0 ) then
                     chmat(irow)(icol:icol) = '0'
                  elseif( a(j).gt.0.0d0 ) then
                     chmat(irow)(icol:icol) = 'd'
                  else
                     chmat(irow)(icol:icol) = 'n'
                  endif
               else
c         off-diagonal entry
                  if( a(j).eq.0.0d0 ) then
                     chmat(irow)(icol:icol) = '0'
                  elseif( a(j).gt.0.0d0 ) then
                     chmat(irow)(icol:icol) = '#'
                  else
                     chmat(irow)(icol:icol) = '*'
                  endif
               endif
            endif
 20      continue
 30   continue
c
c         write out the heading.
      write(iunit,1000) n, nelt, real(nelt)/(n*n)
      write(iunit,1010) (mod(i,10),i=1,nmax)
c
c         write out the character representations matrix elements.
      do 40 irow = 1, nmax
         write(iunit,1020) irow, chmat(irow)(1:nmax)
 40   continue
      return
c
 1000 format(/'**** picture of column slap matrix follows ****'/
     $     ' n, nelt and density = ',2i10,d16.7)
c      the following assumes maxord.le.225.
 1010 format(4x,225(i1))
 1020 format(1x,i3,a)
c------------- last line of dcpplt follows ----------------------------
      end
