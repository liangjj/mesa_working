*deck dbhin
      subroutine dbhin (n, nelt, ia, ja, a, isym, soln, rhs, iunit, job)
c***begin prologue  dbhin
c***purpose  read a sparse linear system in the boeing/harwell format.
c            the matrix is read in and if the right hand side is also
c            present in the input file then it too is read in.  the
c            matrix is then modified to be in the slap column format.
c***library   slatec (slap)
c***category  n1
c***type      double precision (sbhin-s, dbhin-d)
c***keywords  linear system, matrix read, slap sparse
c***author  seager, mark k., (llnl)
c             lawrence livermore national laboratory
c             po box 808, l-60
c             livermore, ca 94550 (510) 423-3141
c             seager@llnl.gov
c***description
c
c *usage:
c     integer n, nelt, ia(nelt), ja(nelt), isym, iunit, job
c     double precision a(nelt), soln(n), rhs(n)
c
c     call dbhin( n, nelt, ia, ja, a, isym, soln, rhs, iunit, job )
c
c *arguments:
c n      :out      integer
c         order of the matrix.
c nelt   :inout    integer.
c         on input nelt is the maximum number of non-zeros that
c         can be stored in the ia, ja, a arrays.
c         on output nelt is the number of non-zeros stored in a.
c ia     :out      integer ia(nelt).
c ja     :out      integer ja(nelt).
c a      :out      double precision a(nelt).
c         on output these arrays hold the matrix a in the slap
c         triad format.  see "description", below.
c isym   :out      integer.
c         flag to indicate symmetric storage format.
c         if isym=0, all non-zero entries of the matrix are stored.
c         if isym=1, the matrix is symmetric, and only the lower
c         triangle of the matrix is stored.
c soln   :out      double precision soln(n).
c         the solution to the linear system, if present.  this array
c         is accessed if and only if job is set to read it in, see
c         below.  if the user requests that soln be read in, but it is
c         not in the file, then it is simply zeroed out.
c rhs    :out      double precision rhs(n).
c         the right hand side vector.  this array is accessed if and
c         only if job is set to read it in, see below.
c         if the user requests that rhs be read in, but it is not in
c         the file, then it is simply zeroed out.
c iunit  :in       integer.
c         fortran logical i/o device unit number to read the matrix
c         from.  this unit must be connected in a system dependent
c         fashion to a file, or you will get a nasty message
c         from the fortran i/o libraries.
c job    :inout    integer.
c         flag indicating what i/o operations to perform.
c         on input job indicates what input operations to try to
c         perform.
c         job = 0 => read only the matrix.
c         job = 1 => read matrix and rhs (if present).
c         job = 2 => read matrix and soln (if present).
c         job = 3 => read matrix, rhs and soln (if present).
c         on output job indicates what operations were actually
c         performed.
c         job = -3 => unable to parse matrix "code" from input file
c                     to determine if only the lower triangle of matrix
c                     is stored.
c         job = -2 => number of non-zeros (nelt) too large.
c         job = -1 => system size (n) too large.
c         job =  0 => read in only the matrix.
c         job =  1 => read in the matrix and rhs.
c         job =  2 => read in the matrix and soln.
c         job =  3 => read in the matrix, rhs and soln.
c         job = 10 => read in only the matrix *structure*, but no
c                     non-zero entries.  hence, a(*) is not referenced
c                     and has the return values the same as the input.
c         job = 11 => read in the matrix *structure* and rhs.
c         job = 12 => read in the matrix *structure* and soln.
c         job = 13 => read in the matrix *structure*, rhs and soln.
c
c *description:
c       the format for the input is as follows.  the first line contains
c       a title to identify the data file.  on the second line (5i4) are
c       counters: nline, npls, nrils, nnvls, nrhsls.
c        nline  number of data lines (after the header) in the file.
c        npls   number of lines for the column pointer data in the file.
c        nrils  number of lines for the row indices in the file.
c        nnvls  number of lines for the matrix elements in the file.
c        nrhsls number of lines for the rhs in the file.
c       the third line (a3,11x,4i4) contains a symmetry code and some
c       additional counters: code, nrow, ncol, nind, nele.
c       on the fourth line (2a16,2a20) are formats to be used to read
c       the following data: pntfnt, rinfmt, nvlfmt, rhsfmt.
c       following that are the blocks of data in the order indicated.
c
c       =================== s l a p triad format ===================
c       this routine requires that the  matrix a be   stored in  the
c       slap  triad format.  in  this format only the non-zeros  are
c       stored.  they may appear in  *any* order.  the user supplies
c       three arrays of  length nelt, where  nelt is  the number  of
c       non-zeros in the matrix: (ia(nelt), ja(nelt), a(nelt)).  for
c       each non-zero the user puts the row and column index of that
c       matrix element  in the ia and  ja arrays.  the  value of the
c       non-zero  matrix  element is  placed   in  the corresponding
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
c *portability:
c         you must make sure that iunit is a valid fortran logical
c         i/o device unit number and that the unit number has been
c         associated with a file or the console.  this is a system
c         dependent function.
c
c *implementation note:
c         soln is not read by this version.  it will simply be
c         zeroed out if job = 2 or 3 and the returned value of
c         job will indicate soln has not been read.
c***references  (none)
c***routines called  (none)
c***revision history  (yymmdd)
c   881107  date written
c   881213  previous revision date
c   890915  made changes requested at july 1989 cml meeting.  (mks)
c   890922  numerous changes to prologue to make closer to slatec
c           standard.  (fnf)
c   890929  numerous changes to reduce sp/dp differences.  (fnf)
c   910411  prologue converted to version 4.0 format.  (bab)
c   911122  added loop to zero out rhs if user wants to read rhs, but
c           it's not in the input file. (mks)
c   911125  minor improvements to prologue.  (fnf)
c   920511  added complete declaration section.  (wrb)
c   921007  corrected description of input format.  (fnf)
c   921208  added implementation note and code to zero out soln.  (fnf)
c   930701  updated category section.  (fnf, wrb)
c***end prologue  dbhin
c     .. scalar arguments ..
      integer isym, iunit, job, n, nelt
c     .. array arguments ..
      double precision a(nelt), rhs(n), soln(n)
      integer ia(nelt), ja(nelt)
c     .. local scalars ..
      double precision temp
      integer i, ibgn, icol, iend, itemp, j, jobret, ncol, nele, nind,
     +        nline, nnvls, npls, nrhsls, nrils, nrow
      character code*3, pntfmt*16, rinfmt*16, nvlfmt*20, rhsfmt*20,
     +          title*80
c     .. intrinsic functions ..
      intrinsic mod
c***first executable statement  dbhin
c
c         read matrices in boeing-harwell format.
c
c title  header line to identify data file.
c nline  number of data lines (after the header) in the file.
c npls   number of lines for the column pointer data in the file.
c nrils  number of lines for the row indices in the data file.
c nnvls  number of lines for the matrix elements in the data file.
c nrhsls number of lines for the rhs in the data file.
c ---- only those variables needed by slap are referenced. ----
c
      read(iunit,9000) title
      read(iunit,9010) nline, npls, nrils, nnvls, nrhsls
      read(iunit,9020) code, nrow, ncol, nind, nele
      read(iunit,9030) pntfmt, rinfmt, nvlfmt, rhsfmt
c
      if( nrow.gt.n ) then
         n = nrow
         jobret = -1
         goto 999
      endif
      if( nind.gt.nelt ) then
         nelt = nind
         jobret = -2
         goto 999
      endif
c
c         set the parameters.
c
      n    = nrow
      nelt = nind
      if( code.eq.'rua' ) then
         isym = 0
      else if( code.eq.'rsa' ) then
         isym = 1
      else
         jobret = -3
         goto 999
      endif
      read(iunit,pntfmt) (ja(i), i = 1, n+1)
      read(iunit,rinfmt) (ia(i), i = 1, nelt)
      jobret = 10
      if( nnvls.gt.0 ) then
         read(iunit,nvlfmt) (a(i),  i = 1, nelt)
         jobret = 0
      endif
      if( mod(job,2).eq.1 ) then
c
c         user requests that the rhs be read in.  if it is in the input
c         file, read it in; otherwise just zero it out.
c
         if( nrhsls.gt.0 ) then
            read(5,rhsfmt) (rhs(i), i = 1, n)
            jobret = jobret + 1
         else
            do 10 i = 1, n
               rhs(i) = 0
 10         continue
         endif
      endif
      if ( (job.eq.2).or.(job.eq.3) ) then
c
c         user requests that the soln be read in.
c         just zero out the array.
c
         do 20 i = 1, n
            soln(i) = 0
 20      continue
      endif
c
c         now loop through the ia array making sure that the diagonal
c         matrix element appears first in the column.  then sort the
c         rest of the column in ascending order.
c
cvd$r noconcur
cvd$r novector
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
c
c         set return flag.
 999  job = jobret
      return
 9000 format( a80 )
 9010 format( 5i14 )
 9020 format( a3, 11x, 4i14 )
 9030 format( 2a16, 2a20 )
c------------- last line of dbhin follows ------------------------------
      end
