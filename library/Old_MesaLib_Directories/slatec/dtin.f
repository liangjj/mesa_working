*deck dtin
      subroutine dtin (n, nelt, ia, ja, a, isym, soln, rhs, iunit, job)
c***begin prologue  dtin
c***purpose  read in slap triad format linear system.
c            routine to read in a slap triad format matrix and right
c            hand side and solution to the system, if known.
c***library   slatec (slap)
c***category  n1
c***type      double precision (stin-s, dtin-d)
c***keywords  diagnostics, linear system, slap sparse
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
c     call dtin( n, nelt, ia, ja, a, isym, soln, rhs, iunit, job )
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
c         is accessed if and only if job to read it in, see below.
c         if the user requests that soln be read in, but it is not in
c         the file, then it is simply zeroed out.
c rhs    :out      double precision rhs(n).
c         the right hand side vector.  this array is accessed if and
c         only if job is set to read it in, see below.
c         if the user requests that rhs be read in, but it is not in
c         the file, then it is simply zeroed out.
c iunit  :in       integer.
c         fortran logical i/o device unit number to write the matrix
c         to.  this unit must be connected in a system dependent fashion
c         to a file or the console or you will get a nasty message
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
c         job = 0 => read in only the matrix.
c         job = 1 => read in the matrix and rhs.
c         job = 2 => read in the matrix and soln.
c         job = 3 => read in the matrix, rhs and soln.
c
c *description:
c       the format for the  input is as follows.  on  the first line
c       are counters and flags: n, nelt, isym, irhs, isoln.  n, nelt
c       and isym are described above.  irhs is  a flag indicating if
c       the rhs was  written out (1 is  yes, 0 is  no).  isoln  is a
c       flag indicating if the soln was written out  (1 is yes, 0 is
c       no).  the format for the fist line is: 5i10.  then comes the
c       nelt triad's ia(i), ja(i) and a(i), i = 1, nelt.  the format
c       for  these lines is   :  1x,i5,1x,i5,1x,d16.7.   then  comes
c       rhs(i), i = 1, n, if irhs = 1.  then  comes soln(i), i  = 1,
c       n, if isoln = 1.  the format for these lines is: 1x,d16.7.
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
c *cautions:
c     this routine will attempt to write to the fortran logical output
c     unit iunit, if iunit .ne. 0.  thus, the user must make sure that
c     this logical unit is attached to a file or terminal before calling
c     this routine with a non-zero value for iunit.  this routine does
c     not check for the validity of a non-zero iunit unit number.
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
c   921007  changed e's to d's in formats.  (fnf)
c   930701  updated category section.  (fnf, wrb)
c***end prologue  dtin
c     .. scalar arguments ..
      integer isym, iunit, job, n, nelt
c     .. array arguments ..
      double precision a(nelt), rhs(n), soln(n)
      integer ia(nelt), ja(nelt)
c     .. local scalars ..
      integer i, irhs, isoln, jobret, neltmx
c     .. intrinsic functions ..
      intrinsic min
c***first executable statement  dtin
c
c         read in the information heading.
c
      neltmx = nelt
      read(iunit,1000) n, nelt, isym, irhs, isoln
      nelt = min( nelt, neltmx )
c
c         read in the matrix non-zeros in triad format.
      do 10 i = 1, nelt
         read(iunit,1010) ia(i), ja(i), a(i)
 10   continue
c
c         if requested, read in the rhs.
      jobret = 0
      if( job.eq.1 .or. job.eq.3 ) then
c
c         check to see if rhs is in the file.
         if( irhs.eq.1 ) then
            jobret = 1
            read(iunit,1020) (rhs(i),i=1,n)
         else
            do 20 i = 1, n
               rhs(i) = 0
 20         continue
         endif
      endif
c
c         if requested, read in the solution.
      if( job.gt.1 ) then
c
c         check to see if solution is in the file.
         if( isoln.eq.1 ) then
            jobret = jobret + 2
            read(iunit,1020) (soln(i),i=1,n)
         else
            do 30 i = 1, n
               soln(i) = 0
 30         continue
         endif
      endif
c
      job = jobret
      return
 1000 format(5i10)
 1010 format(1x,i5,1x,i5,1x,d16.7)
 1020 format(1x,d16.7)
c------------- last line of dtin follows ----------------------------
      end
