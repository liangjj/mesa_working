*deck dsdscl
      subroutine dsdscl (n, nelt, ia, ja, a, isym, x, b, dinv, job,
     +   itol)
c***begin prologue  dsdscl
c***purpose  diagonal scaling of system ax = b.
c            this routine scales (and unscales) the system  ax = b
c            by symmetric diagonal scaling.
c***library   slatec (slap)
c***category  d2e
c***type      double precision (ssdscl-s, dsdscl-d)
c***keywords  diagonal, slap sparse
c***author  greenbaum, anne, (courant institute)
c           seager, mark k., (llnl)
c             lawrence livermore national laboratory
c             po box 808, l-60
c             livermore, ca 94550 (510) 423-3141
c             seager@llnl.gov
c***description
c
c    this routine scales (and unscales) the system ax = b by symmetric
c    diagonal scaling.  the new system is:
c         -1/2  -1/2  1/2      -1/2
c        d    ad    (d   x) = d    b
c    when scaling is selected with the job parameter.  when unscaling
c    is selected this process is reversed.  the true solution is also
c    scaled or unscaled if itol is set appropriately, see below.
c
c *usage:
c     integer n, nelt, ia(nelt), ja(nelt), isym, job, itol
c     double precision a(nelt), x(n), b(n), dinv(n)
c
c     call dsdscl( n, nelt, ia, ja, a, isym, x, b, dinv, job, itol )
c
c *arguments:
c n      :in       integer
c         order of the matrix.
c nelt   :in       integer.
c         number of elements in arrays ia, ja, and a.
c ia     :in       integer ia(nelt).
c ja     :in       integer ja(nelt).
c a      :in       double precision a(nelt).
c         these arrays should hold the matrix a in the slap column
c         format.  see "description", below.
c isym   :in       integer.
c         flag to indicate symmetric storage format.
c         if isym=0, all non-zero entries of the matrix are stored.
c         if isym=1, the matrix is symmetric, and only the upper
c         or lower triangle of the matrix is stored.
c x      :inout    double precision x(n).
c         initial guess that will be later used in the iterative
c         solution.
c         of the scaled system.
c b      :inout    double precision b(n).
c         right hand side vector.
c dinv   :inout    double precision dinv(n).
c         upon return this array holds 1./diag(a).
c         this is an input if job = 0.
c job    :in       integer.
c         flag indicating whether to scale or not.
c         job non-zero means do scaling.
c         job = 0 means do unscaling.
c itol   :in       integer.
c         flag indicating what type of error estimation to do in the
c         iterative method.  when itol = 11 the exact solution from
c         common block dslblk will be used.  when the system is scaled
c         then the true solution must also be scaled.  if itol is not
c         11 then this vector is not referenced.
c
c *common blocks:
c soln    :inout   double precision soln(n).  common block /dslblk/
c         the true solution, soln, is scaled (or unscaled) if itol is
c         set to 11, see above.
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
c       and that the operation dinv = 1.0/diag(a)  will  not  under-
c       flow or overflow. this is done so that the loop  vectorizes.
c       matrices  with zero or near zero or very  large entries will
c       have numerical difficulties  and  must  be fixed before this
c       routine is called.
c
c***see also  dsdcg
c***references  (none)
c***routines called  (none)
c***common blocks    dslblk
c***revision history  (yymmdd)
c   890404  date written
c   890404  previous revision date
c   890915  made changes requested at july 1989 cml meeting.  (mks)
c   890922  numerous changes to prologue to make closer to slatec
c           standard.  (fnf)
c   890929  numerous changes to reduce sp/dp differences.  (fnf)
c   910411  prologue converted to version 4.0 format.  (bab)
c   910502  added c***first executable statement line.  (fnf)
c   920407  common block renamed dslblk.  (wrb)
c   920511  added complete declaration section.  (wrb)
c   921113  corrected c***category line.  (fnf)
c   930701  updated category section.  (fnf, wrb)
c***end prologue  dsdscl
c     .. scalar arguments ..
      integer isym, itol, job, n, nelt
c     .. array arguments ..
      double precision a(nelt), b(n), dinv(n), x(n)
      integer ia(nelt), ja(nelt)
c     .. arrays in common ..
      double precision soln(1)
c     .. local scalars ..
      double precision di
      integer icol, j, jbgn, jend
c     .. intrinsic functions ..
      intrinsic sqrt
c     .. common blocks ..
      common /dslblk/ soln
c***first executable statement  dsdscl
c
c         scaling...
c
      if( job.ne.0 ) then
         do 10 icol = 1, n
            dinv(icol) = 1.0d0/sqrt( a(ja(icol)) )
 10      continue
      else
c
c         unscaling...
c
         do 15 icol = 1, n
            dinv(icol) = 1.0d0/dinv(icol)
 15      continue
      endif
c
      do 30 icol = 1, n
         jbgn = ja(icol)
         jend = ja(icol+1)-1
         di = dinv(icol)
         do 20 j = jbgn, jend
            a(j) = dinv(ia(j))*a(j)*di
 20      continue
 30   continue
c
      do 40 icol = 1, n
         b(icol) = b(icol)*dinv(icol)
         x(icol) = x(icol)/dinv(icol)
 40   continue
c
c         check to see if we need to scale the "true solution" as well.
c
      if( itol.eq.11 ) then
         do 50 icol = 1, n
            soln(icol) = soln(icol)/dinv(icol)
 50      continue
      endif
c
      return
c------------- last line of dsdscl follows ----------------------------
      end
