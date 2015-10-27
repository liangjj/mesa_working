*deck ssdi
      subroutine ssdi (n, b, x, nelt, ia, ja, a, isym, rwork, iwork)
c***begin prologue  ssdi
c***purpose  diagonal matrix vector multiply.
c            routine to calculate the product  x = diag*b, where diag
c            is a diagonal matrix.
c***library   slatec (slap)
c***category  d1b4
c***type      single precision (ssdi-s, dsdi-d)
c***keywords  iterative precondition, linear system solve, slap, sparse
c***author  greenbaum, anne, (courant institute)
c           seager, mark k., (llnl)
c             lawrence livermore national laboratory
c             po box 808, l-60
c             livermore, ca 94550 (510) 423-3141
c             seager@llnl.gov
c***description
c
c *usage:
c     integer  n, nelt, ia(nelt), ja(nelt), isym, iwork(10)
c     real b(n), x(n), a(nelt), rwork(user defined)
c
c     call ssdi (n, b, x, nelt, ia, ja, a, isym, rwork, iwork)
c
c *arguments:
c n      :in       integer
c         order of the matrix.
c b      :in       real b(n).
c         vector to multiply the diagonal by.
c x      :out      real x(n).
c         result of diag*b.
c nelt   :dummy    integer.
c ia     :dummy    integer ia(nelt).
c ja     :dummy    integer ja(nelt).
c a      :dummy    real a(nelt).
c isym   :dummy    integer.
c         these are for compatibility with slap msolve calling sequence.
c rwork  :in       real rwork(user defined).
c         work array holding the diagonal of some matrix to scale
c         b by.  this array must be set by the user or by a call
c         to the slap routine ssds or ssd2s.  the length of rwork
c         must be >= iwork(4)+n.
c iwork  :in       integer iwork(10).
c         iwork(4) holds the offset into rwork for the diagonal matrix
c         to scale b by.  this is usually set up by the slap pre-
c         conditioner setup routines ssds or ssd2s.
c
c *description:
c         this routine is supplied with the slap package to perform
c         the  msolve  operation for iterative drivers that require
c         diagonal  scaling  (e.g., ssdcg, ssdbcg).   it  conforms
c         to the slap msolve calling convention  and hence does not
c         require an interface routine as do some of the other pre-
c         conditioners supplied with slap.
c
c***see also  ssds, ssd2s
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
c***end prologue  ssdi
c     .. scalar arguments ..
      integer isym, n, nelt
c     .. array arguments ..
      real a(nelt), b(n), rwork(*), x(n)
      integer ia(nelt), iwork(10), ja(nelt)
c     .. local scalars ..
      integer i, locd
c***first executable statement  ssdi
c
c         determine where the inverse of the diagonal
c         is in the work array and then scale by it.
c
      locd = iwork(4) - 1
      do 10 i = 1, n
         x(i) = rwork(locd+i)*b(i)
 10   continue
      return
c------------- last line of ssdi follows ----------------------------
      end
