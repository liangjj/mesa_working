*deck dsluti
      subroutine dsluti (n, b, x, nelt, ia, ja, a, isym, rwork, iwork)
c***begin prologue  dsluti
c***purpose  slap mtsolv for ldu factorization.
c            this routine acts as an interface between the slap generic
c            mtsolv calling convention and the routine that actually
c                           -t
c            computes  (ldu)  b = x.
c***library   slatec (slap)
c***category  d2e
c***type      double precision (ssluti-s, dsluti-d)
c***keywords  iterative precondition, linear system solve, slap, sparse
c***author  greenbaum, anne, (courant institute)
c           seager, mark k., (llnl)
c             lawrence livermore national laboratory
c             po box 808, l-60
c             livermore, ca 94550 (510) 423-3141
c             seager@llnl.gov
c***description
c       it is assumed that rwork and iwork have initialized with
c       the information required for dslui4:
c          iwork(1) = starting location of il in iwork.
c          iwork(2) = starting location of jl in iwork.
c          iwork(3) = starting location of iu in iwork.
c          iwork(4) = starting location of ju in iwork.
c          iwork(5) = starting location of l in rwork.
c          iwork(6) = starting location of dinv in rwork.
c          iwork(7) = starting location of u in rwork.
c       see the description of dslui4 for details.
c***references  (none)
c***routines called  dslui4
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
c***end prologue  dsluti
c     .. scalar arguments ..
      integer isym, n, nelt
c     .. array arguments ..
      double precision a(n), b(n), rwork(*), x(n)
      integer ia(nelt), iwork(10), ja(nelt)
c     .. local scalars ..
      integer locdin, locil, lociu, locjl, locju, locl, locu
c     .. external subroutines ..
      external dslui4
c***first executable statement  dsluti
c
c         pull out the pointers to the l, d and u matrices and call
c         the workhorse routine.
c
      locil = iwork(1)
      locjl = iwork(2)
      lociu = iwork(3)
      locju = iwork(4)
      locl = iwork(5)
      locdin = iwork(6)
      locu = iwork(7)
c
      call dslui4(n, b, x, iwork(locil), iwork(locjl), rwork(locl),
     $     rwork(locdin), iwork(lociu), iwork(locju), rwork(locu))
c
      return
c------------- last line of dsluti follows ----------------------------
      end
