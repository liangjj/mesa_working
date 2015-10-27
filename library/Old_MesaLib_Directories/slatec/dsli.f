*deck dsli
      subroutine dsli (n, b, x, nelt, ia, ja, a, isym, rwork, iwork)
c***begin prologue  dsli
c***purpose  slap msolve for lower triangle matrix.
c            this routine acts as an interface between the slap generic
c            msolve calling convention and the routine that actually
c                      -1
c            computes l  b = x.
c***library   slatec (slap)
c***category  d2a3
c***type      double precision (ssli-s, dsli-d)
c***keywords  iterative precondition, linear system solve, slap, sparse
c***author  greenbaum, anne, (courant institute)
c           seager, mark k., (llnl)
c             lawrence livermore national laboratory
c             po box 808, l-60
c             livermore, ca 94550 (510) 423-3141
c             seager@llnl.gov
c***description
c       it is assumed that rwork and iwork have initialized with
c       the information required for dsli2:
c          iwork(1) = nel
c          iwork(2) = starting location of iel in iwork.
c          iwork(3) = starting location of jel in iwork.
c          iwork(4) = starting location of el in rwork.
c       see the description of dsli2 for details.
c***references  (none)
c***routines called  dsli2
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
c***end prologue  dsli
c     .. scalar arguments ..
      integer isym, n, nelt
c     .. array arguments ..
      double precision a(nelt), b(n), rwork(*), x(n)
      integer ia(nelt), iwork(10), ja(nelt)
c     .. local scalars ..
      integer locel, lociel, locjel, nel
c     .. external subroutines ..
      external dsli2
c***first executable statement  dsli
c
      nel = iwork(1)
      lociel = iwork(2)
      locjel = iwork(3)
      locel = iwork(4)
      call dsli2(n, b, x, nel, iwork(lociel), iwork(locjel),
     $     rwork(locel))
c
      return
c------------- last line of dsli follows ----------------------------
      end
