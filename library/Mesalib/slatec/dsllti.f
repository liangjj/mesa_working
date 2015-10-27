*deck dsllti
      subroutine dsllti (n, b, x, nelt, ia, ja, a, isym, rwork, iwork)
c***begin prologue  dsllti
c***purpose  slap msolve for ldl' (ic) factorization.
c            this routine acts as an interface between the slap generic
c            msolve calling convention and the routine that actually
c                           -1
c            computes (ldl')  b = x.
c***library   slatec (slap)
c***category  d2e
c***type      double precision (ssllti-s, dsllti-d)
c***keywords  iterative precondition, linear system solve, slap, sparse
c***author  greenbaum, anne, (courant institute)
c           seager, mark k., (llnl)
c             lawrence livermore national laboratory
c             po box 808, l-60
c             livermore, ca 94550 (510) 423-3141
c             seager@llnl.gov
c***description
c       it is assumed that rwork and iwork have initialized with
c       the information required for dllti2:
c          iwork(1) = nel
c          iwork(2) = starting location of iel in iwork.
c          iwork(3) = starting location of jel in iwork.
c          iwork(4) = starting location of el in rwork.
c          iwork(5) = starting location of dinv in rwork.
c       see the description of dllti2 for details.
c***references  (none)
c***routines called  dllti2
c***revision history  (yymmdd)
c   871119  date written
c   881213  previous revision date
c   890915  made changes requested at july 1989 cml meeting.  (mks)
c   890922  numerous changes to prologue to make closer to slatec
c           standard.  (fnf)
c   890929  numerous changes to reduce sp/dp differences.  (fnf)
c   910411  prologue converted to version 4.0 format.  (bab)
c   910502  corrected conversion error.  (fnf)
c   920511  added complete declaration section.  (wrb)
c   921113  corrected c***category line.  (fnf)
c   930701  updated category section.  (fnf, wrb)
c***end prologue  dsllti
c     .. scalar arguments ..
      integer isym, n, nelt
c     .. array arguments ..
      double precision a(nelt), b(*), rwork(*), x(*)
      integer ia(nelt), iwork(*), ja(nelt)
c     .. local scalars ..
      integer locdin, locel, lociel, locjel, nel
c     .. external subroutines ..
      external dllti2
c***first executable statement  dsllti
      nel = iwork(1)
      lociel = iwork(3)
      locjel = iwork(2)
      locel  = iwork(4)
      locdin = iwork(5)
      call dllti2(n, b, x, nel, iwork(lociel), iwork(locjel),
     $     rwork(locel), rwork(locdin))
c
      return
c------------- last line of dsllti follows ----------------------------
      end
