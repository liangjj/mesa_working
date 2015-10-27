*deck sorth
      subroutine sorth (vnew, v, hes, n, ll, ldhes, kmp, snormw)
c***begin prologue  sorth
c***subsidiary
c***purpose  internal routine for sgmres.
c***library   slatec (slap)
c***category  d2a4, d2b4
c***type      single precision (sorth-s, dorth-d)
c***keywords  generalized minimum residual, iterative precondition,
c             non-symmetric linear system, slap, sparse
c***author  brown, peter, (llnl), pnbrown@llnl.gov
c           hindmarsh, alan, (llnl), alanh@llnl.gov
c           seager, mark k., (llnl), seager@llnl.gov
c             lawrence livermore national laboratory
c             po box 808, l-60
c             livermore, ca 94550 (510) 423-3141
c***description
c        this routine  orthogonalizes  the  vector  vnew  against the
c        previous kmp  vectors in the   v array.  it uses  a modified
c        gram-schmidt   orthogonalization procedure with  conditional
c        reorthogonalization.
c
c *usage:
c      integer n, ll, ldhes, kmp
c      real vnew(n), v(n,ll), hes(ldhes,ll), snormw
c
c      call sorth(vnew, v, hes, n, ll, ldhes, kmp, snormw)
c
c *arguments:
c vnew   :inout    real vnew(n)
c         on input, the vector of length n containing a scaled
c         product of the jacobian and the vector v(*,ll).
c         on output, the new vector orthogonal to v(*,i0) to v(*,ll),
c         where i0 = max(1, ll-kmp+1).
c v      :in       real v(n,ll)
c         the n x ll array containing the previous ll
c         orthogonal vectors v(*,1) to v(*,ll).
c hes    :inout    real hes(ldhes,ll)
c         on input, an ll x ll upper hessenberg matrix containing,
c         in hes(i,k), k.lt.ll, the scaled inner products of
c         a*v(*,k) and v(*,i).
c         on return, column ll of hes is filled in with
c         the scaled inner products of a*v(*,ll) and v(*,i).
c n      :in       integer
c         the order of the matrix a, and the length of vnew.
c ll     :in       integer
c         the current order of the matrix hes.
c ldhes  :in       integer
c         the leading dimension of the hes array.
c kmp    :in       integer
c         the number of previous vectors the new vector vnew
c         must be made orthogonal to (kmp .le. maxl).
c snormw :out      real
c         scalar containing the l-2 norm of vnew.
c
c***see also  sgmres
c***routines called  saxpy, sdot, snrm2
c***revision history  (yymmdd)
c   871001  date written
c   881213  previous revision date
c   890915  made changes requested at july 1989 cml meeting.  (mks)
c   890922  numerous changes to prologue to make closer to slatec
c           standard.  (fnf)
c   890929  numerous changes to reduce sp/dp differences.  (fnf)
c   910411  prologue converted to version 4.0 format.  (bab)
c   910506  made subsidiary to sgmres.  (fnf)
c   920511  added complete declaration section.  (wrb)
c***end prologue  sorth
c         the following is for optimized compilation on llnl/ltss crays.
clll. optimize
c     .. scalar arguments ..
      real snormw
      integer kmp, ldhes, ll, n
c     .. array arguments ..
      real hes(ldhes,*), v(n,*), vnew(*)
c     .. local scalars ..
      real arg, sumdsq, tem, vnrm
      integer i, i0
c     .. external functions ..
      real sdot, snrm2
      external sdot, snrm2
c     .. external subroutines ..
      external saxpy
c     .. intrinsic functions ..
      intrinsic max, sqrt
c***first executable statement  sorth
c
c         get norm of unaltered vnew for later use.
c
      vnrm = snrm2(n, vnew, 1)
c   -------------------------------------------------------------------
c         perform the modified gram-schmidt procedure on vnew =a*v(ll).
c         scaled inner products give new column of hes.
c         projections of earlier vectors are subtracted from vnew.
c   -------------------------------------------------------------------
      i0 = max(1,ll-kmp+1)
      do 10 i = i0,ll
         hes(i,ll) = sdot(n, v(1,i), 1, vnew, 1)
         tem = -hes(i,ll)
         call saxpy(n, tem, v(1,i), 1, vnew, 1)
 10   continue
c   -------------------------------------------------------------------
c         compute snormw = norm of vnew.  if vnew is small compared
c         to its input value (in norm), then reorthogonalize vnew to
c         v(*,1) through v(*,ll).  correct if relative correction
c         exceeds 1000*(unit roundoff).  finally, correct snormw using
c         the dot products involved.
c   -------------------------------------------------------------------
      snormw = snrm2(n, vnew, 1)
      if (vnrm + 0.001e0*snormw .ne. vnrm) return
      sumdsq = 0
      do 30 i = i0,ll
         tem = -sdot(n, v(1,i), 1, vnew, 1)
         if (hes(i,ll) + 0.001e0*tem .eq. hes(i,ll)) go to 30
         hes(i,ll) = hes(i,ll) - tem
         call saxpy(n, tem, v(1,i), 1, vnew, 1)
         sumdsq = sumdsq + tem**2
 30   continue
      if (sumdsq .eq. 0.0e0) return
      arg = max(0.0e0,snormw**2 - sumdsq)
      snormw = sqrt(arg)
c
      return
c------------- last line of sorth follows ----------------------------
      end
