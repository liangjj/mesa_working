*deck dorth
      subroutine dorth (vnew, v, hes, n, ll, ldhes, kmp, snormw)
c***begin prologue  dorth
c***subsidiary
c***purpose  internal routine for dgmres.
c***library   slatec (slap)
c***category  d2a4, d2b4
c***type      double precision (sorth-s, dorth-d)
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
c      double precision vnew(n), v(n,ll), hes(ldhes,ll), snormw
c
c      call dorth(vnew, v, hes, n, ll, ldhes, kmp, snormw)
c
c *arguments:
c vnew   :inout    double precision vnew(n)
c         on input, the vector of length n containing a scaled
c         product of the jacobian and the vector v(*,ll).
c         on output, the new vector orthogonal to v(*,i0) to v(*,ll),
c         where i0 = max(1, ll-kmp+1).
c v      :in       double precision v(n,ll)
c         the n x ll array containing the previous ll
c         orthogonal vectors v(*,1) to v(*,ll).
c hes    :inout    double precision hes(ldhes,ll)
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
c snormw :out      double precision
c         scalar containing the l-2 norm of vnew.
c
c***see also  dgmres
c***routines called  daxpy, ddot, dnrm2
c***revision history  (yymmdd)
c   890404  date written
c   890404  previous revision date
c   890915  made changes requested at july 1989 cml meeting.  (mks)
c   890922  numerous changes to prologue to make closer to slatec
c           standard.  (fnf)
c   890929  numerous changes to reduce sp/dp differences.  (fnf)
c   910411  prologue converted to version 4.0 format.  (bab)
c   910506  made subsidiary to dgmres.  (fnf)
c   920511  added complete declaration section.  (wrb)
c***end prologue  dorth
c         the following is for optimized compilation on llnl/ltss crays.
clll. optimize
c     .. scalar arguments ..
      double precision snormw
      integer kmp, ldhes, ll, n
c     .. array arguments ..
      double precision hes(ldhes,*), v(n,*), vnew(*)
c     .. local scalars ..
      double precision arg, sumdsq, tem, vnrm
      integer i, i0
c     .. external functions ..
      double precision ddot, dnrm2
      external ddot, dnrm2
c     .. external subroutines ..
      external daxpy
c     .. intrinsic functions ..
      intrinsic max, sqrt
c***first executable statement  dorth
c
c         get norm of unaltered vnew for later use.
c
      vnrm = dnrm2(n, vnew, 1)
c   -------------------------------------------------------------------
c         perform the modified gram-schmidt procedure on vnew =a*v(ll).
c         scaled inner products give new column of hes.
c         projections of earlier vectors are subtracted from vnew.
c   -------------------------------------------------------------------
      i0 = max(1,ll-kmp+1)
      do 10 i = i0,ll
         hes(i,ll) = ddot(n, v(1,i), 1, vnew, 1)
         tem = -hes(i,ll)
         call daxpy(n, tem, v(1,i), 1, vnew, 1)
 10   continue
c   -------------------------------------------------------------------
c         compute snormw = norm of vnew.  if vnew is small compared
c         to its input value (in norm), then reorthogonalize vnew to
c         v(*,1) through v(*,ll).  correct if relative correction
c         exceeds 1000*(unit roundoff).  finally, correct snormw using
c         the dot products involved.
c   -------------------------------------------------------------------
      snormw = dnrm2(n, vnew, 1)
      if (vnrm + 0.001d0*snormw .ne. vnrm) return
      sumdsq = 0
      do 30 i = i0,ll
         tem = -ddot(n, v(1,i), 1, vnew, 1)
         if (hes(i,ll) + 0.001d0*tem .eq. hes(i,ll)) go to 30
         hes(i,ll) = hes(i,ll) - tem
         call daxpy(n, tem, v(1,i), 1, vnew, 1)
         sumdsq = sumdsq + tem**2
 30   continue
      if (sumdsq .eq. 0.0d0) return
      arg = max(0.0d0,snormw**2 - sumdsq)
      snormw = sqrt(arg)
c
      return
c------------- last line of dorth follows ----------------------------
      end
