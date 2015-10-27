*deck drlcal
      subroutine drlcal (n, kmp, ll, maxl, v, q, rl, snormw, prod,
     +   r0nrm)
c***begin prologue  drlcal
c***subsidiary
c***purpose  internal routine for dgmres.
c***library   slatec (slap)
c***category  d2a4, d2b4
c***type      double precision (srlcal-s, drlcal-d)
c***keywords  generalized minimum residual, iterative precondition,
c             non-symmetric linear system, slap, sparse
c***author  brown, peter, (llnl), pnbrown@llnl.gov
c           hindmarsh, alan, (llnl), alanh@llnl.gov
c           seager, mark k., (llnl), seager@llnl.gov
c             lawrence livermore national laboratory
c             po box 808, l-60
c             livermore, ca 94550 (510) 423-3141
c***description
c         this routine calculates the scaled residual rl from the
c         v(i)'s.
c *usage:
c      integer n, kmp, ll, maxl
c      double precision v(n,ll), q(2*maxl), rl(n), snormw, prod, r0norm
c
c      call drlcal(n, kmp, ll, maxl, v, q, rl, snormw, prod, r0nrm)
c
c *arguments:
c n      :in       integer
c         the order of the matrix a, and the lengths
c         of the vectors sr, sz, r0 and z.
c kmp    :in       integer
c         the number of previous v vectors the new vector vnew
c         must be made orthogonal to. (kmp .le. maxl)
c ll     :in       integer
c         the current dimension of the krylov subspace.
c maxl   :in       integer
c         the maximum dimension of the krylov subspace.
c v      :in       double precision v(n,ll)
c         the n x ll array containing the orthogonal vectors
c         v(*,1) to v(*,ll).
c q      :in       double precision q(2*maxl)
c         a double precision array of length 2*maxl containing the
c         components of the givens rotations used in the qr
c         decomposition of hes.  it is loaded in dheqr and used in
c         dhels.
c rl     :out      double precision rl(n)
c         the residual vector rl.  this is either sb*(b-a*xl) if
c         not preconditioning or preconditioning on the right,
c         or sb*(m-inverse)*(b-a*xl) if preconditioning on the
c         left.
c snormw :in       double precision
c         scale factor.
c prod   :in       double precision
c         the product s1*s2*...*sl = the product of the sines of the
c         givens rotations used in the qr factorization of
c         the hessenberg matrix hes.
c r0nrm  :in       double precision
c         the scaled norm of initial residual r0.
c
c***see also  dgmres
c***routines called  dcopy, dscal
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
c***end prologue  drlcal
c         the following is for optimized compilation on llnl/ltss crays.
clll. optimize
c     .. scalar arguments ..
      double precision prod, r0nrm, snormw
      integer kmp, ll, maxl, n
c     .. array arguments ..
      double precision q(*), rl(n), v(n,*)
c     .. local scalars ..
      double precision c, s, tem
      integer i, i2, ip1, k, llm1, llp1
c     .. external subroutines ..
      external dcopy, dscal
c***first executable statement  drlcal
      if (kmp .eq. maxl) then
c
c         calculate rl.  start by copying v(*,1) into rl.
c
         call dcopy(n, v(1,1), 1, rl, 1)
         llm1 = ll - 1
         do 20 i = 1,llm1
            ip1 = i + 1
            i2 = i*2
            s = q(i2)
            c = q(i2-1)
            do 10 k = 1,n
               rl(k) = s*rl(k) + c*v(k,ip1)
 10         continue
 20      continue
         s = q(2*ll)
         c = q(2*ll-1)/snormw
         llp1 = ll + 1
         do 30 k = 1,n
            rl(k) = s*rl(k) + c*v(k,llp1)
 30      continue
      endif
c
c         when kmp < maxl, rl vector already partially calculated.
c         scale rl by r0nrm*prod to obtain the residual rl.
c
      tem = r0nrm*prod
      call dscal(n, tem, rl, 1)
      return
c------------- last line of drlcal follows ----------------------------
      end
