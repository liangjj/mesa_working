*deck dhels
      subroutine dhels (a, lda, n, q, b)
c***begin prologue  dhels
c***subsidiary
c***purpose  internal routine for dgmres.
c***library   slatec (slap)
c***category  d2a4, d2b4
c***type      double precision (shels-s, dhels-d)
c***keywords  generalized minimum residual, iterative precondition,
c             non-symmetric linear system, slap, sparse
c***author  brown, peter, (llnl), pnbrown@llnl.gov
c           hindmarsh, alan, (llnl), alanh@llnl.gov
c           seager, mark k., (llnl), seager@llnl.gov
c             lawrence livermore national laboratory
c             po box 808, l-60
c             livermore, ca 94550 (510) 423-3141
c***description
c        this routine is extracted from the linpack routine sgesl with
c        changes due to the fact that a is an upper hessenberg matrix.
c
c        dhels solves the least squares problem:
c
c                   min(b-a*x,b-a*x)
c
c        using the factors computed by dheqr.
c
c *usage:
c      integer lda, n
c      double precision a(lda,n), q(2*n), b(n+1)
c
c      call dhels(a, lda, n, q, b)
c
c *arguments:
c a       :in       double precision a(lda,n)
c          the output from dheqr which contains the upper
c          triangular factor r in the qr decomposition of a.
c lda     :in       integer
c          the leading dimension of the array a.
c n       :in       integer
c          a is originally an (n+1) by n matrix.
c q       :in       double precision q(2*n)
c          the coefficients of the n givens rotations
c          used in the qr factorization of a.
c b       :inout    double precision b(n+1)
c          on input, b is the right hand side vector.
c          on output, b is the solution vector x.
c
c***see also  dgmres
c***routines called  daxpy
c***revision history  (yymmdd)
c   890404  date written
c   890404  previous revision date
c   890915  made changes requested at july 1989 cml meeting.  (mks)
c   890922  numerous changes to prologue to make closer to slatec
c           standard.  (fnf)
c   890929  numerous changes to reduce sp/dp differences.  (fnf)
c   910411  prologue converted to version 4.0 format.  (bab)
c   910502  added c***first executable statement line.  (fnf)
c   910506  made subsidiary to dgmres.  (fnf)
c   920511  added complete declaration section.  (wrb)
c***end prologue  dhels
c         the following is for optimized compilation on llnl/ltss crays.
clll. optimize
c     .. scalar arguments ..
      integer lda, n
c     .. array arguments ..
      double precision a(lda,*), b(*), q(*)
c     .. local scalars ..
      double precision c, s, t, t1, t2
      integer iq, k, kb, kp1
c     .. external subroutines ..
      external daxpy
c***first executable statement  dhels
c
c         minimize(b-a*x,b-a*x).  first form q*b.
c
      do 20 k = 1, n
         kp1 = k + 1
         iq = 2*(k-1) + 1
         c = q(iq)
         s = q(iq+1)
         t1 = b(k)
         t2 = b(kp1)
         b(k) = c*t1 - s*t2
         b(kp1) = s*t1 + c*t2
 20   continue
c
c         now solve  r*x = q*b.
c
      do 40 kb = 1, n
         k = n + 1 - kb
         b(k) = b(k)/a(k,k)
         t = -b(k)
         call daxpy(k-1, t, a(1,k), 1, b(1), 1)
 40   continue
      return
c------------- last line of dhels follows ----------------------------
      end
