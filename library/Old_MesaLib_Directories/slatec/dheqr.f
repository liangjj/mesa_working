*deck dheqr
      subroutine dheqr (a, lda, n, q, info, ijob)
c***begin prologue  dheqr
c***subsidiary
c***purpose  internal routine for dgmres.
c***library   slatec (slap)
c***category  d2a4, d2b4
c***type      double precision (sheqr-s, dheqr-d)
c***keywords  generalized minimum residual, iterative precondition,
c             non-symmetric linear system, slap, sparse
c***author  brown, peter, (llnl), pnbrown@llnl.gov
c           hindmarsh, alan, (llnl), alanh@llnl.gov
c           seager, mark k., (llnl), seager@llnl.gov
c             lawrence livermore national laboratory
c             po box 808, l-60
c             livermore, ca 94550 (510) 423-3141
c***description
c        this   routine  performs  a qr   decomposition  of an  upper
c        hessenberg matrix a using givens  rotations.  there  are two
c        options  available: 1)  performing  a fresh decomposition 2)
c        updating the qr factors by adding a row and  a column to the
c        matrix a.
c
c *usage:
c      integer lda, n, info, ijob
c      double precision a(lda,n), q(2*n)
c
c      call dheqr(a, lda, n, q, info, ijob)
c
c *arguments:
c a      :inout    double precision a(lda,n)
c         on input, the matrix to be decomposed.
c         on output, the upper triangular matrix r.
c         the factorization can be written q*a = r, where
c         q is a product of givens rotations and r is upper
c         triangular.
c lda    :in       integer
c         the leading dimension of the array a.
c n      :in       integer
c         a is an (n+1) by n hessenberg matrix.
c q      :out      double precision q(2*n)
c         the factors c and s of each givens rotation used
c         in decomposing a.
c info   :out      integer
c         = 0  normal value.
c         = k  if  a(k,k) .eq. 0.0 .  this is not an error
c           condition for this subroutine, but it does
c           indicate that dhels will divide by zero
c           if called.
c ijob   :in       integer
c         = 1     means that a fresh decomposition of the
c                 matrix a is desired.
c         .ge. 2  means that the current decomposition of a
c                 will be updated by the addition of a row
c                 and a column.
c
c***see also  dgmres
c***routines called  (none)
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
c***end prologue  dheqr
c         the following is for optimized compilation on llnl/ltss crays.
clll. optimize
c     .. scalar arguments ..
      integer ijob, info, lda, n
c     .. array arguments ..
      double precision a(lda,*), q(*)
c     .. local scalars ..
      double precision c, s, t, t1, t2
      integer i, iq, j, k, km1, kp1, nm1
c     .. intrinsic functions ..
      intrinsic abs, sqrt
c***first executable statement  dheqr
      if (ijob .gt. 1) go to 70
c   -------------------------------------------------------------------
c         a new factorization is desired.
c   -------------------------------------------------------------------
c         qr decomposition without pivoting.
c
      info = 0
      do 60 k = 1, n
         km1 = k - 1
         kp1 = k + 1
c
c           compute k-th column of r.
c           first, multiply the k-th column of a by the previous
c           k-1 givens rotations.
c
         if (km1 .lt. 1) go to 20
         do 10 j = 1, km1
            i = 2*(j-1) + 1
            t1 = a(j,k)
            t2 = a(j+1,k)
            c = q(i)
            s = q(i+1)
            a(j,k) = c*t1 - s*t2
            a(j+1,k) = s*t1 + c*t2
 10      continue
c
c         compute givens components c and s.
c
 20      continue
         iq = 2*km1 + 1
         t1 = a(k,k)
         t2 = a(kp1,k)
         if( t2.eq.0.0d0 ) then
            c = 1
            s = 0
         elseif( abs(t2).ge.abs(t1) ) then
            t = t1/t2
            s = -1.0d0/sqrt(1.0d0+t*t)
            c = -s*t
         else
            t = t2/t1
            c = 1.0d0/sqrt(1.0d0+t*t)
            s = -c*t
         endif
         q(iq) = c
         q(iq+1) = s
         a(k,k) = c*t1 - s*t2
         if( a(k,k).eq.0.0d0 ) info = k
 60   continue
      return
c   -------------------------------------------------------------------
c         the old factorization of a will be updated.  a row and a
c         column has been added to the matrix a.  n by n-1 is now
c         the old size of the matrix.
c   -------------------------------------------------------------------
 70   continue
      nm1 = n - 1
c   -------------------------------------------------------------------
c         multiply the new column by the n previous givens rotations.
c   -------------------------------------------------------------------
      do 100 k = 1,nm1
         i = 2*(k-1) + 1
         t1 = a(k,n)
         t2 = a(k+1,n)
         c = q(i)
         s = q(i+1)
         a(k,n) = c*t1 - s*t2
         a(k+1,n) = s*t1 + c*t2
 100  continue
c   -------------------------------------------------------------------
c         complete update of decomposition by forming last givens
c         rotation, and multiplying it times the column
c         vector(a(n,n),a(np1,n)).
c   -------------------------------------------------------------------
      info = 0
      t1 = a(n,n)
      t2 = a(n+1,n)
      if ( t2.eq.0.0d0 ) then
         c = 1
         s = 0
      elseif( abs(t2).ge.abs(t1) ) then
         t = t1/t2
         s = -1.0d0/sqrt(1.0d0+t*t)
         c = -s*t
      else
         t = t2/t1
         c = 1.0d0/sqrt(1.0d0+t*t)
         s = -c*t
      endif
      iq = 2*n - 1
      q(iq) = c
      q(iq+1) = s
      a(n,n) = c*t1 - s*t2
      if (a(n,n) .eq. 0.0d0) info = n
      return
c------------- last line of dheqr follows ----------------------------
      end
