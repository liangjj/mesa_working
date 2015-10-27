*deck sxlcal
      subroutine sxlcal (n, lgmr, x, xl, zl, hes, maxlp1, q, v, r0nrm,
     +   wk, sz, jscal, jpre, msolve, nmsl, rpar, ipar, nelt, ia, ja, a,
     +   isym)
c***begin prologue  sxlcal
c***subsidiary
c***purpose  internal routine for sgmres.
c***library   slatec (slap)
c***category  d2a4, d2b4
c***type      single precision (sxlcal-s, dxlcal-d)
c***keywords  generalized minimum residual, iterative precondition,
c             non-symmetric linear system, slap, sparse
c***author  brown, peter, (llnl), pnbrown@llnl.gov
c           hindmarsh, alan, (llnl), alanh@llnl.gov
c           seager, mark k., (llnl), seager@llnl.gov
c             lawrence livermore national laboratory
c             po box 808, l-60
c             livermore, ca 94550 (510) 423-3141
c***description
c        this  routine computes the solution  xl,  the current sgmres
c        iterate, given the  v(i)'s and  the  qr factorization of the
c        hessenberg  matrix hes.   this routine  is  only called when
c        itol=11.
c
c *usage:
c      integer n, lgmr, maxlp1, jscal, jpre, nmsl, ipar(user defined)
c      integer nelt, ia(nelt), ja(nelt), isym
c      real x(n), xl(n), zl(n), hes(maxlp1,maxl), q(2*maxl),
c     $     v(n,maxlp1), r0nrm, wk(n), sz(n), rpar(user defined),
c     $     a(nelt)
c      external msolve
c
c      call sxlcal(n, lgmr, x, xl, zl, hes, maxlp1, q, v, r0nrm,
c     $     wk, sz, jscal, jpre, msolve, nmsl, rpar, ipar,
c     $     nelt, ia, ja, a, isym)
c
c *arguments:
c n      :in       integer
c         the order of the matrix a, and the lengths
c         of the vectors sr, sz, r0 and z.
c lgmr   :in       integer
c         the number of iterations performed and
c         the current order of the upper hessenberg
c         matrix hes.
c x      :in       real x(n)
c         the current approximate solution as of the last restart.
c xl     :out      real xl(n)
c         an array of length n used to hold the approximate
c         solution x(l).
c         warning: xl and zl are the same array in the calling routine.
c zl     :in       real zl(n)
c         an array of length n used to hold the approximate
c         solution z(l).
c hes    :in       real hes(maxlp1,maxl)
c         the upper triangular factor of the qr decomposition
c         of the (lgmr+1) by lgmr upper hessenberg matrix whose
c         entries are the scaled inner-products of a*v(*,i) and v(*,k).
c maxlp1 :in       integer
c         maxlp1 = maxl + 1, used for dynamic dimensioning of hes.
c         maxl is the maximum allowable order of the matrix hes.
c q      :in       real q(2*maxl)
c         a real array of length 2*maxl containing the components
c         of the givens rotations used in the qr decomposition
c         of hes.  it is loaded in sheqr.
c v      :in       real v(n,maxlp1)
c         the n by(lgmr+1) array containing the lgmr
c         orthogonal vectors v(*,1) to v(*,lgmr).
c r0nrm  :in       real
c         the scaled norm of the initial residual for the
c         current call to spigmr.
c wk     :in       real wk(n)
c         a real work array of length n.
c sz     :in       real sz(n)
c         a vector of length n containing the non-zero
c         elements of the diagonal scaling matrix for z.
c jscal  :in       integer
c         a flag indicating whether arrays sr and sz are used.
c         jscal=0 means sr and sz are not used and the
c                 algorithm will perform as if all
c                 sr(i) = 1 and sz(i) = 1.
c         jscal=1 means only sz is used, and the algorithm
c                 performs as if all sr(i) = 1.
c         jscal=2 means only sr is used, and the algorithm
c                 performs as if all sz(i) = 1.
c         jscal=3 means both sr and sz are used.
c jpre   :in       integer
c         the preconditioner type flag.
c msolve :ext      external.
c         name of the routine which solves a linear system mz = r for
c         z given r with the preconditioning matrix m (m is supplied via
c         rpar and ipar arrays.  the name of the msolve routine must
c         be declared external in the calling program.  the calling
c         sequence to msolve is:
c             call msolve(n, r, z, nelt, ia, ja, a, isym, rpar, ipar)
c         where n is the number of unknowns, r is the right-hand side
c         vector and z is the solution upon return.  nelt, ia, ja, a and
c         isym are defined as below.  rpar is a real array that can be
c         used to pass necessary preconditioning information and/or
c         workspace to msolve.  ipar is an integer work array for the
c         same purpose as rpar.
c nmsl   :in       integer
c         the number of calls to msolve.
c rpar   :in       real rpar(user defined)
c         real workspace passed directly to the msolve routine.
c ipar   :in       integer ipar(user defined)
c         integer workspace passed directly to the msolve routine.
c nelt   :in       integer
c         the length of arrays ia, ja and a.
c ia     :in       integer ia(nelt)
c         an integer array of length nelt containing matrix data.
c         it is passed directly to the matvec and msolve routines.
c ja     :in       integer ja(nelt)
c         an integer array of length nelt containing matrix data.
c         it is passed directly to the matvec and msolve routines.
c a      :in       real a(nelt)
c         a real array of length nelt containing matrix data.
c         it is passed directly to the matvec and msolve routines.
c isym   :in       integer
c         a flag to indicate symmetric matrix storage.
c         if isym=0, all non-zero entries of the matrix are
c         stored.  if isym=1, the matrix is symmetric and
c         only the upper or lower triangular part is stored.
c
c***see also  sgmres
c***routines called  saxpy, scopy, shels
c***revision history  (yymmdd)
c   871001  date written
c   881213  previous revision date
c   890915  made changes requested at july 1989 cml meeting.  (mks)
c   890922  numerous changes to prologue to make closer to slatec
c           standard.  (fnf)
c   890929  numerous changes to reduce sp/dp differences.  (fnf)
c   910411  prologue converted to version 4.0 format.  (bab)
c   910502  removed msolve from routines called list.  (fnf)
c   910506  made subsidiary to sgmres.  (fnf)
c   920511  added complete declaration section.  (wrb)
c***end prologue  sxlcal
c         the following is for optimized compilation on llnl/ltss crays.
clll. optimize
c     .. scalar arguments ..
      real r0nrm
      integer isym, jpre, jscal, lgmr, maxlp1, n, nelt, nmsl
c     .. array arguments ..
      real a(nelt), hes(maxlp1,*), q(*), rpar(*), sz(*), v(n,*), wk(n),
     +     x(n), xl(n), zl(n)
      integer ia(nelt), ipar(*), ja(nelt)
c     .. subroutine arguments ..
      external msolve
c     .. local scalars ..
      integer i, k, ll, llp1
c     .. external subroutines ..
      external saxpy, scopy, shels
c***first executable statement  sxlcal
      ll = lgmr
      llp1 = ll + 1
      do 10 k = 1,llp1
         wk(k) = 0
 10   continue
      wk(1) = r0nrm
      call shels(hes, maxlp1, ll, q, wk)
      do 20 k = 1,n
         zl(k) = 0
 20   continue
      do 30 i = 1,ll
         call saxpy(n, wk(i), v(1,i), 1, zl, 1)
 30   continue
      if ((jscal .eq. 1) .or.(jscal .eq. 3)) then
         do 40 k = 1,n
            zl(k) = zl(k)/sz(k)
 40      continue
      endif
      if (jpre .gt. 0) then
         call scopy(n, zl, 1, wk, 1)
         call msolve(n, wk, zl, nelt, ia, ja, a, isym, rpar, ipar)
         nmsl = nmsl + 1
      endif
c         calculate xl from x and zl.
      do 50 k = 1,n
         xl(k) = x(k) + zl(k)
 50   continue
      return
c------------- last line of sxlcal follows ----------------------------
      end
