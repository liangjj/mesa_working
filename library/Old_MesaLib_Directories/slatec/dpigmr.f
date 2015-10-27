*deck dpigmr
      subroutine dpigmr (n, r0, sr, sz, jscal, maxl, maxlp1, kmp, nrsts,
     +   jpre, matvec, msolve, nmsl, z, v, hes, q, lgmr, rpar, ipar, wk,
     +   dl, rhol, nrmax, b, bnrm, x, xl, itol, tol, nelt, ia, ja, a,
     +   isym, iunit, iflag, err)
c***begin prologue  dpigmr
c***subsidiary
c***purpose  internal routine for dgmres.
c***library   slatec (slap)
c***category  d2a4, d2b4
c***type      double precision (spigmr-s, dpigmr-d)
c***keywords  generalized minimum residual, iterative precondition,
c             non-symmetric linear system, slap, sparse
c***author  brown, peter, (llnl), pnbrown@llnl.gov
c           hindmarsh, alan, (llnl), alanh@llnl.gov
c           seager, mark k., (llnl), seager@llnl.gov
c             lawrence livermore national laboratory
c             po box 808, l-60
c             livermore, ca 94550 (510) 423-3141
c***description
c         this routine solves the linear system a * z = r0 using a
c         scaled preconditioned version of the generalized minimum
c         residual method.  an initial guess of z = 0 is assumed.
c
c *usage:
c      integer n, jscal, maxl, maxlp1, kmp, nrsts, jpre, nmsl, lgmr
c      integer ipar(user defined), nrmax, itol, nelt, ia(nelt), ja(nelt)
c      integer isym, iunit, iflag
c      double precision r0(n), sr(n), sz(n), z(n), v(n,maxlp1),
c     $                 hes(maxlp1,maxl), q(2*maxl), rpar(user defined),
c     $                 wk(n), dl(n), rhol, b(n), bnrm, x(n), xl(n),
c     $                 tol, a(nelt), err
c      external matvec, msolve
c
c      call dpigmr(n, r0, sr, sz, jscal, maxl, maxlp1, kmp,
c     $     nrsts, jpre, matvec, msolve, nmsl, z, v, hes, q, lgmr,
c     $     rpar, ipar, wk, dl, rhol, nrmax, b, bnrm, x, xl,
c     $     itol, tol, nelt, ia, ja, a, isym, iunit, iflag, err)
c
c *arguments:
c n      :in       integer
c         the order of the matrix a, and the lengths
c         of the vectors sr, sz, r0 and z.
c r0     :in       double precision r0(n)
c         r0 = the right hand side of the system a*z = r0.
c         r0 is also used as workspace when computing
c         the final approximation.
c         (r0 is the same as v(*,maxl+1) in the call to dpigmr.)
c sr     :in       double precision sr(n)
c         sr is a vector of length n containing the non-zero
c         elements of the diagonal scaling matrix for r0.
c sz     :in       double precision sz(n)
c         sz is a vector of length n containing the non-zero
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
c maxl   :in       integer
c         the maximum allowable order of the matrix h.
c maxlp1 :in       integer
c         maxpl1 = maxl + 1, used for dynamic dimensioning of hes.
c kmp    :in       integer
c         the number of previous vectors the new vector vnew
c         must be made orthogonal to.  (kmp .le. maxl)
c nrsts  :in       integer
c         counter for the number of restarts on the current
c         call to dgmres.  if nrsts .gt. 0, then the residual
c         r0 is already scaled, and so scaling of it is
c         not necessary.
c jpre   :in       integer
c         preconditioner type flag.
c matvec :ext      external.
c         name of a routine which performs the matrix vector multiply
c         y = a*x given a and x.  the name of the matvec routine must
c         be declared external in the calling program.  the calling
c         sequence to matvec is:
c             call matvec(n, x, y, nelt, ia, ja, a, isym)
c         where n is the number of unknowns, y is the product a*x
c         upon return, x is an input vector, and nelt is the number of
c         non-zeros in the slap ia, ja, a storage for the matrix a.
c         isym is a flag which, if non-zero, denotes that a is
c         symmetric and only the lower or upper triangle is stored.
c msolve :ext      external.
c         name of the routine which solves a linear system mz = r for
c         z given r with the preconditioning matrix m (m is supplied via
c         rpar and ipar arrays.  the name of the msolve routine must
c         be declared external in the calling program.  the calling
c         sequence to msolve is:
c             call msolve(n, r, z, nelt, ia, ja, a, isym, rpar, ipar)
c         where n is the number of unknowns, r is the right-hand side
c         vector and z is the solution upon return.  nelt, ia, ja, a and
c         isym are defined as below.  rpar is a double precision array
c         that can be used to pass necessary preconditioning information
c         and/or workspace to msolve.  ipar is an integer work array
c         for the same purpose as rpar.
c nmsl   :out      integer
c         the number of calls to msolve.
c z      :out      double precision z(n)
c         the final computed approximation to the solution
c         of the system a*z = r0.
c v      :out      double precision v(n,maxlp1)
c         the n by (lgmr+1) array containing the lgmr
c         orthogonal vectors v(*,1) to v(*,lgmr).
c hes    :out      double precision hes(maxlp1,maxl)
c         the upper triangular factor of the qr decomposition
c         of the (lgmr+1) by lgmr upper hessenberg matrix whose
c         entries are the scaled inner-products of a*v(*,i)
c         and v(*,k).
c q      :out      double precision q(2*maxl)
c         a double precision array of length 2*maxl containing the
c         components of the givens rotations used in the qr
c         decomposition of hes.  it is loaded in dheqr and used in
c         dhels.
c lgmr   :out      integer
c         the number of iterations performed and
c         the current order of the upper hessenberg
c         matrix hes.
c rpar   :in       double precision rpar(user defined)
c         double precision workspace passed directly to the msolve
c         routine.
c ipar   :in       integer ipar(user defined)
c         integer workspace passed directly to the msolve routine.
c wk     :in       double precision wk(n)
c         a double precision work array of length n used by routines
c         matvec and msolve.
c dl     :inout    double precision dl(n)
c         on input, a double precision work array of length n used for
c         calculation of the residual norm rho when the method is
c         incomplete (kmp.lt.maxl), and/or when using restarting.
c         on output, the scaled residual vector rl.  it is only loaded
c         when performing restarts of the krylov iteration.
c rhol   :out      double precision
c         a double precision scalar containing the norm of the final
c         residual.
c nrmax  :in       integer
c         the maximum number of restarts of the krylov iteration.
c         nrmax .gt. 0 means restarting is active, while
c         nrmax = 0 means restarting is not being used.
c b      :in       double precision b(n)
c         the right hand side of the linear system a*x = b.
c bnrm   :in       double precision
c         the scaled norm of b.
c x      :in       double precision x(n)
c         the current approximate solution as of the last
c         restart.
c xl     :in       double precision xl(n)
c         an array of length n used to hold the approximate
c         solution x(l) when itol=11.
c itol   :in       integer
c         a flag to indicate the type of convergence criterion
c         used.  see the driver for its description.
c tol    :in       double precision
c         the tolerance on residuals r0-a*z in scaled norm.
c nelt   :in       integer
c         the length of arrays ia, ja and a.
c ia     :in       integer ia(nelt)
c         an integer array of length nelt containing matrix data.
c         it is passed directly to the matvec and msolve routines.
c ja     :in       integer ja(nelt)
c         an integer array of length nelt containing matrix data.
c         it is passed directly to the matvec and msolve routines.
c a      :in       double precision a(nelt)
c         a double precision array of length nelt containing matrix
c         data. it is passed directly to the matvec and msolve routines.
c isym   :in       integer
c         a flag to indicate symmetric matrix storage.
c         if isym=0, all non-zero entries of the matrix are
c         stored.  if isym=1, the matrix is symmetric and
c         only the upper or lower triangular part is stored.
c iunit  :in       integer
c         the i/o unit number for writing intermediate residual
c         norm values.
c iflag  :out      integer
c         an integer error flag..
c         0 means convergence in lgmr iterations, lgmr.le.maxl.
c         1 means the convergence test did not pass in maxl
c           iterations, but the residual norm is .lt. norm(r0),
c           and so z is computed.
c         2 means the convergence test did not pass in maxl
c           iterations, residual .ge. norm(r0), and z = 0.
c err    :out      double precision.
c         error estimate of error in final approximate solution, as
c         defined by itol.
c
c *cautions:
c     this routine will attempt to write to the fortran logical output
c     unit iunit, if iunit .ne. 0.  thus, the user must make sure that
c     this logical unit is attached to a file or terminal before calling
c     this routine with a non-zero value for iunit.  this routine does
c     not check for the validity of a non-zero iunit unit number.
c
c***see also  dgmres
c***routines called  daxpy, dcopy, dhels, dheqr, dnrm2, dorth, drlcal,
c                    dscal, isdgmr
c***revision history  (yymmdd)
c   890404  date written
c   890404  previous revision date
c   890915  made changes requested at july 1989 cml meeting.  (mks)
c   890922  numerous changes to prologue to make closer to slatec
c           standard.  (fnf)
c   890929  numerous changes to reduce sp/dp differences.  (fnf)
c   910411  prologue converted to version 4.0 format.  (bab)
c   910502  removed matvec and msolve from routines called list.  (fnf)
c   910506  made subsidiary to dgmres.  (fnf)
c   920511  added complete declaration section.  (wrb)
c***end prologue  dpigmr
c         the following is for optimized compilation on llnl/ltss crays.
clll. optimize
c     .. scalar arguments ..
      double precision bnrm, err, rhol, tol
      integer iflag, isym, itol, iunit, jpre, jscal, kmp, lgmr, maxl,
     +        maxlp1, n, nelt, nmsl, nrmax, nrsts
c     .. array arguments ..
      double precision a(nelt), b(*), dl(*), hes(maxlp1,*), q(*), r0(*),
     +                 rpar(*), sr(*), sz(*), v(n,*), wk(*), x(*),
     +                 xl(*), z(*)
      integer ia(nelt), ipar(*), ja(nelt)
c     .. subroutine arguments ..
      external matvec, msolve
c     .. local scalars ..
      double precision c, dlnrm, prod, r0nrm, rho, s, snormw, tem
      integer i, i2, info, ip1, iter, itmax, j, k, ll, llp1
c     .. external functions ..
      double precision dnrm2
      integer isdgmr
      external dnrm2, isdgmr
c     .. external subroutines ..
      external daxpy, dcopy, dhels, dheqr, dorth, drlcal, dscal
c     .. intrinsic functions ..
      intrinsic abs
c***first executable statement  dpigmr
c
c         zero out the z array.
c
      do 5 i = 1,n
         z(i) = 0
 5    continue
c
      iflag = 0
      lgmr = 0
      nmsl = 0
c         load itmax, the maximum number of iterations.
      itmax =(nrmax+1)*maxl
c   -------------------------------------------------------------------
c         the initial residual is the vector r0.
c         apply left precon. if jpre < 0 and this is not a restart.
c         apply scaling to r0 if jscal = 2 or 3.
c   -------------------------------------------------------------------
      if ((jpre .lt. 0) .and.(nrsts .eq. 0)) then
         call dcopy(n, r0, 1, wk, 1)
         call msolve(n, wk, r0, nelt, ia, ja, a, isym, rpar, ipar)
         nmsl = nmsl + 1
      endif
      if (((jscal.eq.2) .or.(jscal.eq.3)) .and.(nrsts.eq.0)) then
         do 10 i = 1,n
            v(i,1) = r0(i)*sr(i)
 10      continue
      else
         do 20 i = 1,n
            v(i,1) = r0(i)
 20      continue
      endif
      r0nrm = dnrm2(n, v, 1)
      iter = nrsts*maxl
c
c         call stopping routine isdgmr.
c
      if (isdgmr(n, b, x, xl, nelt, ia, ja, a, isym, msolve,
     $    nmsl, itol, tol, itmax, iter, err, iunit, v(1,1), z, wk,
     $    rpar, ipar, r0nrm, bnrm, sr, sz, jscal,
     $    kmp, lgmr, maxl, maxlp1, v, q, snormw, prod, r0nrm,
     $    hes, jpre) .ne. 0) return
      tem = 1.0d0/r0nrm
      call dscal(n, tem, v(1,1), 1)
c
c         zero out the hes array.
c
      do 50 j = 1,maxl
         do 40 i = 1,maxlp1
            hes(i,j) = 0
 40      continue
 50   continue
c   -------------------------------------------------------------------
c         main loop to compute the vectors v(*,2) to v(*,maxl).
c         the running product prod is needed for the convergence test.
c   -------------------------------------------------------------------
      prod = 1
      do 90 ll = 1,maxl
         lgmr = ll
c   -------------------------------------------------------------------
c        unscale  the  current v(ll)  and store  in wk.  call routine
c        msolve    to   compute(m-inverse)*wk,   where    m   is  the
c        preconditioner matrix.  save the answer in z.   call routine
c        matvec to compute  vnew  = a*z,  where  a is  the the system
c        matrix.  save the answer in  v(ll+1).  scale v(ll+1).   call
c        routine dorth  to  orthogonalize the    new vector vnew   =
c        v(*,ll+1).  call routine dheqr to update the factors of hes.
c   -------------------------------------------------------------------
        if ((jscal .eq. 1) .or.(jscal .eq. 3)) then
           do 60 i = 1,n
              wk(i) = v(i,ll)/sz(i)
 60        continue
        else
           call dcopy(n, v(1,ll), 1, wk, 1)
        endif
        if (jpre .gt. 0) then
           call msolve(n, wk, z, nelt, ia, ja, a, isym, rpar, ipar)
           nmsl = nmsl + 1
           call matvec(n, z, v(1,ll+1), nelt, ia, ja, a, isym)
        else
           call matvec(n, wk, v(1,ll+1), nelt, ia, ja, a, isym)
        endif
        if (jpre .lt. 0) then
           call dcopy(n, v(1,ll+1), 1, wk, 1)
           call msolve(n,wk,v(1,ll+1),nelt,ia,ja,a,isym,rpar,ipar)
           nmsl = nmsl + 1
        endif
        if ((jscal .eq. 2) .or.(jscal .eq. 3)) then
           do 65 i = 1,n
              v(i,ll+1) = v(i,ll+1)*sr(i)
 65        continue
        endif
        call dorth(v(1,ll+1), v, hes, n, ll, maxlp1, kmp, snormw)
        hes(ll+1,ll) = snormw
        call dheqr(hes, maxlp1, ll, q, info, ll)
        if (info .eq. ll) go to 120
c   -------------------------------------------------------------------
c         update rho, the estimate of the norm of the residual r0-a*zl.
c         if kmp <  maxl, then the vectors v(*,1),...,v(*,ll+1) are not
c         necessarily orthogonal for ll > kmp.  the vector dl must then
c         be computed, and its norm used in the calculation of rho.
c   -------------------------------------------------------------------
        prod = prod*q(2*ll)
        rho = abs(prod*r0nrm)
        if ((ll.gt.kmp) .and.(kmp.lt.maxl)) then
           if (ll .eq. kmp+1) then
              call dcopy(n, v(1,1), 1, dl, 1)
              do 75 i = 1,kmp
                 ip1 = i + 1
                 i2 = i*2
                 s = q(i2)
                 c = q(i2-1)
                 do 70 k = 1,n
                    dl(k) = s*dl(k) + c*v(k,ip1)
 70              continue
 75           continue
           endif
           s = q(2*ll)
           c = q(2*ll-1)/snormw
           llp1 = ll + 1
           do 80 k = 1,n
              dl(k) = s*dl(k) + c*v(k,llp1)
 80        continue
           dlnrm = dnrm2(n, dl, 1)
           rho = rho*dlnrm
        endif
        rhol = rho
c   -------------------------------------------------------------------
c         test for convergence.  if passed, compute approximation zl.
c         if failed and ll < maxl, then continue iterating.
c   -------------------------------------------------------------------
        iter = nrsts*maxl + lgmr
        if (isdgmr(n, b, x, xl, nelt, ia, ja, a, isym, msolve,
     $      nmsl, itol, tol, itmax, iter, err, iunit, dl, z, wk,
     $      rpar, ipar, rhol, bnrm, sr, sz, jscal,
     $      kmp, lgmr, maxl, maxlp1, v, q, snormw, prod, r0nrm,
     $      hes, jpre) .ne. 0) go to 200
        if (ll .eq. maxl) go to 100
c   -------------------------------------------------------------------
c         rescale so that the norm of v(1,ll+1) is one.
c   -------------------------------------------------------------------
        tem = 1.0d0/snormw
        call dscal(n, tem, v(1,ll+1), 1)
 90   continue
 100  continue
      if (rho .lt. r0nrm) go to 150
 120  continue
      iflag = 2
c
c         load approximate solution with zero.
c
      do 130 i = 1,n
         z(i) = 0
 130  continue
      return
 150  iflag = 1
c
c         tolerance not met, but residual norm reduced.
c
      if (nrmax .gt. 0) then
c
c        if performing restarting (nrmax > 0)  calculate the residual
c        vector rl and  store it in the dl  array.  if the incomplete
c        version is being used (kmp < maxl) then dl has  already been
c        calculated up to a scaling factor.   use drlcal to calculate
c        the scaled residual vector.
c
         call drlcal(n, kmp, maxl, maxl, v, q, dl, snormw, prod,
     $        r0nrm)
      endif
c   -------------------------------------------------------------------
c         compute the approximation zl to the solution.  since the
c         vector z was used as workspace, and the initial guess
c         of the linear iteration is zero, z must be reset to zero.
c   -------------------------------------------------------------------
 200  continue
      ll = lgmr
      llp1 = ll + 1
      do 210 k = 1,llp1
         r0(k) = 0
 210  continue
      r0(1) = r0nrm
      call dhels(hes, maxlp1, ll, q, r0)
      do 220 k = 1,n
         z(k) = 0
 220  continue
      do 230 i = 1,ll
         call daxpy(n, r0(i), v(1,i), 1, z, 1)
 230  continue
      if ((jscal .eq. 1) .or.(jscal .eq. 3)) then
         do 240 i = 1,n
            z(i) = z(i)/sz(i)
 240     continue
      endif
      if (jpre .gt. 0) then
         call dcopy(n, z, 1, wk, 1)
         call msolve(n, wk, z, nelt, ia, ja, a, isym, rpar, ipar)
         nmsl = nmsl + 1
      endif
      return
c------------- last line of dpigmr follows ----------------------------
      end
