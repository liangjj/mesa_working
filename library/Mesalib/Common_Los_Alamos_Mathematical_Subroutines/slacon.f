      subroutine slacon( n, v, x, isgn, est, kase )
*
*  -- lapack auxiliary routine (version 1.1) --
*     univ. of tennessee, univ. of california berkeley, nag ltd.,
*     courant institute, argonne national lab, and rice university
*     february 29, 1992
*
*     .. scalar arguments ..
      integer            kase, n
      real*8             est
*     ..
*     .. array arguments ..
      integer            isgn( * )
      real*8             v( * ), x( * )
*     ..
*
*  purpose
*  =======
*
*  slacon estimates the 1-norm of a square, real matrix a.
*  reverse communication is used for evaluating matrix-vector products.
*
*  arguments
*  =========
*
*  n      (input) integer
*         the order of the matrix.  n >= 1.
*
*  v      (workspace) real array, dimension (n)
*         on the final return, v = a*w,  where  est = norm(v)/norm(w)
*         (w is not returned).
*
*  x      (input/output) real array, dimension (n)
*         on an intermediate return, x should be overwritten by
*               a * x,   if kase=1,
*               a' * x,  if kase=2,
*         and slacon must be re-called with all the other parameters
*         unchanged.
*
*  isgn   (workspace) integer array, dimension (n)
*
*  est    (output) real
*         an estimate (a lower bound) for norm(a).
*
*  kase   (input/output) integer
*         on the initial call to slacon, kase should be 0.
*         on an intermediate return, kase will be 1 or 2, indicating
*         whether x should be overwritten by a * x  or a' * x.
*         on the final return from slacon, kase will again be 0.
*
*  further details
*  ======= =======
*
*  contributed by nick higham, university of manchester.
*  originally named sonest, dated march 16, 1988.
*
*  reference: n.j. higham, "fortran codes for estimating the one-norm of
*  a real or complex matrix, with applications to condition estimation",
*  acm trans. math. soft., vol. 14, no. 4, pp. 381-396, december 1988.
*
*  =====================================================================
*
*     .. parameters ..
      integer            itmax
      parameter          ( itmax = 5 )
      real*8             zero, one, two
      parameter          ( zero = 0.0d+0, one = 1.0d+0, two = 2.0d+0 )
*     ..
*     .. local scalars ..
      integer            i, iter, j, jlast, jump
      real*8             altsgn, estold, temp
*     ..
*     .. external functions ..
      integer            isamax
      real*8             sasum
      external           isamax, sasum
*     ..
*     .. external subroutines ..
      external           scopy
*     ..
*     .. intrinsic functions ..
      intrinsic          abs, nint, real, sign
*     ..
*     .. save statement ..
      save
*     ..
*     .. executable statements ..
*
      if( kase.eq.0 ) then
         do 10 i = 1, n
            x( i ) = one / real( n )
   10    continue
         kase = 1
         jump = 1
         return
      end if
*
      go to ( 20, 40, 70, 110, 140 )jump
*
*     ................ entry   (jump = 1)
*     first iteration.  x has been overwritten by a*x.
*
   20 continue
      if( n.eq.1 ) then
         v( 1 ) = x( 1 )
         est = abs( v( 1 ) )
*        ... quit
         go to 150
      end if
      est = sasum( n, x, 1 )
*
      do 30 i = 1, n
         x( i ) = sign( one, x( i ) )
         isgn( i ) = nint( x( i ) )
   30 continue
      kase = 2
      jump = 2
      return
*
*     ................ entry   (jump = 2)
*     first iteration.  x has been overwritten by transpose(a)*x.
*
   40 continue
      j = isamax( n, x, 1 )
      iter = 2
*
*     main loop - iterations 2,3,...,itmax.
*
   50 continue
      do 60 i = 1, n
         x( i ) = zero
   60 continue
      x( j ) = one
      kase = 1
      jump = 3
      return
*
*     ................ entry   (jump = 3)
*     x has been overwritten by a*x.
*
   70 continue
      call scopy( n, x, 1, v, 1 )
      estold = est
      est = sasum( n, v, 1 )
      do 80 i = 1, n
         if( nint( sign( one, x( i ) ) ).ne.isgn( i ) )
     $      go to 90
   80 continue
*     repeated sign vector detected, hence algorithm has converged.
      go to 120
*
   90 continue
*     test for cycling.
      if( est.le.estold )
     $   go to 120
*
      do 100 i = 1, n
         x( i ) = sign( one, x( i ) )
         isgn( i ) = nint( x( i ) )
  100 continue
      kase = 2
      jump = 4
      return
*
*     ................ entry   (jump = 4)
*     x has been overwritten by transpose(a)*x.
*
  110 continue
      jlast = j
      j = isamax( n, x, 1 )
      if( ( x( jlast ).ne.abs( x( j ) ) ) .and. ( iter.lt.itmax ) ) then
         iter = iter + 1
         go to 50
      end if
*
*     iteration complete.  final stage.
*
  120 continue
      altsgn = one
      do 130 i = 1, n
         x( i ) = altsgn*( one+real( i-1 ) / real( n-1 ) )
         altsgn = -altsgn
  130 continue
      kase = 1
      jump = 5
      return
*
*     ................ entry   (jump = 5)
*     x has been overwritten by a*x.
*
  140 continue
      temp = two*( sasum( n, x, 1 ) / real( 3*n ) )
      if( temp.gt.est ) then
         call scopy( n, x, 1, v, 1 )
         est = temp
      end if
*
  150 continue
      kase = 0
      return
*
*     end of slacon
*
      end
