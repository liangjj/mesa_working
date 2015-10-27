      SUBROUTINE ZLARGV( N, X, INCX, Y, INCY, C, INCC )
*
*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      INTEGER            INCC, INCX, INCY, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   C( * )
      COMPLEX*16         X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  ZLARGV generates a vector of complex plane rotations with real
*  cosines, determined by elements of the complex vectors x and y.
*  For i = 1,2,...,n
*
*     (        c(i)   s(i) ) ( x(i) ) = ( a(i) )
*     ( -conjg(s(i))  c(i) ) ( y(i) ) = (   0  )
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of plane rotations to be generated.
*
*  X       (input/output) COMPLEX*16 array, dimension (1+(N-1)*INCX)
*          On entry, the vector x.
*          On exit, x(i) is overwritten by a(i), for i = 1,...,n.
*
*  INCX    (input) INTEGER
*          The increment between elements of X. INCX > 0.
*
*  Y       (input/output) COMPLEX*16 array, dimension (1+(N-1)*INCY)
*          On entry, the vector y.
*          On exit, the sines of the plane rotations.
*
*  INCY    (input) INTEGER
*          The increment between elements of Y. INCY > 0.
*
*  C       (output) DOUBLE PRECISION array, dimension (1+(N-1)*INCC)
*          The cosines of the plane rotations.
*
*  INCC    (input) INTEGER
*          The increment between elements of C. INCC > 0.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IC, IX, IY
      DOUBLE PRECISION   ABSX, ABSY, TT, W
      COMPLEX*16         T, XI, YI
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DCONJG, MAX, SQRT
*     ..
*     .. Executable Statements ..
*
      IX = 1
      IY = 1
      IC = 1
      DO 10 I = 1, N
         XI = X( IX )
         YI = Y( IY )
         ABSX = ABS( XI )
         IF( ABSX.EQ.ZERO ) THEN
            C( IC ) = ZERO
            Y( IY ) = ONE
            X( IX ) = YI
         ELSE
            ABSY = ABS( YI )
            W = MAX( ABSX, ABSY )
            T = XI / ABSX
            ABSX = ABSX / W
            ABSY = ABSY / W
            TT = SQRT( ABSX*ABSX+ABSY*ABSY )
            C( IC ) = ABSX / TT
            Y( IY ) = ( T*DCONJG( YI ) ) / ( W*TT )
            X( IX ) = T*( W*TT )
         END IF
         IX = IX + INCX
         IY = IY + INCY
         IC = IC + INCC
   10 CONTINUE
      RETURN
*
*     End of ZLARGV
*
      END
