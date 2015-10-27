!
      SUBROUTINE DZTPSV ( UPLO, TRANS, DIAG, N, AP, X, INCX )
!     .. Scalar Arguments ..
      INTEGER                  :: INCX
      INTEGER                  :: N
      CHARACTER (LEN=1)        :: DIAG
      CHARACTER (LEN=1)        :: TRANS
      CHARACTER (LEN=1)        :: UPLO
!     .. Array Arguments ..
      REAL*8, DIMENSION(:)     :: AP
      COMPLEX*16, DIMENSION(:) ::  X
!     ..
!
!  Purpose
!  =======
!
!  DZTPSV  solves one of the systems of equations
!
!     A*x = b,   or   A'*x = b,
!
!  where b and x are n element vectors and A is an n by n unit, or
!  non-unit, upper or lower triangular matrix, supplied in packed form.
!
!  No test for singularity or near-singularity is included in this
!  routine. Such tests must be performed before calling this routine.
!
!  Parameters
!  ==========
!
!  UPLO   - CHARACTER.
!           On entry, UPLO specifies whether the matrix is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANS  - CHARACTER
!           On entry, TRANS specifies the equations to be solved as
!           follows:
!
!              TRANS = 'N' or 'n'   Ax = b.
!
!              TRANS = 'T' or 't'   A'x = b.
!
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER
!           On entry, DIAG specifies whether or not A is unit
!           triangular as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  AP     - DOUBLE PRECISION array of DIMENSION at least
!           ( ( n*( n + 1 ) )/2 ).
!           Before entry with  UPLO = 'U' or 'u', the array AP must
!           contain the upper triangular matrix packed sequentially,
!           column by column, so that AP( 1 ) contains a( 1, 1 ),
!           AP( 2 ) and AP( 3 ) contain a( 1, 2 ) and a( 2, 2 )
!           respectively, and so on.
!           Before entry with UPLO = 'L' or 'l', the array AP must
!           contain the lower triangular matrix packed sequentially,
!           column by column, so that AP( 1 ) contains a( 1, 1 ),
!           AP( 2 ) and AP( 3 ) contain a( 2, 1 ) and a( 3, 1 )
!           respectively, and so on.
!           Note that when  DIAG = 'U' or 'u', the diagonal elements of
!           A are not referenced, but are assumed to be unity.
!           Unchanged on exit.
!
!  X      - COMPLEX*16 array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element right-hand side vector b. On exit, X is overwritten
!           with the solution vector x.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      COMPLEX*16   ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
!     .. Local Scalars ..
      COMPLEX*16   TEMP
      INTEGER            I, INFO, IX, J, JX, K, KK, KX
      LOGICAL            NOUNIT
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF     ( .NOT.LSAME( UPLO , 'U' ).AND.        &
               .NOT.LSAME( UPLO , 'L' )      ) THEN
         INFO = 1
      ELSE IF( .NOT.LSAME( TRANS, 'N' ).AND.        &
               .NOT.LSAME( TRANS, 'T' ).AND.        &
               .NOT.LSAME( TRANS, 'C' )      ) THEN
         INFO = 2
      ELSE IF( .NOT.LSAME( DIAG , 'U' ).AND.        &
               .NOT.LSAME( DIAG , 'N' )      ) THEN
         INFO = 3
      ELSE IF( N.LT.0 )THEN
         INFO = 4
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 7
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DTPSV ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( N.EQ.0 ) THEN
          RETURN
      END IF
!
      NOUNIT = LSAME( DIAG, 'N' )
!
!     Set up the start point in X if the increment is not unity. This
!     will be  ( N - 1 )*INCX  too small for descending loops.
!
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
!
!     Start the operations. In this version the elements of AP are
!     accessed sequentially with one pass through AP.
!
      IF( LSAME( TRANS, 'N' ) )THEN
!
!        Form  x := inv( A )*x.
!
         IF( LSAME( UPLO, 'U' ) )THEN
            KK = ( N*( N + 1 ) )/2
            IF( INCX.EQ.1 )THEN
               DO J = N, 1, -1
                  IF( X( J ).NE.ZERO )THEN
                     IF( NOUNIT ) THEN
                         X( J ) = X( J )/AP( KK )
                     END IF
                     TEMP = X( J )
                     K    = KK     - 1
                     DO I = J - 1, 1, -1
                        X( I ) = X( I ) - TEMP*AP( K )
                        K      = K      - 1
                     END DO
                  END IF
                  KK = KK - J
               END DO
            ELSE
               JX = KX + ( N - 1 )*INCX
               DO J = N, 1, -1
                  IF( X( JX ).NE.ZERO )THEN
                     IF( NOUNIT ) THEN
                         X( JX ) = X( JX )/AP( KK )
                     END IF
                     TEMP = X( JX )
                     IX   = JX
                     DO K = KK - 1, KK - J + 1, -1
                        IX      = IX      - INCX
                        X( IX ) = X( IX ) - TEMP*AP( K )
                     END DO
                  END IF
                  JX = JX - INCX
                  KK = KK - J
               END DO
            END IF
         ELSE
            KK = 1
            IF( INCX.EQ.1 )THEN
               DO J = 1, N
                  IF( X( J ).NE.ZERO )THEN
                     IF( NOUNIT ) THEN
                         X( J ) = X( J )/AP( KK )
                     END IF
                     TEMP = X( J )
                     K    = KK     + 1
                     DO I = J + 1, N
                        X( I ) = X( I ) - TEMP*AP( K )
                        K      = K      + 1
                     END DO
                  END IF
                  KK = KK + ( N - J + 1 )
               END DO
            ELSE
               JX = KX
               DO J = 1, N
                  IF( X( JX ).NE.ZERO )THEN
                     IF( NOUNIT ) THEN
                         X( JX ) = X( JX )/AP( KK )
                     END IF
                     TEMP = X( JX )
                     IX   = JX
                     DO K = KK + 1, KK + N - J
                        IX      = IX      + INCX
                        X( IX ) = X( IX ) - TEMP*AP( K )
                     END DO
                  END IF
                  JX = JX + INCX
                  KK = KK + ( N - J + 1 )
   80          END DO
            END IF
         END IF
      ELSE
!
!        Form  x := inv( A' )*x.
!
         IF( LSAME( UPLO, 'U' ) )THEN
            KK = 1
            IF( INCX.EQ.1 )THEN
               DO J = 1, N
                  TEMP = X( J )
                  K    = KK
                  DO I = 1, J - 1
                     TEMP = TEMP - AP( K )*X( I )
                     K    = K    + 1
                  END DO
                  IF( NOUNIT ) THEN
                      TEMP = TEMP/AP( KK + J - 1 )
                  END IF
                  X( J ) = TEMP
                  KK     = KK   + J
               END DO
            ELSE
               JX = KX
               DO J = 1, N
                  TEMP = X( JX )
                  IX   = KX
                  DO K = KK, KK + J - 2
                     TEMP = TEMP - AP( K )*X( IX )
                     IX   = IX   + INCX
                  END DO
                  IF( NOUNIT ) THEN
                      TEMP = TEMP/AP( KK + J - 1 )
                  END IF
                  X( JX ) = TEMP
                  JX      = JX   + INCX
                  KK      = KK   + J
               END DO
            END IF
         ELSE
            KK = ( N*( N + 1 ) )/2
            IF( INCX.EQ.1 )THEN
               DO J = N, 1, -1
                  TEMP = X( J )
                  K = KK
                  DO I = N, J + 1, -1
                     TEMP = TEMP - AP( K )*X( I )
                     K    = K    - 1
                  END DO
                  IF( NOUNIT ) THEN
                      TEMP = TEMP/AP( KK - N + J )
                  END IF
                  X( J ) = TEMP
                  KK     = KK   - ( N - J + 1 )
               END DO
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO J = N, 1, -1
                  TEMP = X( JX )
                  IX   = KX
                  DO K = KK, KK - ( N - ( J + 1 ) ), -1
                     TEMP = TEMP - AP( K )*X( IX )
                     IX   = IX   - INCX
                  END DO
                  IF( NOUNIT ) THEN
                      TEMP = TEMP/AP( KK - N + J )
                  END IF
                  X( JX ) = TEMP
                  JX      = JX   - INCX
                  KK      = KK   - (N - J + 1 )
               END DO
            END IF
         END IF
      END IF
!
      RETURN
!
!
  END SUBROUTINE DZTPSV
