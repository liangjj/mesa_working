      SUBROUTINE ZSGT01( ITYPE, UPLO, N, A, LDA, B, LDB, Z, LDZ, D,
     $                   WORK, RWORK, RESULT )
*
*  -- LAPACK test routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            ITYPE, LDA, LDB, LDZ, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), RESULT( * ), RWORK( * )
      COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * ),
     $                   Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  CDGT01 checks a decomposition of the form
*
*     A =  B Z D Z' or
*     A B =  Z D Z' or
*     B A =  Z D Z'
*
*  where ' means conjugate transpose, A is a Hermitian matrix, B is
*  Hermitian positive definite, Z is orthogonal, and D is diagonal.
*
*  One of the following test ratios is computed:
*
*  ITYPE = 1:  RESULT(1) = | A - B Z D Z' | / ( |A| n ulp )
*
*  ITYPE = 2:  RESULT(1) = | A B - Z D Z' | / ( |A| n ulp )
*
*  ITYPE = 3:  RESULT(1) = | B A - Z D Z' | / ( |A| n ulp )
*
*  Arguments
*  =========
*
*  ITYPE   (input) INTEGER
*          The form of the Hermitian generalized eigenproblem.
*          = 1:  A*z = (lambda)*B*z
*          = 2:  A*B*z = (lambda)*z
*          = 3:  B*A*z = (lambda)*z
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          Hermitian matrices A and B is stored.
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input) COMPLEX*16 array, dimension (LDA, N)
*          The original Hermitian matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  B       (input) COMPLEX*16 array, dimension (LDB, N)
*          The original Hermitian positive definite matrix B.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  Z       (input) COMPLEX*16 array, dimension (LDZ, N)
*          The eigenvectors of the generalized eigenproblem, as computed
*          by ZHEGV.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= max(1,N).
*
*  D       (input) DOUBLE PRECISION array, dimension (N)
*          The eigenvalues of the generalized eigenproblem, as computed
*          by ZHEGV.
*
*  WORK    (workspace) COMPLEX*16 array, dimension (N*N)
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension (N)
*
*  RESULT  (output) DOUBLE PRECISION array, dimension (1)
*          The test ratio as described above.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = 0.0D0, CONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   ANORM, ULP
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, ZLANGE, ZLANSY
      EXTERNAL           DLAMCH, ZLANGE, ZLANSY
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZDSCAL, ZHEMM
*     ..
*     .. Executable Statements ..
*
      RESULT( 1 ) = ZERO
      IF( N.LE.0 )
     $   RETURN
*
      ULP = DLAMCH( 'Epsilon' )
*
*     Compute 1-norm of A.
*
      ANORM = ZLANSY( '1', UPLO, N, A, LDA, RWORK )*
     $        ZLANGE( '1', N, N, Z, LDZ, RWORK )
      IF( ANORM.EQ.ZERO )
     $   ANORM = ONE
*
      IF( ITYPE.EQ.1 ) THEN
*
*        Norm of AZ - BZD
*
         CALL ZHEMM( 'Left', UPLO, N, N, CONE, A, LDA, Z, LDZ, CZERO,
     $               WORK, N )
         DO 10 I = 1, N
            CALL ZDSCAL( N, D( I ), Z( 1, I ), 1 )
   10    CONTINUE
*
         CALL ZHEMM( 'Left', UPLO, N, N, CONE, B, LDB, Z, LDZ, -CONE,
     $               WORK, N )
*
         RESULT( 1 ) = ( ZLANGE( '1', N, N, WORK, N, RWORK ) / ANORM ) /
     $                 ( N*ULP )
*
      ELSE IF( ITYPE.EQ.2 ) THEN
*
*        Norm of ABZ - ZD
*
         CALL ZHEMM( 'Left', UPLO, N, N, CONE, B, LDB, Z, LDZ, CZERO,
     $               WORK, N )
         DO 20 I = 1, N
            CALL ZDSCAL( N, D( I ), Z( 1, I ), 1 )
   20    CONTINUE
         CALL ZHEMM( 'Left', UPLO, N, N, CONE, A, LDA, WORK, N, -CONE,
     $               Z, LDZ )
*
         RESULT( 1 ) = ( ZLANGE( '1', N, N, Z, LDZ, RWORK ) / ANORM ) /
     $                 ( N*ULP )
*
      ELSE IF( ITYPE.EQ.3 ) THEN
*
*        Norm of BAZ - ZD
*
         CALL ZHEMM( 'Left', UPLO, N, N, CONE, A, LDA, Z, LDZ, CZERO,
     $               WORK, N )
         DO 30 I = 1, N
            CALL ZDSCAL( N, D( I ), Z( 1, I ), 1 )
   30    CONTINUE
         CALL ZHEMM( 'Left', UPLO, N, N, CONE, B, LDB, WORK, N, -CONE,
     $               Z, LDZ )
*
         RESULT( 1 ) = ( ZLANGE( '1', N, N, Z, LDZ, RWORK ) / ANORM ) /
     $                 ( N*ULP )
      END IF
*
      RETURN
*
*     End of CDGT01
*
      END
