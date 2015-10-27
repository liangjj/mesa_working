      SUBROUTINE ZGBSVX( FACT, TRANS, N, KL, KU, NRHS, AB, LDAB, AFB,
     $                   LDAFB, IPIV, EQUED, R, C, B, LDB, X, LDX,
     $                   RCOND, FERR, BERR, WORK, RWORK, INFO )
*
*  -- LAPACK driver routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      CHARACTER          EQUED, FACT, TRANS
      INTEGER            INFO, KL, KU, LDAB, LDAFB, LDB, LDX, N, NRHS
      DOUBLE PRECISION   RCOND
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   BERR( * ), C( * ), FERR( * ), R( * ),
     $                   RWORK( * )
      COMPLEX*16         AB( LDAB, * ), AFB( LDAFB, * ), B( LDB, * ),
     $                   WORK( * ), X( LDX, * )
*     ..
*
*  Purpose
*  =======
*
*  ZGBSVX uses the LU factorization to compute the solution to a complex
*  system of linear equations A * X = B, A**T * X = B, or A**H * X = B,
*  where A is a band matrix of order N with KL subdiagonals and KU
*  superdiagonals, and X and B are N-by-NRHS matrices.
*
*  Error bounds on the solution and a condition estimate are also
*  provided.
*
*  Description
*  ===========
*
*  The following steps are performed by this subroutine:
*
*  1. If FACT = 'E', real scaling factors are computed to equilibrate
*     the system:
*        TRANS = 'N':  diag(R)*A*diag(C)     *inv(diag(C))*X = diag(R)*B
*        TRANS = 'T': (diag(R)*A*diag(C))**T *inv(diag(R))*X = diag(C)*B
*        TRANS = 'C': (diag(R)*A*diag(C))**H *inv(diag(R))*X = diag(C)*B
*     Whether or not the system will be equilibrated depends on the
*     scaling of the matrix A, but if equilibration is used, A is
*     overwritten by diag(R)*A*diag(C) and B by diag(R)*B (if TRANS='N')
*     or diag(C)*B (if TRANS = 'T' or 'C').
*
*  2. If FACT = 'N' or 'E', the LU decomposition is used to factor the
*     matrix A (after equilibration if FACT = 'E') as
*        A = L * U,
*     where L is a product of permutation and unit lower triangular
*     matrices with KL subdiagonals, and U is upper triangular with
*     KL+KU superdiagonals.
*
*  3. The factored form of A is used to estimate the condition number
*     of the matrix A.  If the reciprocal of the condition number is
*     less than machine precision, steps 4-6 are skipped.
*
*  4. The system of equations is solved for X using the factored form
*     of A.
*
*  5. Iterative refinement is applied to improve the computed solution
*     matrix and calculate error bounds and backward error estimates
*     for it.
*
*  6. If equilibration was used, the matrix X is premultiplied by
*     diag(C) (if TRANS = 'N') or diag(R) (if TRANS = 'T' or 'C') so
*     that it solves the original system before equilibration.
*
*  Arguments
*  =========
*
*  FACT    (input) CHARACTER*1
*          Specifies whether or not the factored form of the matrix A is
*          supplied on entry, and if not, whether the matrix A should be
*          equilibrated before it is factored.
*          = 'F':  On entry, AFB and IPIV contain the factored form of
*                  A.  If EQUED is not 'N', the matrix A has been
*                  equilibrated with scaling factors given by R and C.
*                  AB, AFB, and IPIV are not modified.
*          = 'N':  The matrix A will be copied to AFB and factored.
*          = 'E':  The matrix A will be equilibrated if necessary, then
*                  copied to AFB and factored.
*
*  TRANS   (input) CHARACTER*1
*          Specifies the form of the system of equations.
*          = 'N':  A * X = B     (No transpose)
*          = 'T':  A**T * X = B  (Transpose)
*          = 'C':  A**H * X = B  (Conjugate transpose)
*
*  N       (input) INTEGER
*          The number of linear equations, i.e., the order of the
*          matrix A.  N >= 0.
*
*  KL      (input) INTEGER
*          The number of subdiagonals within the band of A.  KL >= 0.
*
*  KU      (input) INTEGER
*          The number of superdiagonals within the band of A.  KU >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right-hand sides, i.e., the number of columns
*          of the matrices B and X.  NRHS >= 0.
*
*  AB      (input/output) COMPLEX*16 array, dimension (LDAB,N)
*          On entry, the matrix A in band storage, in rows 1 to KL+KU+1.
*          The j-th column of A is stored in the j-th column of the
*          array AB as follows:
*          AB(KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+kl)
*
*          If FACT = 'F' and EQUED is not 'N', then A must have been
*          equilibrated by the scaling factors in R and/or C.  AB is not
*          modified if FACT = 'F' or 'N', or if FACT = 'E' and
*          EQUED = 'N' on exit.
*
*          On exit, if EQUED .ne. 'N', A is scaled as follows:
*          EQUED = 'R':  A := diag(R) * A
*          EQUED = 'C':  A := A * diag(C)
*          EQUED = 'B':  A := diag(R) * A * diag(C).
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= KL+KU+1.
*
*  AFB     (input or output) COMPLEX*16 array, dimension (LDAFB,N)
*          If FACT = 'F', then AFB is an input argument and on entry
*          contains details of the LU factorization of the band matrix
*          A, as computed by ZGBTRF.  U is stored as an upper triangular
*          band matrix with KL+KU superdiagonals in rows 1 to KL+KU+1,
*          and the multipliers used during the factorization are stored
*          in rows KL+KU+2 to 2*KL+KU+1.  If EQUED .ne. 'N', then AFB is
*          the factored form of the equilibrated matrix A.
*
*          If FACT = 'N', then AFB is an output argument and on exit
*          returns details of the LU factorization of A.
*
*          If FACT = 'E', then AFB is an output argument and on exit
*          returns details of the LU factorization of the equilibrated
*          matrix A (see the description of AB for the form of the
*          equilibrated matrix).
*
*  LDAFB   (input) INTEGER
*          The leading dimension of the array AFB.  LDAFB >= 2*KL+KU+1.
*
*  IPIV    (input or output) INTEGER array, dimension (N)
*          If FACT = 'F', then IPIV is an input argument and on entry
*          contains the pivot indices from the factorization A = L*U
*          as computed by ZGBTRF; row i of the matrix was interchanged
*          with row IPIV(i).
*
*          If FACT = 'N', then IPIV is an output argument and on exit
*          contains the pivot indices from the factorization A = L*U
*          of the original matrix A.
*
*          If FACT = 'E', then IPIV is an output argument and on exit
*          contains the pivot indices from the factorization A = L*U
*          of the equilibrated matrix A.
*
*  EQUED   (input/output) CHARACTER*1
*          Specifies the form of equilibration that was done.
*          = 'N':  No equilibration (always true if FACT = 'N').
*          = 'R':  Row equilibration, i.e., A has been premultiplied by
*                  diag(R).
*          = 'C':  Column equilibration, i.e., A has been postmultiplied
*                  by diag(C).
*          = 'B':  Both row and column equilibration, i.e., A has been
*                  replaced by diag(R) * A * diag(C).
*          EQUED is an input variable if FACT = 'F'; otherwise, it is an
*          output variable.
*
*  R       (input/output) DOUBLE PRECISION array, dimension (N)
*          The row scale factors for A.  If EQUED = 'R' or 'B', A is
*          multiplied on the left by diag(R); if EQUED = 'N' or 'C', R
*          is not accessed.  R is an input variable if FACT = 'F';
*          otherwise, R is an output variable.  If FACT = 'F' and
*          EQUED = 'R' or 'B', each element of R must be positive.
*
*  C       (input/output) DOUBLE PRECISION array, dimension (N)
*          The column scale factors for A.  If EQUED = 'C' or 'B', A is
*          multiplied on the right by diag(C); if EQUED = 'N' or 'R', C
*          is not accessed.  C is an input variable if FACT = 'F';
*          otherwise, C is an output variable.  If FACT = 'F' and
*          EQUED = 'C' or 'B', each element of C must be positive.
*
*  B       (input/output) COMPLEX*16 array, dimension (LDB,NRHS)
*          On entry, the right-hand side matrix B.
*          On exit, if EQUED = 'N', B is not modified; if TRANS = 'N'
*          and EQUED = 'R' or 'B', B is overwritten by diag(R)*B; if
*          TRANS = 'T' or 'C' and EQUED = 'C' or 'B', B is overwritten
*          by diag(C)*B.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  X       (output) COMPLEX*16 array, dimension (LDX,NRHS)
*          If INFO = 0, the n-by-nrhs solution matrix X to the original
*          system of equations.  Note that A and B are modified on exit
*          if EQUED .ne. 'N', and the solution to the equilibrated
*          system is inv(diag(C))*X if TRANS = 'N' and EQUED = 'C' or
*          'B', or inv(diag(R))*X if TRANS = 'T' or 'C' and EQUED = 'R'
*          or 'B'.
*
*  LDX     (input) INTEGER
*          The leading dimension of the array X.  LDX >= max(1,N).
*
*  RCOND   (output) DOUBLE PRECISION
*          The estimate of the reciprocal condition number of the matrix
*          A after equilibration (if done).  If RCOND is less than the
*          machine precision (in particular, if RCOND = 0), the matrix
*          is singular to working precision.  This condition is
*          indicated by a return code of INFO > 0, and the solution and
*          error bounds are not computed.
*
*  FERR    (output) DOUBLE PRECISION array, dimension (NRHS)
*          The estimated forward error bounds for each solution vector
*          X(j) (the j-th column of the solution matrix X).
*          If XTRUE is the true solution, FERR(j) bounds the magnitude
*          of the largest entry in (X(j) - XTRUE) divided by the
*          magnitude of the largest entry in X(j).  The quality of the
*          error bound depends on the quality of the estimate of
*          norm(inv(A)) computed in the code; if the estimate of
*          norm(inv(A)) is accurate, the error bound is guaranteed.
*
*  BERR    (output) DOUBLE PRECISION array, dimension (NRHS)
*          The componentwise relative backward error of each solution
*          vector X(j) (i.e., the smallest relative change in
*          any entry of A or B that makes X(j) an exact solution).
*
*  WORK    (workspace) COMPLEX*16 array, dimension (2*N)
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, and i is
*                <= N:  U(i,i) is exactly zero.  The factorization
*                       has been completed, but the factor U is exactly
*                       singular, so the solution and error bounds
*                       could not be computed.
*                = N+1: RCOND is less than machine precision.  The
*                       factorization has been completed, but the
*                       matrix A is singular to working precision, and
*                       the solution and error bounds have not been
*                       computed.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            COLEQU, EQUIL, NOFACT, NOTRAN, ROWEQU
      CHARACTER          NORM
      INTEGER            I, INFEQU, J, J1, J2
      DOUBLE PRECISION   AMAX, ANORM, BIGNUM, COLCND, RCMAX, RCMIN,
     $                   ROWCND, SMLNUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, ZLANGB
      EXTERNAL           LSAME, DLAMCH, ZLANGB
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZCOPY, ZGBCON, ZGBEQU, ZGBRFS, ZGBTRF,
     $                   ZGBTRS, ZLACPY, ZLAQGB
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      NOFACT = LSAME( FACT, 'N' )
      EQUIL = LSAME( FACT, 'E' )
      NOTRAN = LSAME( TRANS, 'N' )
      IF( NOFACT .OR. EQUIL ) THEN
         EQUED = 'N'
         ROWEQU = .FALSE.
         COLEQU = .FALSE.
      ELSE
         ROWEQU = LSAME( EQUED, 'R' ) .OR. LSAME( EQUED, 'B' )
         COLEQU = LSAME( EQUED, 'C' ) .OR. LSAME( EQUED, 'B' )
         SMLNUM = DLAMCH( 'Safe minimum' )
         BIGNUM = ONE / SMLNUM
      END IF
*
*     Test the input parameters.
*
      IF( .NOT.NOFACT .AND. .NOT.EQUIL .AND. .NOT.LSAME( FACT, 'F' ) )
     $     THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     $         LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( KL.LT.0 ) THEN
         INFO = -4
      ELSE IF( KU.LT.0 ) THEN
         INFO = -5
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -6
      ELSE IF( LDAB.LT.KL+KU+1 ) THEN
         INFO = -8
      ELSE IF( LDAFB.LT.2*KL+KU+1 ) THEN
         INFO = -10
      ELSE IF( LSAME( FACT, 'F' ) .AND. .NOT.
     $         ( ROWEQU .OR. COLEQU .OR. LSAME( EQUED, 'N' ) ) ) THEN
         INFO = -12
      ELSE
         IF( ROWEQU ) THEN
            RCMIN = BIGNUM
            RCMAX = ZERO
            DO 10 J = 1, N
               RCMIN = MIN( RCMIN, R( J ) )
               RCMAX = MAX( RCMAX, R( J ) )
   10       CONTINUE
            IF( RCMIN.LE.ZERO ) THEN
               INFO = -13
            ELSE IF( N.GT.0 ) THEN
               ROWCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
            ELSE
               ROWCND = ONE
            END IF
         END IF
         IF( COLEQU .AND. INFO.EQ.0 ) THEN
            RCMIN = BIGNUM
            RCMAX = ZERO
            DO 20 J = 1, N
               RCMIN = MIN( RCMIN, C( J ) )
               RCMAX = MAX( RCMAX, C( J ) )
   20       CONTINUE
            IF( RCMIN.LE.ZERO ) THEN
               INFO = -14
            ELSE IF( N.GT.0 ) THEN
               COLCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
            ELSE
               COLCND = ONE
            END IF
         END IF
         IF( INFO.EQ.0 ) THEN
            IF( LDB.LT.MAX( 1, N ) ) THEN
               INFO = -16
            ELSE IF( LDX.LT.MAX( 1, N ) ) THEN
               INFO = -18
            END IF
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGBSVX', -INFO )
         RETURN
      END IF
*
      IF( EQUIL ) THEN
*
*        Compute row and column scalings to equilibrate the matrix A.
*
         CALL ZGBEQU( N, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND,
     $                AMAX, INFEQU )
         IF( INFEQU.EQ.0 ) THEN
*
*           Equilibrate the matrix.
*
            CALL ZLAQGB( N, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND,
     $                   AMAX, EQUED )
            ROWEQU = LSAME( EQUED, 'R' ) .OR. LSAME( EQUED, 'B' )
            COLEQU = LSAME( EQUED, 'C' ) .OR. LSAME( EQUED, 'B' )
         END IF
      END IF
*
*     Scale the right-hand side.
*
      IF( NOTRAN ) THEN
         IF( ROWEQU ) THEN
            DO 40 J = 1, NRHS
               DO 30 I = 1, N
                  B( I, J ) = R( I )*B( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
      ELSE IF( COLEQU ) THEN
         DO 60 J = 1, NRHS
            DO 50 I = 1, N
               B( I, J ) = C( I )*B( I, J )
   50       CONTINUE
   60    CONTINUE
      END IF
*
      IF( NOFACT .OR. EQUIL ) THEN
*
*        Compute the LU factorization of the band matrix A.
*
         DO 70 J = 1, N
            J1 = MAX( J-KU, 1 )
            J2 = MIN( J+KL, N )
            CALL ZCOPY( J2-J1+1, AB( KU+1-J+J1, J ), 1,
     $                  AFB( KL+KU+1-J+J1, J ), 1 )
   70    CONTINUE
*
         CALL ZGBTRF( N, N, KL, KU, AFB, LDAFB, IPIV, INFO )
*
*        Return if INFO is non-zero.
*
         IF( INFO.NE.0 ) THEN
            IF( INFO.GT.0 )
     $         RCOND = ZERO
            RETURN
         END IF
      END IF
*
*     Compute the norm of the matrix A.
*
      IF( NOTRAN ) THEN
         NORM = '1'
      ELSE
         NORM = 'I'
      END IF
      ANORM = ZLANGB( NORM, N, KL, KU, AB, LDAB, RWORK )
*
*     Compute the reciprocal of the condition number of A.
*
      CALL ZGBCON( NORM, N, KL, KU, AFB, LDAFB, IPIV, ANORM, RCOND,
     $             WORK, RWORK, INFO )
*
*     Return if the matrix is singular to working precision.
*
      IF( RCOND.LT.DLAMCH( 'Epsilon' ) ) THEN
         INFO = N + 1
         RETURN
      END IF
*
*     Compute the solution matrix X.
*
      CALL ZLACPY( 'Full', N, NRHS, B, LDB, X, LDX )
      CALL ZGBTRS( TRANS, N, KL, KU, NRHS, AFB, LDAFB, IPIV, X, LDX,
     $             INFO )
*
*     Use iterative refinement to improve the computed solution and
*     compute error bounds and backward error estimates for it.
*
      CALL ZGBRFS( TRANS, N, KL, KU, NRHS, AB, LDAB, AFB, LDAFB, IPIV,
     $             B, LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO )
*
*     Transform the solution matrix X to a solution of the original
*     system.
*
      IF( NOTRAN ) THEN
         IF( COLEQU ) THEN
            DO 90 J = 1, NRHS
               DO 80 I = 1, N
                  X( I, J ) = C( I )*X( I, J )
   80          CONTINUE
   90       CONTINUE
            DO 100 J = 1, NRHS
               FERR( J ) = FERR( J ) / COLCND
  100       CONTINUE
         END IF
      ELSE IF( ROWEQU ) THEN
         DO 120 J = 1, NRHS
            DO 110 I = 1, N
               X( I, J ) = R( I )*X( I, J )
  110       CONTINUE
  120    CONTINUE
         DO 130 J = 1, NRHS
            FERR( J ) = FERR( J ) / ROWCND
  130    CONTINUE
      END IF
*
      RETURN
*
*     End of ZGBSVX
*
      END