      SUBROUTINE ZGGLSE( M, N, P, A, LDA, B, LDB, C, D, X, WORK, LWORK,
     $                   INFO )
*
*  -- LAPACK driver routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDB, LWORK, M, N, P
*     ..
*     .. Array Arguments ..
*
*  Purpose
*  =======
*
*  ZGGLSE solves the linear equality constrained least squares (LSE)
*  problem:
*
*          minimize || A*x - c ||_2   subject to B*x = d
*
*  using a generalized RQ factorization of matrices A and B, where A is
*  M-by-N, B is P-by-N, assume P <= N <= M+P, and ||.||_2 denotes vector
*  2-norm.  It is assumed that
*
*                       rank(B) = P                                  (1)
*
*  and the null spaces of A and B intersect only trivially, i.e.,
*
*   intersection of Null(A) and Null(B) = {0} <=> rank( ( A ) ) = N  (2)
*                                                     ( ( B ) )
*
*  where N(A) denotes the null space of matrix A. Conditions (1) and (2)
*  ensure that the problem LSE has a unique solution.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrices A and B. N >= 0.
*          assume that P <= N <= M+P.
*
*  P       (input) INTEGER
*          The number of rows of the matrix B.  P >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the P-by-M matrix A.
*          On exit, A is destroyed.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,M).
*
*  B       (input/output) COMPLEX*16 array, dimension (LDB,N)
*          On entry, the P-by-N matrix B.
*          On exit, B is destroyed.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B. LDB >= max(1,P).
*
*  C       (input/output) COMPLEX*16 array, dimension (M)
*          On entry, C contains the right hand side vector for the
*          least squares part of the LSE problem.
*          On exit, the residual sum of squares for the solution
*          is given by the sum of squares of elements N-P+1 to M of
*          vector C.
*
*  D       (input/output) COMPLEX*16 array, dimension (P)
*          On entry, D contains the right hand side vector for the
*          constrained equation.
*          On exit, D is destroyed.
*
*  X       (output) COMPLEX*16 array, dimension (N)
*          On exit, X is the solution of the problem LSE.
*
*  WORK    (workspace) COMPLEX*16 array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK. LWORK >= N+P+max(N,M,P).
*          For optimum performance LWORK >=
*          N+P+max(M,P,N)*max(NB1,NB2), where NB1 is the optimal
*          blocksize for the QR factorization of M-by-N matrix A.
*          NB2 is the optimal blocksize for the RQ factorization of
*          P-by-N matrix B.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  =====================================================================
*
      COMPLEX*16         A( LDA, * ), B( LDB, * ), C( * ), D( * ),
     $                   WORK( * ), X( * )
*     ..
*     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            MN
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZAXPY, ZCOPY, ZGEMV, ZGGRQF, ZTRMV,
     $                   ZTRSV, ZUNMQR, ZUNMRQ
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test input parameters
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( P.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.P .OR. N.GT.M+P ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, P ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGGLSE', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. MAX( M, P ).EQ.0 )
     $   RETURN
*
*     Compute the GRQ factorization of matrices B and A:
*
*                  B*Q' = ( 0  T ) P
*                          N-P P
*     and if M >= N,
*
*               Z'*A*Q' = ( R11 R12 ) N-P,
*                         (  0  R22 ) P
*                         (  0   0  ) M-N
*                           N-P  P
*     or if  if M < N,
*
*               Z'*A*Q' = ( R11   R12   R13 ) N-P  ,
*                         (  0    R22   R23 ) M+P-N
*                           N-P  M+P-N  N-M
*
*     where T, R11 and R22 are upper triangular matrices, Q and Z are
*     unitary matrix.
*
      MN = MIN( M, N )
      CALL ZGGRQF( P, M, N, B, LDB, WORK, A, LDA, WORK( P+1 ),
     $             WORK( P+MN+1 ), LWORK-P-MN, INFO )
*
*     Update c = Z'*c
*
      CALL ZUNMQR( 'Left', 'Conjugate Transpose', M, 1, MIN( M, N ), A,
     $             LDA, WORK( P+1 ), C, MAX( 1, M ), WORK( P+MN+1 ),
     $             LWORK-P-MN, INFO )
*
*     Solve T*x2 = d for x2
*
      CALL ZTRSV( 'Upper', 'No transpose', 'Non unit', P, B( 1, N-P+1 ),
     $            LDB, D, 1 )
*
*     Update c(1:N-P)
*
      CALL ZGEMV( 'No transpose', N-P, P, -ONE, A( 1, N-P+1 ), LDA, D,
     $            1, ONE, C, 1 )
*
*     Sovle R11*x1 = c1 for x1
*
      CALL ZTRSV( 'Upper', 'No transpose', 'Non unit', N-P, A, LDA, C,
     $            1 )
*
*     Put the solutions in X
*
      CALL ZCOPY( N-P, C, 1, X, 1 )
      CALL ZCOPY( P, D, 1, X( N-P+1 ), 1 )
*
*     Compute the residual vector:
*
      IF( M.GE.N ) THEN
         CALL ZTRMV( 'Upper', 'No transpose', 'Non unit', P,
     $               A( N-P+1, N-P+1 ), LDA, D, 1 )
         CALL ZAXPY( P, -ONE, D, 1, C( N-P+1 ), 1 )
      ELSE
         CALL ZTRMV( 'Upper', 'No transpose', 'Non unit', M+P-N,
     $               A( N-P+1, N-P+1 ), LDA, D, 1 )
         CALL ZAXPY( M+P-N, -ONE, D, 1, C( N-P+1 ), 1 )
         CALL ZGEMV( 'No transpose', M+P-N, N-M, -ONE, A( N-P+1, M+1 ),
     $               LDA, D( M+P-N+1 ), 1, ONE, C( N-P+1 ), 1 )
      END IF
*
*     Backward transformation x := Q'*x
*
      CALL ZUNMRQ( 'Left', 'Conjugate Transpose', N, 1, P, B, LDB,
     $             WORK( 1 ), X, N, WORK( P+MN+1 ), LWORK-P-MN, INFO )
*
      RETURN
*
*     End of ZGGLSE
*
      END
