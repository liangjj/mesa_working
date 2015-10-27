      SUBROUTINE CTGSJA( JOBU, JOBV, JOBQ, M, P, N, K, L, A, LDA, B,
     $                   LDB, TOLA, TOLB, ALPHA, BETA, U, LDU, V, LDV,
     $                   Q, LDQ, WORK, NCYCLE, INFO )
*
*  -- LAPACK routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      CHARACTER          JOBQ, JOBU, JOBV
      INTEGER            INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N,
     $                   NCYCLE, P
      REAL               TOLA, TOLB
*     ..
*     .. Array Arguments ..
      REAL               ALPHA( * ), BETA( * )
      COMPLEX            A( LDA, * ), B( LDB, * ), Q( LDQ, * ),
     $                   U( LDU, * ), V( LDV, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  CTGSJA computes the generalized singular value decomposition (GSVD)
*  of two complex upper triangular (or trapezoidal) matrices A and B.
*
*  On entry, it is assumed that matrices A and B have the following
*  forms, which may be obtained by the preprocessing subroutine CGGSVP
*  for two general M-by-N matrix A and P-by-N matrix B:
*
*  If M-K-L >= 0
*
*     A = ( 0    A12  A13 ) K    ,  B = ( 0     0   B13 ) L
*         ( 0     0   A23 ) L           ( 0     0    0  ) P-L
*         ( 0     0    0  ) M-K-L        N-K-L  K    L
*          N-K-L  K    L
*
*  if M-K-L < 0
*
*      A = ( 0    A12  A13 ) K  ,   B = ( 0     0   B13 ) L
*          ( 0     0   A23 ) M-K        ( 0     0    0  ) P-L
*           N-K-L  K    L                N-K-L  K    L
*
*  where K-by-K matrix A12 and L-by-L matrix B13 are nonsingular upper
*  triangular. A23 is L-by-L upper triangular if M-K-L > 0, otherwise
*  A23 is L-by-(M-K) upper trapezoidal.
*
*  On exit,
*
*         U'*A*Q = D1*( 0 R ),    V'*B*Q = D2*( 0 R ),
*
*  where U, V and Q are unitary matrices, Z' denotes the conjugate
*  transpose of Z, R is a nonsingular upper triangular matrix, and D1
*  and D2 are ``diagonal'' matrices, which are of the following
*  structures:
*
*  If M-K-L >= 0,
*
*     U'*A*Q = D1*( 0 R )
*
*            = K     ( I  0 ) * (  0   R11  R12 ) K
*              L     ( 0  C )   (  0    0   R22 ) L
*              M-K-L ( 0  0 )    N-K-L  K    L
*                      K  L
*
*     V'*B*Q = D2*( 0 R )
*
*            = L     ( 0  S ) * (  0   R11  R12 ) K
*              P-L   ( 0  0 )   (  0    0   R22 ) L
*                      K  L      N-K-L  K    L
*  where
*
*    C = diag( ALPHA(K+1), ... , ALPHA(K+L) ),
*    S = diag( BETA(K+1),  ... , BETA(K+L) ),
*    C**2 + S**2 = I.
*
*    The nonsingular triangular matrix R = ( R11 R12 ) is stored
*                                          (  0  R22 )
*    in A(1:K+L,N-K-L+1:N) on exit.
*
*  If M-K-L < 0,
*
*  U'*A*Q = D1*( 0 R )
*
*         = K   ( I  0    0   ) * ( 0    R11  R12  R13  ) K
*           M-K ( 0  C    0   )   ( 0     0   R22  R23  ) M-K
*                 K M-K K+L-M     ( 0     0    0   R33  ) K+L-M
*                                  N-K-L  K   M-K  K+L-M
*
*  V'*B*Q = D2*( 0 R )
*
*         = M-K   ( 0  S    0   ) * ( 0    R11  R12  R13  ) K
*           K+L-M ( 0  0    I   )   ( 0     0   R22  R23  ) M-K
*           P-L   ( 0  0    0   )   ( 0     0    0   R33  ) K+L-M
*                   K M-K K+L-M      N-K-L  K   M-K  K+L-M
*
*  where
*  C = diag( ALPHA(K+1), ... , ALPHA(M) ),
*  S = diag( BETA(K+1),  ... , BETA(M) ),
*  C**2 + S**2 = I.
*
*  R = ( R11 R12 R13 ) is a nonsingular upper triangular matrix, the
*      (  0  R22 R23 )
*      (  0   0  R33 )
*  first M rows of R are stored in A(1:M, N-K-L+1:N) and R33 is stored
*  in B(M-K+1:L,N+M-K-L+1:N) on exit.
*
*  The computations of the unitary transformation matrices U, V and Q
*  are optional and may also be applied to the input unitary matrices U,
*  V and Q.
*
*  Arguments
*  =========
*
*  JOBU    (input) CHARACTER*1
*          = 'U':  U is overwritten on the input unitary matrix U;
*          = 'I':  U is initialized to the identity matrix;
*          = 'N':  U is not computed.
*
*  JOBV    (input) CHARACTER*1
*          = 'V':  V is overwritten on the input unitary matrix V;
*          = 'I':  V is initialized to the identity matrix;
*          = 'N':  V is not computed.
*
*  JOBQ    (input) CHARACTER*1
*          = 'Q':  Q is overwritten on the input unitary matrix Q;
*          = 'I':  Q is initialized to the identity matrix;
*          = 'N':  Q is not computed.
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  P       (input) INTEGER
*          The number of rows of the matrix B.  P >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrices A and B.  N >= 0.
*
*  K       (input) INTEGER
*  L       (input) INTEGER
*          K and L specify the subblocks in the input matrices A and B:
*          A23 = A(K+1:MIN(K+L,M),N-L+1:N) and B13 = B(1:L,,N-L+1:N)
*          of A and B, whose GSVD is going to be computed by CTGSJA.
*          See Further details.
*
*  A       (input/output) COMPLEX array, dimension (LDA,N)
*          On entry, the M-by-N matrix A.
*          On exit, A(N-K+1:N,1:MIN(K+L,M) ) contains the triangular
*          matrix R or part of R.  See Purpose for details.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,M).
*
*  B       (input/output) COMPLEX array, dimension (LDB,N)
*          On entry, the P-by-N matrix B.
*          On exit, if necessary, B(M-K+1:L,N+M-K-L+1:N) contains
*          a part of R.  See Purpose for details.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B. LDB >= max(1,P).
*
*  TOLA    (input) REAL
*  TOLB    (input) REAL
*          TOLA and TOLB are the convergence criteria for the Jacobi-
*          Kogbetliantz iteration procedure. Generally, they are the
*          same as used in the preprocessing step, say
*              TOLA = MAX(M,N)*norm(A)*MACHEPS,
*              TOLB = MAX(P,N)*norm(B)*MACHEPS.
*
*  ALPHA   (output) REAL array, dimension (N)
*  BETA    (output) REAL array, dimension (N)
*          On exit, ALPHA and BETA contain the generalized singular
*          value pairs of A and B;
*          If M-K-L >= 0,
*            ALPHA(1:K) = ONE,  ALPHA(K+1:K+L) = diag(C),
*            BETA(1:K)  = ZERO, BETA(K+1:K+L)  = diag(S),
*          and if M-K-L < 0,
*            ALPHA(1:K)= ONE,  ALPHA(K+1:M)= C, ALPHA(M+1:K+L)= ZERO
*            BETA(1:K) = ZERO, BETA(K+1:M) = S, BETA(M+1:K+L) = ONE.
*          Furthermore, if K+L < N,
*            ALPHA(K+L+1:N) = ZERO
*            BETA(K+L+1:N)  = ZERO.
*
*  U       (input/output) COMPLEX array, dimension (LDU,M)
*          On entry, if JOBU = 'U', U contains the unitary matrix U,
*          On exit, if JOBU = 'U', U  is overwritten on the input
*          unitary matrix U. If JOBU = 'I', U is first set to the
*          identity matrix.  If JOBU = 'N', U is not referenced.
*
*  LDU     (input) INTEGER
*          The leading dimension of the array U. LDU >= max(1,M).
*
*  V       (input/output) COMPLEX array, dimension (LDV,P)
*          On entry, if JOBV = 'V', V contains the unitary matrix V.
*          On exit, if JOBV = 'V', V  is overwritten on the input
*          unitary matrix V. If JOBV = 'I', U is first set to the
*          identity matrix.  If JOBV = 'N', V is not referenced.
*
*  LDV     (input) INTEGER
*          The leading dimension of the array V. LDV >= max(1,P).
*
*  Q       (input/output) COMPLEX array, dimension (LDQ,N)
*          On entry, if JOBQ = 'Q', Q contains the unitary matrix Q.
*          On exit, if JOBQ = 'Q', Q  is overwritten on the input
*          unitary matrix Q. If JOBQ = 'I', Q is first set to the
*          identity matrix.  If JOBQ = 'N', Q is not referenced.
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q. LDQ >= MAX(1,N).
*
*  WORK    (workspace) COMPLEX array, dimension (2*N)
*
*  NCYCLE  (output) INTEGER
*          The number of cycles required for convergence.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          = 1:  the procedure does not converge after MAXIT cycles.
*
*  Internal Parameters
*  ===================
*
*  MAXIT   INTEGER
*          MAXIT specifies the total loops that the iterative procedure
*          may take. If after MAXIT cycles, the routine fails to
*          converge, we return INFO = 1.
*
*  Further Details
*  ===============
*
*  CTGSJA essentially uses a variant of Kogbetliantz algorithm to reduce
*  min(L,M-K)-by-L triangular (or trapezoidal) matrix A23 and L-by-L
*  matrix B13 to the form:
*
*           U1'*A13*Q1 = C1*R1; V1'*B13*Q1 = S1*R1,
*
*  where U1, V1 and Q1 are unitary matrix, and Z' is the conjugate
*  transpose of Z.  C1 and S1 are diagonal matrices satisfying
*
*                C1**2 + S1**2 = I,
*
*  and R1 is an L-by-L nonsingular upper triangular matrix.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 40 )
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      COMPLEX            CZERO, CONE
      PARAMETER          ( CZERO = 0.0E+0, CONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
*
      LOGICAL            INITQ, INITU, INITV, UPPER, WANTQ, WANTU, WANTV
      INTEGER            I, J, KCYCLE
      REAL               A1, A3, B1, B3, CSQ, CSU, CSV, ERROR, GAMMA,
     $                   RWK, SSMIN
      COMPLEX            A2, B2, SNQ, SNU, SNV
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           CCOPY, CLAGS2, CLAPLL, CLAZRO, CROT, CSSCAL,
     $                   SLARTG, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, CONJG, MAX, MIN, REAL
*     ..
*     .. Executable Statements ..
*
*     Decode and test the input parameters
*
      INITU = LSAME( JOBU, 'I' )
      WANTU = INITU .OR. LSAME( JOBU, 'U' )
*
      INITV = LSAME( JOBV, 'I' )
      WANTV = INITV .OR. LSAME( JOBV, 'V' )
*
      INITQ = LSAME( JOBQ, 'I' )
      WANTQ = INITQ .OR. LSAME( JOBQ, 'Q' )
*
      INFO = 0
      IF( .NOT.( INITU .OR. WANTU .OR. LSAME( JOBU, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( INITV .OR. WANTV .OR. LSAME( JOBV, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF( .NOT.( INITQ .OR. WANTQ .OR. LSAME( JOBQ, 'N' ) ) ) THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( P.LT.0 ) THEN
         INFO = -5
      ELSE IF( N.LT.0 ) THEN
         INFO = -6
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -8
      ELSE IF( LDB.LT.MAX( 1, P ) ) THEN
         INFO = -10
      ELSE IF( LDU.LT.1 .OR. ( WANTU .AND. LDU.LT.M ) ) THEN
         INFO = -16
      ELSE IF( LDV.LT.1 .OR. ( WANTV .AND. LDV.LT.P ) ) THEN
         INFO = -18
      ELSE IF( LDQ.LT.1 .OR. ( WANTQ .AND. LDQ.LT.N ) ) THEN
         INFO = -20
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CTGSJA', -INFO )
         RETURN
      END IF
*
*     Initialize U, V and Q, if necessary
*
      IF( INITU )
     $   CALL CLAZRO( M, M, CZERO, CONE, U, LDU )
      IF( INITV )
     $   CALL CLAZRO( P, P, CZERO, CONE, V, LDV )
      IF( INITQ )
     $   CALL CLAZRO( N, N, CZERO, CONE, Q, LDQ )
*
*     Loop until convergence
*
      UPPER = .FALSE.
      DO 40 KCYCLE = 1, MAXIT
*
         UPPER = .NOT.UPPER
*
         DO 20 I = 1, L - 1
            DO 10 J = I + 1, L
*
               A1 = ZERO
               A2 = CZERO
               A3 = ZERO
               IF( K+I.LE.M )
     $            A1 = A( K+I, N-L+I )
               IF( K+J.LE.M )
     $            A3 = A( K+J, N-L+J )
*
               B1 = B( I, N-L+I )
               B3 = B( J, N-L+J )
*
               IF( UPPER ) THEN
                  IF( K+I.LE.M )
     $               A2 = A( K+I, N-L+J )
                  B2 = B( I, N-L+J )
               ELSE
                  IF( K+J.LE.M )
     $               A2 = A( K+J, N-L+I )
                  B2 = B( J, N-L+I )
               END IF
*
               CALL CLAGS2( UPPER, A1, A2, A3, B1, B2, B3, CSU, SNU,
     $                      CSV, SNV, CSQ, SNQ )
*
*              Update (K+I)-th and (K+J)-th rows of matrix A: U'*A
*
               IF( K+J.LE.M )
     $            CALL CROT( L, A( K+J, N-L+1 ), LDA, A( K+I, N-L+1 ),
     $                       LDA, CSU, CONJG( SNU ) )
*
*              Update I-th and J-th rows of matrix B: V'*B
*
               CALL CROT( L, B( J, N-L+1 ), LDB, B( I, N-L+1 ), LDB,
     $                    CSV, CONJG( SNV ) )
*
*              Update (N-L+I)-th and (N-L+J)-th columns of matrices
*              A and B: A*Q and B*Q
*
               CALL CROT( MIN( K+L, M ), A( 1, N-L+J ), 1,
     $                    A( 1, N-L+I ), 1, CSQ, SNQ )
*
               CALL CROT( L, B( 1, N-L+J ), 1, B( 1, N-L+I ), 1, CSQ,
     $                    SNQ )
*
               IF( UPPER ) THEN
                  IF( K+I.LE.M )
     $               A( K+I, N-L+J ) = CZERO
                  B( I, N-L+J ) = CZERO
               ELSE
                  IF( K+J.LE.M )
     $               A( K+J, N-L+I ) = CZERO
                  B( J, N-L+I ) = CZERO
               END IF
*
*              Ensure that the diagonal elements of A and B are real.
*
               IF( K+I.LE.M )
     $            A( K+I, N-L+I ) = REAL( A( K+I, N-L+I ) )
               IF( K+J.LE.M )
     $            A( K+J, N-L+J ) = REAL( A( K+J, N-L+J ) )
               B( I, N-L+I ) = REAL( B( I, N-L+I ) )
               B( J, N-L+J ) = REAL( B( J, N-L+J ) )
*
*              Update unitary matrices U, V, Q, if desired.
*
               IF( WANTU .AND. K+J.LE.M )
     $            CALL CROT( M, U( 1, K+J ), 1, U( 1, K+I ), 1, CSU,
     $                       SNU )
*
               IF( WANTV )
     $            CALL CROT( P, V( 1, J ), 1, V( 1, I ), 1, CSV, SNV )
*
               IF( WANTQ )
     $            CALL CROT( N, Q( 1, N-L+J ), 1, Q( 1, N-L+I ), 1, CSQ,
     $                       SNQ )
*
   10       CONTINUE
   20    CONTINUE
*
         IF( .NOT.UPPER ) THEN
*
*           The matrices A13 and B13 were lower triangular at the start
*           of the cycle, and are now upper triangular.
*
*           Convergence test: test the parallelism of the corresponding
*           rows of A and B.
*
            ERROR = ZERO
            DO 30 I = 1, MIN( L, M-K )
               CALL CCOPY( L-I+1, A( K+I, N-L+I ), LDA, WORK, 1 )
               CALL CCOPY( L-I+1, B( I, N-L+I ), LDB, WORK( L+1 ), 1 )
               CALL CLAPLL( L-I+1, WORK, 1, WORK( L+1 ), 1, SSMIN )
               ERROR = MAX( ERROR, SSMIN )
   30       CONTINUE
*
            IF( ABS( ERROR ).LE.REAL( N )*MIN( TOLA, TOLB ) )
     $         GO TO 50
         END IF
*
*        End of cycle loop
*
   40 CONTINUE
*
*     The algorithm has not converged after MAXIT cycles.
*
      INFO = 1
      GO TO 90
*
   50 CONTINUE
*
*     If ERROR <= N*MIN(TOLA,TOLB), then the algorithm has converged.
*     Compute the generalized singular value pairs (ALPHA, BETA), and
*     set the triangular matrix R to array A.
*
      DO 60 I = 1, K
         ALPHA( I ) = ONE
         BETA( I ) = ZERO
   60 CONTINUE
*
      DO 70 I = 1, MIN( L, M-K )
*
         A1 = A( K+I, N-L+I )
         B1 = B( I, N-L+I )
*
         IF( A1.NE.ZERO ) THEN
            GAMMA = B1 / A1
*
            IF( GAMMA.LT.ZERO )THEN
               CALL CSSCAL( L-I+1, -ONE, B( I, N-L+I ), LDB )
               IF( WANTV )
     $            CALL CSSCAL( P, -ONE, V( 1, I ), 1 )
            END IF
*
            CALL SLARTG( ABS( GAMMA ), ONE, BETA( K+I ), ALPHA( K+I ),
     $                   RWK )
*
            IF( ALPHA( K+I ).GE. BETA( K+I ) )THEN
               CALL CSSCAL( L-I+1, ONE / ALPHA( K+I ), A( K+I, N-L+I ),
     $                      LDA )
            ELSE
               CALL CSSCAL( L-I+1, ONE / BETA( K+I ), B( I, N-L+I ),
     $                      LDB )
               CALL CCOPY( L-I+1, B( I, N-L+I ), LDB, A( K+I, N-L+I ),
     $                      LDA )
            END IF
*
         ELSE
            ALPHA( K+I ) = ZERO
            BETA( K+I ) = ONE
            CALL CCOPY( L-I+1, B( I, N-L+I ), LDB, A( K+I, N-L+I ),
     $                  LDA )
         END IF
   70 CONTINUE
*
*     Post-assignment
*
      DO 80 I = M + 1, K + L
         ALPHA( I ) = ZERO
         BETA( I ) = ONE
   80 CONTINUE
*
      IF( K+L.LT.N )THEN
         DO 85 I = K + L + 1, N
            ALPHA( I ) = ZERO
            BETA( I ) = ZERO
   85    CONTINUE
      END IF
*
   90 CONTINUE
      NCYCLE = KCYCLE
*
      RETURN
*
*     End of CTGSJA
*
      END
