      SUBROUTINE ZHSEIN( JOB, EIGSRC, INITV, SELECT, N, H, LDH, W, VL,
     $                   LDVL, VR, LDVR, MM, M, WORK, RWORK, INFO )
*
*  -- LAPACK routine (instrumented to count operations, version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          EIGSRC, INITV, JOB
      INTEGER            INFO, LDH, LDVL, LDVR, M, MM, N
*     ..
*     .. Array Arguments ..
      LOGICAL            SELECT( * )
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         H( LDH, * ), VL( LDVL, * ), VR( LDVR, * ),
     $                   W( * ), WORK( * )
*     ..
*     Common block to return operation count.
*     .. Common blocks ..
      COMMON             / LATIME / OPS, ITCNT
*     ..
*     .. Scalars in Common ..
      DOUBLE PRECISION   ITCNT, OPS
*     ..
*
*  Purpose
*  =======
*
*  ZHSEIN uses inverse iteration to find specified right and/or left
*  eigenvectors of a complex upper Hessenberg matrix H.
*
*  The right eigenvector x and the left eigenvector y of the matrix H
*  corresponding to an eigenvalue w are defined by:
*
*               H x = w x,     y' H = w y'
*
*  where y' denotes the conjugate transpose of the vector y.
*
*  Arguments
*  =========
*
*  JOB     (input) CHARACTER*1
*          = 'R': compute right eigenvectors only;
*          = 'L': compute left eigenvectors only;
*          = 'B': compute both right and left eigenvectors.
*
*  EIGSRC  (input) CHARACTER*1
*          Specifies the source of eigenvalues supplied in W:
*          = 'Q': the eigenvalues were found using ZHSEQR; thus, if
*                 H has zero sub-diagonal entries, and so is
*                 block-triangular, then the j-th eigenvalue can be
*                 assumed to be an eigenvalue of the block containing
*                 the j-th row/column.  This property allows ZHSEIN to
*                 perform inverse iteration on just one diagonal block.
*          = 'N': no assumptions are made on the correspondence
*                 between eigenvalues and diagonal blocks.  In this
*                 case, ZHSEIN must always perform inverse iteration
*                 using the whole matrix H.
*
*  INITV   (input) CHARACTER*1
*          Specifies whether initial starting vectors are supplied for
*          inverse iteration:
*          = 'N': no initial vectors are supplied;
*          = 'U': user-supplied initial vectors are stored in the arrays
*                 VL and/or VR.
*
*  SELECT  (input) LOGICAL array, dimension (N)
*          Specifies the eigenvectors to be computed. To select the
*          eigenvector corresponding to the eigenvalue W(j),
*          SELECT(j) must be set to .TRUE..
*
*  N       (input) INTEGER
*          The order of the matrix H.  N >= 0.
*
*  H       (input) COMPLEX*16 array, dimension (LDH,N)
*          The upper Hessenberg matrix H.
*
*  LDH     (input) INTEGER
*          The leading dimension of the array H.  LDH >= max(1,N).
*
*  W       (input/output) COMPLEX*16 array, dimension (N)
*          On entry, the eigenvalues of H.
*          On exit, the real parts of W may have been altered since
*          close eigenvalues are perturbed slightly in searching for
*          independent eigenvectors.
*
*  VL      (input/output) COMPLEX*16 array, dimension (LDVL,MM)
*          On entry, if INITV = 'U' and JOB = 'L' or 'B', VL must
*          contain starting vectors for the inverse iteration for the
*          left eigenvectors; the starting vector for each eigenvector
*          must be in the same column in which the eigenvector will be
*          stored.
*          On exit, if JOB = 'L' or 'B', the left eigenvectors
*          specified by SELECT will be stored consecutively in the
*          columns of VL, in the same order as their eigenvalues.
*          If JOB = 'R', VL is not referenced.
*
*  LDVL    (input) INTEGER
*          The leading dimension of the array VL.
*          LDVL >= max(1,N) if JOB = 'L' or 'B'; LDVL >= 1 otheerwise.
*
*  VR      (input/output) COMPLEX*16 array, dimension (LDVR,MM)
*          On entry, if INITV = 'U' and JOB = 'R' or 'B', VR must
*          contain starting vectors for the inverse iteration for the
*          right eigenvectors; the starting vector for each eigenvector
*          must be in the same column in which the eigenvector will be
*          stored.
*          On exit, if JOB = 'R' or 'B', the right eigenvectors
*          specified by SELECT will be stored consecutively in the
*          columns of VR, in the same order as their eigenvalues.
*          If JOB = 'L', VR is not referenced.
*
*  LDVR    (input) INTEGER
*          The leading dimension of the array VR.
*          LDVR >= max(1,N) if JOB = 'R' or 'B'; LDVR >= 1 otherwise.
*
*  MM      (input) INTEGER
*          The number of columns in the arrays VL and/or VR. MM >= M.
*
*  M       (output) INTEGER
*          The number of columns in the arrays VL and/or VR required to
*          store the eigenvectors (= the number of .TRUE. elements in
*          SELECT).
*
*  WORK    (workspace) COMPLEX*16 array, dimension (N*N)
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          > 0: INFO is the number of eigenvectors which failed to
*               converge; the corresponding columns of VL and/or VR are
*               set to zero.
*
*  Further Details
*  ===============
*
*  Each eigenvector is normalized so that the element of largest
*  magnitude has magnitude 1; here the magnitude of a complex number
*  (x,y) is taken to be |x|+|y|.
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
      DOUBLE PRECISION   RZERO
      PARAMETER          ( RZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            BOTHV, FROMQR, LEFTV, NOINIT, RIGHTV
      INTEGER            I, IINFO, K, KL, KR, KS, LDWORK
      DOUBLE PRECISION   EPS3, HNORM, OPST, SMLNUM, ULP, UNFL
      COMPLEX*16         CDUM, WK
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, ZLANHS
      EXTERNAL           LSAME, DLAMCH, ZLANHS
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLAEIN
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DIMAG, MAX
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
*     ..
*     .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( DBLE( CDUM ) ) + ABS( DIMAG( CDUM ) )
*     ..
*     .. Executable Statements ..
*
*     Decode and test the input parameters.
*
      BOTHV = LSAME( JOB, 'B' )
      RIGHTV = LSAME( JOB, 'R' ) .OR. BOTHV
      LEFTV = LSAME( JOB, 'L' ) .OR. BOTHV
*
      FROMQR = LSAME( EIGSRC, 'Q' )
*
      NOINIT = LSAME( INITV, 'N' )
*
*     Set M to the number of columns required to store the selected
*     eigenvectors.
*
      M = 0
      DO 10 K = 1, N
         IF( SELECT( K ) )
     $      M = M + 1
   10 CONTINUE
*
      INFO = 0
      IF( .NOT.RIGHTV .AND. .NOT.LEFTV ) THEN
         INFO = -1
      ELSE IF( .NOT.FROMQR .AND. .NOT.LSAME( EIGSRC, 'N' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.NOINIT .AND. .NOT.LSAME( INITV, 'U' ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDH.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDVL.LT.1 .OR. ( LEFTV .AND. LDVL.LT.N ) ) THEN
         INFO = -10
      ELSE IF( LDVR.LT.1 .OR. ( RIGHTV .AND. LDVR.LT.N ) ) THEN
         INFO = -12
      ELSE IF( MM.LT.M ) THEN
         INFO = -13
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZHSEIN', -INFO )
         RETURN
      END IF
***
*     Initialize
      OPST = 0
***
*
*     Quick return if possible.
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Set machine-dependent constants.
*
      UNFL = DLAMCH( 'Safe minimum' )
      ULP = DLAMCH( 'Precision' )
      SMLNUM = UNFL*( N / ULP )
*
*     Compute infinity-norm of matrix H.
*
      HNORM = ZLANHS( 'I', N, H, LDH, RWORK )
***
*     Increment op count for computing the norm of the matrix
      OPS = OPS + 5*N*( N+1 ) / 2 + 4*( N-1 )
***
      IF( HNORM.GT.RZERO ) THEN
         EPS3 = HNORM*ULP
      ELSE
         EPS3 = SMLNUM
      END IF
*
      LDWORK = N
*
      KL = 1
      IF( FROMQR ) THEN
         KR = 0
      ELSE
         KR = N
      END IF
      KS = 1
*
      DO 100 K = 1, N
         IF( SELECT( K ) ) THEN
*
*           Compute eigenvector(s) corresponding to W(K).
*
            IF( FROMQR ) THEN
*
*              If affiliation of eigenvalues is known, check whether
*              the matrix splits.
*
*              Determine KL and KR such that 1 <= KL <= K <= KR <= N
*              and H(KL,KL-1) and H(KR+1,KR) are zero (or KL = 1 or
*              KR = N).
*
*              Then inverse iteration can be performed with the
*              submatrix H(KL:N,KL:N) for a left eigenvector, and with
*              the submatrix H(1:KR,1:KR) for a right eigenvector.
*
               DO 20 I = K, KL + 1, -1
                  IF( H( I, I-1 ).EQ.ZERO )
     $               GO TO 30
   20          CONTINUE
   30          CONTINUE
               KL = I
               IF( K.GT.KR ) THEN
                  DO 40 I = K, N - 1
                     IF( H( I+1, I ).EQ.ZERO )
     $                  GO TO 50
   40             CONTINUE
   50             CONTINUE
                  KR = I
               END IF
            END IF
*
*           Perturb eigenvalue if it is close to any previous
*           selected eigenvalues affiliated to the submatrix
*           H(KL:KR,KL:KR). Close roots are modified by EPS3.
*
            WK = W( K )
   60       CONTINUE
            DO 70 I = K - 1, KL, -1
               IF( SELECT( I ) .AND. CABS1( W( I )-WK ).LT.EPS3 ) THEN
                  WK = WK + EPS3
                  GO TO 60
               END IF
   70       CONTINUE
            W( K ) = WK
***
*        Increment opcount for loop 70
            OPST = OPST + 2*( K-1 )
***
*
            IF( LEFTV ) THEN
*
*              Compute left eigenvector.
*
               CALL ZLAEIN( .FALSE., NOINIT, N-KL+1, H( KL, KL ), LDH,
     $                      WK, VL( KL, KS ), WORK, LDWORK, RWORK, EPS3,
     $                      SMLNUM, IINFO )
               IF( IINFO.GT.0 )
     $            INFO = INFO + 1
               DO 80 I = 1, KL - 1
                  VL( I, KS ) = ZERO
   80          CONTINUE
            END IF
            IF( RIGHTV ) THEN
*
*              Compute right eigenvector.
*
               CALL ZLAEIN( .TRUE., NOINIT, KR, H, LDH, WK, VR( 1, KS ),
     $                      WORK, LDWORK, RWORK, EPS3, SMLNUM, IINFO )
               IF( IINFO.GT.0 )
     $            INFO = INFO + 1
               DO 90 I = KR + 1, N
                  VR( I, KS ) = ZERO
   90          CONTINUE
            END IF
            KS = KS + 1
         END IF
  100 CONTINUE
*
***
*     Compute final op count
      OPS = OPS + OPST
***
      RETURN
*
*     End of ZHSEIN
*
      END