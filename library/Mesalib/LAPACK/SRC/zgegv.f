      SUBROUTINE ZGEGV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHA, BETA,
     $                  VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO )
*
*  -- LAPACK driver routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      CHARACTER          JOBVL, JOBVR
      INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         A( LDA, * ), ALPHA( * ), B( LDB, * ),
     $                   BETA( * ), VL( LDVL, * ), VR( LDVR, * ),
     $                   WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  For a pair of N-by-N complex nonsymmetric matrices A, B:
*
*     compute the generalized eigenvalues (alpha, beta)
*     compute the left and/or right generalized eigenvectors
*             (VL and VR)
*
*  The second action is optional -- see the description of JOBVL and
*  JOBVR below.
*
*  A generalized eigenvalue for a pair of matrices (A,B) is, roughly
*  speaking, a scalar w or a ratio  alpha/beta = w, such that  A - w*B
*  is singular.  It is usually represented as the pair (alpha,beta),
*  as there is a reasonable interpretation for beta=0, and even for
*  both being zero.  A good beginning reference is the book, "Matrix
*  Computations", by G. Golub & C. van Loan (Johns Hopkins U. Press)
*
*  A right generalized eigenvector corresponding to a generalized
*  eigenvalue  w  for a pair of matrices (A,B) is a vector  r  such
*  that  (A - w B) r = 0 .  A left generalized eigenvector is a vector
*  l  such that  (A - w B)**H l = 0 .
*
*  Note: this routine performs "full balancing" on A and B -- see
*  "Further Details", below.
*
*  Arguments
*  =========
*
*  JOBVL   (input) CHARACTER*1
*          = 'N':  do not compute the left generalized eigenvectors;
*          = 'V':  compute the left generalized eigenvectors.
*
*  JOBVR   (input) CHARACTER*1
*          = 'N':  do not compute the right generalized eigenvectors;
*          = 'V':  compute the right generalized eigenvectors.
*
*  N       (input) INTEGER
*          The number of rows and columns in the matrices A, B, VL, and
*          VR.  N >= 0.
*
*  A       (input/workspace) COMPLEX*16 array, dimension (LDA, N)
*          On entry, the first of the pair of matrices whose
*          generalized eigenvalues and (optionally) generalized
*          eigenvectors are to be computed.
*          On exit, the contents will have been destroyed.  (For a
*          description of the contents of A on exit, see "Further
*          Details", below.)
*
*  LDA     (input) INTEGER
*          The leading dimension of A.  LDA >= max(1,N).
*
*  B       (input/workspace) COMPLEX*16 array, dimension (LDB, N)
*          On entry, the second of the pair of matrices whose
*          generalized eigenvalues and (optionally) generalized
*          eigenvectors are to be computed.
*          On exit, the contents will have been destroyed.  (For a
*          description of the contents of B on exit, see "Further
*          Details", below.)
*
*  LDB     (input) INTEGER
*          The leading dimension of B.  LDB >= max(1,N).
*
*  ALPHA   (output) COMPLEX*16 array, dimension (N)
*  BETA    (output) COMPLEX*16 array, dimension (N)
*          On exit, ALPHA(j)/BETA(j), j=1,...,N, will be the
*          generalized eigenvalues.
*
*          Note: the quotients ALPHA(j)/BETA(j) may easily over- or
*          underflow, and BETA(j) may even be zero.  Thus, the user
*          should avoid naively computing the ratio alpha/beta.
*          However, ALPHA will be always less than and usually
*          comparable with norm(A) in magnitude, and BETA always less
*          than and usually comparable with norm(B).
*
*  VL      (output) COMPLEX*16 array, dimension (LDVL,N)
*          If JOBVL = 'V', the left generalized eigenvectors.  (See
*          "Purpose", above.)
*          Each eigenvector will be scaled so the largest component
*          will have abs(real part) + abs(imag. part) = 1, *except*
*          that for eigenvalues with alpha=beta=0, a zero vector will
*          be returned as the corresponding eigenvector.
*          Not referenced if JOBVL = 'N'.
*
*  LDVL    (input) INTEGER
*          The leading dimension of the matrix VL. LDVL >= 1, and
*          if JOBVL = 'V', LDVL >= N.
*
*  VR      (output) COMPLEX*16 array, dimension (LDVR,N)
*          If JOBVL = 'V', the right generalized eigenvectors.  (See
*          "Purpose", above.)
*          Each eigenvector will be scaled so the largest component
*          will have abs(real part) + abs(imag. part) = 1, *except*
*          that for eigenvalues with alpha=beta=0, a zero vector will
*          be returned as the corresponding eigenvector.
*          Not referenced if JOBVR = 'N'.
*
*  LDVR    (input) INTEGER
*          The leading dimension of the matrix VR. LDVR >= 1, and
*          if JOBVR = 'V', LDVR >= N.
*
*  WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= max(1,2*N).
*          For good performance, LWORK must generally be larger.
*          To compute the optimal value of LWORK, call ILAENV to get
*          blocksizes (for ZGEQRF, ZUNMQR, and CUNGQR.)  Then compute:
*          NB  -- MAX of the blocksizes for ZGEQRF, ZUNMQR, and CUNGQR;
*          The optimal LWORK is  MAX( 2*N, N*(NB+1) ).
*
*  RWORK   (workspace/output) DOUBLE PRECISION array, dimension (8*N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          =1,...,N:
*                The QZ iteration failed.  No eigenvectors have been
*                calculated, but ALPHA(j) and BETA(j) should be
*                correct for j=INFO+1,...,N.
*          > N:  errors that usually indicate LAPACK problems:
*                =N+1: error return from ZGGBAL
*                =N+2: error return from ZGEQRF
*                =N+3: error return from ZUNMQR
*                =N+4: error return from ZUNGQR
*                =N+5: error return from ZGGHRD
*                =N+6: error return from ZHGEQZ (other than failed
*                                               iteration)
*                =N+7: error return from ZTGEVC
*                =N+8: error return from ZGGBAK (computing VL)
*                =N+9: error return from ZGGBAK (computing VR)
*                =N+10: error return from ZLASCL (various calls)
*
*  Further Details
*  ===============
*
*  Balancing
*  ---------
*
*  This driver calls ZGGBAL to both permute and scale rows and columns
*  of A and B.  The permutations PL and PR are chosen so that PL*A*PR
*  and PL*B*R will be upper triangular except for the diagonal blocks
*  A(i:j,i:j) and B(i:j,i:j), with i and j as close together as
*  possible.  The diagonal scaling matrices DL and DR are chosen so
*  that the pair  DL*PL*A*PR*DR, DL*PL*B*PR*DR have entries close to
*  one (except for the entries that start out zero.)
*
*  After the eigenvalues and eigenvectors of the balanced matrices
*  have been computed, ZGGBAK transforms the eigenvectors back to what
*  they would have been (in perfect arithmetic) if they had not been
*  balanced.
*
*  Contents of A and B on Exit
*  -------- -- - --- - -- ----
*
*  If any eigenvectors are computed (either JOBVL='V' or JOBVR='V' or
*  both), then on exit the arrays A and B will contain the complex Schur
*  form[*] of the "balanced" versions of A and B.  If no eigenvectors
*  are computed, then only the diagonal blocks will be correct.
*
*  [*] In other words, upper triangular form.
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
      LOGICAL            ILIMIT, ILV, ILVL, ILVR
      CHARACTER          CHTEMP
      INTEGER            ICOLS, IHI, IINFO, IJOBVL, IJOBVR, ILEFT, ILO,
     $                   IN, IRIGHT, IROWS, IRWORK, ITAU, IWORK, JC, JR,
     $                   LWKMIN, LWKOPT
      DOUBLE PRECISION   ABSAI, ABSAR, ABSB, ANRM, ANRM1, ANRM2, BNRM,
     $                   BNRM1, BNRM2, EPS, SAFMAX, SAFMIN, SALFAI,
     $                   SALFAR, SBETA, SCALE, TEMP
      COMPLEX*16         X
*     ..
*     .. Local Arrays ..
      LOGICAL            LDUMMA( 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGEQRF, ZGGBAK, ZGGBAL, ZGGHRD, ZHGEQZ,
     $                   ZLACPY, ZLASCL, ZLAZRO, ZTGEVC, ZUNGQR, ZUNMQR
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, ZLANGE
      EXTERNAL           LSAME, DLAMCH, ZLANGE
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCMPLX, DIMAG, INT, MAX
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   ABS1
*     ..
*     .. Statement Function definitions ..
      ABS1( X ) = ABS( DBLE( X ) ) + ABS( DIMAG( X ) )
*     ..
*     .. Executable Statements ..
*
*     Decode the input arguments
*
      IF( LSAME( JOBVL, 'N' ) ) THEN
         IJOBVL = 1
         ILVL = .FALSE.
      ELSE IF( LSAME( JOBVL, 'V' ) ) THEN
         IJOBVL = 2
         ILVL = .TRUE.
      ELSE
         IJOBVL = -1
         ILVL = .FALSE.
      END IF
*
      IF( LSAME( JOBVR, 'N' ) ) THEN
         IJOBVR = 1
         ILVR = .FALSE.
      ELSE IF( LSAME( JOBVR, 'V' ) ) THEN
         IJOBVR = 2
         ILVR = .TRUE.
      ELSE
         IJOBVR = -1
         ILVR = .FALSE.
      END IF
      ILV = ILVL .OR. ILVR
*
*     Test the input arguments
*
      LWKMIN = MAX( 2*N, 1 )
      LWKOPT = LWKMIN
      INFO = 0
      IF( IJOBVL.LE.0 ) THEN
         INFO = -1
      ELSE IF( IJOBVR.LE.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDVL.LT.1 .OR. ( ILVL .AND. LDVL.LT.N ) ) THEN
         INFO = -11
      ELSE IF( LDVR.LT.1 .OR. ( ILVR .AND. LDVR.LT.N ) ) THEN
         INFO = -13
      ELSE IF( LWORK.LT.LWKMIN ) THEN
         INFO = -15
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGEGV ', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      WORK( 1 ) = LWKOPT
      IF( N.EQ.0 )
     $   RETURN
*
*     Get machine constants
*
      EPS = DLAMCH( 'E' )*DLAMCH( 'B' )
      SAFMIN = DLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
*
*     Scale A
*
      ANRM = ZLANGE( 'M', N, N, A, LDA, RWORK )
      ANRM1 = ANRM
      ANRM2 = ONE
      IF( ANRM.LT.ONE ) THEN
         IF( SAFMAX*ANRM.LT.ONE ) THEN
            ANRM1 = SAFMIN
            ANRM2 = SAFMAX*ANRM
         END IF
      END IF
*
      IF( ANRM.GT.ZERO ) THEN
         CALL ZLASCL( 'G', -1, -1, ANRM, ONE, N, N, A, LDA, IINFO )
         IF( IINFO.NE.0 ) THEN
            INFO = N + 10
            RETURN
         END IF
      END IF
*
*     Scale B
*
      BNRM = ZLANGE( 'M', N, N, B, LDB, RWORK )
      BNRM1 = BNRM
      BNRM2 = ONE
      IF( BNRM.LT.ONE ) THEN
         IF( SAFMAX*BNRM.LT.ONE ) THEN
            BNRM1 = SAFMIN
            BNRM2 = SAFMAX*BNRM
         END IF
      END IF
*
      IF( BNRM.GT.ZERO ) THEN
         CALL ZLASCL( 'G', -1, -1, BNRM, ONE, N, N, B, LDB, IINFO )
         IF( IINFO.NE.0 ) THEN
            INFO = N + 10
            RETURN
         END IF
      END IF
*
*     Permute the matrix to make it more nearly triangular
*     Also "balance" the matrix.
*
      ILEFT = 1
      IRIGHT = N + 1
      IRWORK = IRIGHT + N
      CALL ZGGBAL( 'B', N, A, LDA, B, LDB, ILO, IHI, RWORK( ILEFT ),
     $             RWORK( IRIGHT ), RWORK( IRWORK ), IINFO )
      IF( IINFO.NE.0 ) THEN
         INFO = N + 1
         GO TO 80
      END IF
*
*     Reduce B to triangular form, and initialize VL and/or VR
*
      IROWS = IHI + 1 - ILO
      IF( ILV ) THEN
         ICOLS = N + 1 - ILO
      ELSE
         ICOLS = IROWS
      END IF
      ITAU = 1
      IWORK = ITAU + IROWS
      CALL ZGEQRF( IROWS, ICOLS, B( ILO, ILO ), LDB, WORK( ITAU ),
     $             WORK( IWORK ), LWORK+1-IWORK, IINFO )
      IF( IINFO.GE.0 )
     $   LWKOPT = MAX( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 )
      IF( IINFO.NE.0 ) THEN
         INFO = N + 2
         GO TO 80
      END IF
*
      CALL ZUNMQR( 'L', 'C', IROWS, ICOLS, IROWS, B( ILO, ILO ), LDB,
     $             WORK( ITAU ), A( ILO, ILO ), LDA, WORK( IWORK ),
     $             LWORK+1-IWORK, IINFO )
      IF( IINFO.GE.0 )
     $   LWKOPT = MAX( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 )
      IF( IINFO.NE.0 ) THEN
         INFO = N + 3
         GO TO 80
      END IF
*
      IF( ILVL ) THEN
         CALL ZLAZRO( N, N, CZERO, CONE, VL, LDVL )
         CALL ZLACPY( 'L', IROWS-1, IROWS-1, B( ILO+1, ILO ), LDB,
     $                VL( ILO+1, ILO ), LDVL )
         CALL ZUNGQR( IROWS, IROWS, IROWS, VL( ILO, ILO ), LDVL,
     $                WORK( ITAU ), WORK( IWORK ), LWORK+1-IWORK,
     $                IINFO )
         IF( IINFO.GE.0 )
     $      LWKOPT = MAX( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 )
         IF( IINFO.NE.0 ) THEN
            INFO = N + 4
            GO TO 80
         END IF
      END IF
*
      IF( ILVR )
     $   CALL ZLAZRO( N, N, CZERO, CONE, VR, LDVR )
*
*     Reduce to generalized Hessenberg form
*
      IF( ILV ) THEN
*
*        Eigenvectors requested -- work on whole matrix.
*
         CALL ZGGHRD( JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB, VL,
     $                LDVL, VR, LDVR, IINFO )
      ELSE
         CALL ZGGHRD( 'N', 'N', IROWS, 1, IROWS, A( ILO, ILO ), LDA,
     $                B( ILO, ILO ), LDB, VL, LDVL, VR, LDVR, IINFO )
      END IF
      IF( IINFO.NE.0 ) THEN
         INFO = N + 5
         GO TO 80
      END IF
*
*     Perform QZ algorithm
*
      IWORK = ITAU
      IF( ILV ) THEN
         CHTEMP = 'S'
      ELSE
         CHTEMP = 'E'
      END IF
      CALL ZHGEQZ( CHTEMP, JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB,
     $             ALPHA, BETA, VL, LDVL, VR, LDVR, WORK( IWORK ),
     $             LWORK+1-IWORK, RWORK( IRWORK ), IINFO )
      IF( IINFO.GE.0 )
     $   LWKOPT = MAX( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 )
      IF( IINFO.NE.0 ) THEN
         IF( IINFO.GT.0 .AND. IINFO.LE.N ) THEN
            INFO = IINFO
         ELSE IF( IINFO.GT.N .AND. IINFO.LE.2*N ) THEN
            INFO = IINFO - N
         ELSE
            INFO = N + 6
         END IF
         GO TO 80
      END IF
*
      IF( ILV ) THEN
*
*        Compute Eigenvectors
*
         IF( ILVL ) THEN
            IF( ILVR ) THEN
               CHTEMP = 'B'
            ELSE
               CHTEMP = 'L'
            END IF
         ELSE
            CHTEMP = 'R'
         END IF
*
         CALL ZTGEVC( 'B', CHTEMP, LDUMMA, N, A, LDA, B, LDB, VL, LDVL,
     $                VR, LDVR, N, IN, WORK( IWORK ), RWORK( IRWORK ),
     $                IINFO )
         IF( IINFO.NE.0 ) THEN
            INFO = N + 7
            GO TO 80
         END IF
*
*        Undo balancing on VL and VR, rescale
*
         IF( ILVL ) THEN
            CALL ZGGBAK( 'B', 'L', N, ILO, IHI, RWORK( ILEFT ),
     $                   RWORK( IRIGHT ), N, VL, LDVL, IINFO )
            IF( IINFO.NE.0 ) THEN
               INFO = N + 8
               GO TO 80
            END IF
            DO 30 JC = 1, N
               TEMP = ZERO
               DO 10 JR = 1, N
                  TEMP = MAX( TEMP, ABS1( VL( JR, JC ) ) )
   10          CONTINUE
               IF( TEMP.LT.SAFMIN )
     $            GO TO 30
               TEMP = ONE / TEMP
               DO 20 JR = 1, N
                  VL( JR, JC ) = VL( JR, JC )*TEMP
   20          CONTINUE
   30       CONTINUE
         END IF
         IF( ILVR ) THEN
            CALL ZGGBAK( 'B', 'R', N, ILO, IHI, RWORK( ILEFT ),
     $                   RWORK( IRIGHT ), N, VR, LDVR, IINFO )
            IF( IINFO.NE.0 ) THEN
               INFO = N + 9
               GO TO 80
            END IF
            DO 60 JC = 1, N
               TEMP = ZERO
               DO 40 JR = 1, N
                  TEMP = MAX( TEMP, ABS1( VR( JR, JC ) ) )
   40          CONTINUE
               IF( TEMP.LT.SAFMIN )
     $            GO TO 60
               TEMP = ONE / TEMP
               DO 50 JR = 1, N
                  VR( JR, JC ) = VR( JR, JC )*TEMP
   50          CONTINUE
   60       CONTINUE
         END IF
*
*        End of eigenvector calculation
*
      END IF
*
*     Undo scaling in alpha, beta
*
*     Note: this does not give the alpha and beta for the unscaled
*     problem.
*
*     Un-scaling is limited to avoid underflow in alpha and beta
*     if they are significant.
*
      DO 70 JC = 1, N
         ABSAR = ABS( DBLE( ALPHA( JC ) ) )
         ABSAI = ABS( DIMAG( ALPHA( JC ) ) )
         ABSB = ABS( DBLE( BETA( JC ) ) )
         SALFAR = ANRM*DBLE( ALPHA( JC ) )
         SALFAI = ANRM*DIMAG( ALPHA( JC ) )
         SBETA = BNRM*DBLE( BETA( JC ) )
         ILIMIT = .FALSE.
         SCALE = ONE
*
*        Check for significant underflow in imaginary part of ALPHA
*
         IF( ABS( SALFAI ).LT.SAFMIN .AND. ABSAI.GE.
     $       MAX( SAFMIN, EPS*ABSAR, EPS*ABSB ) ) THEN
            ILIMIT = .TRUE.
            SCALE = ( SAFMIN / ANRM1 ) / MAX( SAFMIN, ANRM2*ABSAI )
         END IF
*
*        Check for significant underflow in real part of ALPHA
*
         IF( ABS( SALFAR ).LT.SAFMIN .AND. ABSAR.GE.
     $       MAX( SAFMIN, EPS*ABSAI, EPS*ABSB ) ) THEN
            ILIMIT = .TRUE.
            SCALE = MAX( SCALE, ( SAFMIN / ANRM1 ) /
     $              MAX( SAFMIN, ANRM2*ABSAR ) )
         END IF
*
*        Check for significant underflow in BETA
*
         IF( ABS( SBETA ).LT.SAFMIN .AND. ABSB.GE.
     $       MAX( SAFMIN, EPS*ABSAR, EPS*ABSAI ) ) THEN
            ILIMIT = .TRUE.
            SCALE = MAX( SCALE, ( SAFMIN / BNRM1 ) /
     $              MAX( SAFMIN, BNRM2*ABSB ) )
         END IF
*
*        Check for possible overflow when limiting scaling
*
         IF( ILIMIT ) THEN
            TEMP = ( SCALE*SAFMIN )*MAX( ABS( SALFAR ), ABS( SALFAI ),
     $             ABS( SBETA ) )
            IF( TEMP.GT.ONE )
     $         SCALE = SCALE / TEMP
            IF( SCALE.LT.ONE )
     $         ILIMIT = .FALSE.
         END IF
*
*        Recompute un-scaled ALPHA, BETA if necessary.
*
         IF( ILIMIT ) THEN
            SALFAR = ( SCALE*DBLE( ALPHA( JC ) ) )*ANRM
            SALFAI = ( SCALE*DIMAG( ALPHA( JC ) ) )*ANRM
            SBETA = ( SCALE*BETA( JC ) )*BNRM
         END IF
         ALPHA( JC ) = DCMPLX( SALFAR, SALFAI )
         BETA( JC ) = SBETA
   70 CONTINUE
*
   80 CONTINUE
      WORK( 1 ) = LWKOPT
*
      RETURN
*
*     End of ZGEGV
*
      END