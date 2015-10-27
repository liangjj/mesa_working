      SUBROUTINE SGEGV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI,
     $                  BETA, VL, LDVL, VR, LDVR, WORK, LWORK,
     $                  INFO )
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
      REAL               A( LDA, * ), ALPHAI( * ), ALPHAR( * ),
     $                   B( LDB, * ), BETA( * ), VL( LDVL, * ),
     $                   VR( LDVR, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  For a pair of N-by-N real nonsymmetric matrices A, B:
*
*     compute the generalized eigenvalues (alphar +/- alphai*i, beta)
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
*                         H
*  l  such that  (A - w B) l = 0 .
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
*  A       (input/workspace) REAL array, dimension (LDA, N)
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
*  B       (input/workspace) REAL array, dimension (LDB, N)
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
*  ALPHAR  (output) REAL array, dimension (N)
*  ALPHAI  (output) REAL array, dimension (N)
*  BETA    (output) REAL array, dimension (N)
*
*          On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will
*          be the generalized eigenvalues.  If ALPHAI(j) is zero, then
*          the j-th eigenvalue is real; if positive, then the j-th and
*          (j+1)-st eigenvalues are a complex conjugate pair, with
*          ALPHAI(j+1) negative.
*
*          Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j)
*          may easily over- or underflow, and BETA(j) may even be zero.
*          Thus, the user should avoid naively computing the ratio
*          alpha/beta.  However, ALPHAR and ALPHAI will be always less
*          than and usually comparable with norm(A) in magnitude, and
*          BETA always less than and usually comparable with norm(B).
*
*  VL      (output) REAL array, dimension (LDVL,N)
*          If JOBVL = 'V', the left generalized eigenvectors.  (See
*          "Purpose", above.)  Real eigenvectors take one column,
*          complex take two columns, the first for the real part and
*          the second for the imaginary part.  Complex eigenvectors
*          correspond to an eigenvalue with positive imaginary part.
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
*  VR      (output) REAL array, dimension (LDVR,N)
*          If JOBVL = 'V', the right generalized eigenvectors.  (See
*          "Purpose", above.)  Real eigenvectors take one column,
*          complex take two columns, the first for the real part and
*          the second for the imaginary part.  Complex eigenvectors
*          correspond to an eigenvalue with positive imaginary part.
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
*  WORK    (workspace/output) REAL array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= max(1,8*N).
*          For good performance, LWORK must generally be larger.
*          To compute the optimal value of LWORK, call ILAENV to get
*          blocksizes (for SGEQRF, SORMQR, and SORGQR.)  Then compute:
*          NB  -- MAX of the blocksizes for SGEQRF, SORMQR, and SORGQR;
*          The optimal LWORK is:
*              2*N + MAX( 6*N, N*(NB+1) ).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          = 1,...,N:
*                The QZ iteration failed.  No eigenvectors have been
*                calculated, but ALPHAR(j), ALPHAI(j), and BETA(j)
*                should be correct for j=INFO+1,...,N.
*          > N:  errors that usually indicate LAPACK problems:
*                =N+1: error return from SGGBAL
*                =N+2: error return from SGEQRF
*                =N+3: error return from SORMQR
*                =N+4: error return from SORGQR
*                =N+5: error return from SGGHRD
*                =N+6: error return from SHGEQZ (other than failed
*                                                iteration)
*                =N+7: error return from STGEVC
*                =N+8: error return from SGGBAK (computing VL)
*                =N+9: error return from SGGBAK (computing VR)
*                =N+10: error return from SLASCL (various calls)
*
*  Further Details
*  ===============
*
*  Balancing
*  ---------
*
*  This driver calls SGGBAL to both permute and scale rows and columns
*  of A and B.  The permutations PL and PR are chosen so that PL*A*PR
*  and PL*B*R will be upper triangular except for the diagonal blocks
*  A(i:j,i:j) and B(i:j,i:j), with i and j as close together as
*  possible.  The diagonal scaling matrices DL and DR are chosen so
*  that the pair  DL*PL*A*PR*DR, DL*PL*B*PR*DR have entries close to
*  one (except for the entries that start out zero.)
*
*  After the eigenvalues and eigenvectors of the balanced matrices
*  have been computed, SGGBAK transforms the eigenvectors back to what
*  they would have been (in perfect arithmetic) if they had not been
*  balanced.
*
*  Contents of A and B on Exit
*  -------- -- - --- - -- ----
*
*  If any eigenvectors are computed (either JOBVL='V' or JOBVR='V' or
*  both), then on exit the arrays A and B will contain the real Schur
*  form[*] of the "balanced" versions of A and B.  If no eigenvectors
*  are computed, then only the diagonal blocks will be correct.
*
*  [*] See SHGEQZ, SGEGS, or read the book "Matrix Computations",
*      by Golub & van Loan, pub. by Johns Hopkins U. Press.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ILV, ILVL, ILVR, ILIMIT
      CHARACTER          CHTEMP
      INTEGER            ICOLS, IJOBVL, IJOBVR, IHI, IINFO, ILEFT, ILO,
     $                   IN, IRIGHT, IROWS, ITAU, IWORK, JC, JR, LWKMIN,
     $                   LWKOPT
      REAL               ABSAI, ABSAR, ABSB, ANRM, ANRM1, ANRM2, BNRM,
     $                   BNRM1, BNRM2, EPS, ONEPLS, SAFMAX, SAFMIN,
     $                   SALFAI, SALFAR, SBETA, SCALE, TEMP
*     ..
*     .. Local Arrays ..
      LOGICAL            LDUMMA( 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGEQRF, SGGBAK, SGGBAL, SGGHRD, SHGEQZ, SLACPY,
     $                   SLASCL, SLAZRO, SORGQR, SORMQR, STGEVC, XERBLA
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SLAMCH, SLANGE
      EXTERNAL           LSAME, SLAMCH, SLANGE
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, MAX
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
      LWKMIN = MAX( 8*N, 1 )
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
         INFO = -12
      ELSE IF( LDVR.LT.1 .OR. ( ILVR .AND. LDVR.LT.N ) ) THEN
         INFO = -14
      ELSE IF( LWORK.LT.LWKMIN ) THEN
         INFO = -16
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGEGV ', -INFO )
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
      EPS = SLAMCH( 'E' )*SLAMCH( 'B' )
      SAFMIN = SLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
      ONEPLS = ONE + ( 4*EPS )
*
*     Scale A
*
      ANRM = SLANGE( 'M', N, N, A, LDA, WORK )
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
         CALL SLASCL( 'G', -1, -1, ANRM, ONE, N, N, A, LDA, IINFO )
         IF( IINFO.NE.0 ) THEN
            INFO = N + 10
            RETURN
         END IF
      END IF
*
*     Scale B
*
      BNRM = SLANGE( 'M', N, N, B, LDB, WORK )
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
         CALL SLASCL( 'G', -1, -1, BNRM, ONE, N, N, B, LDB, IINFO )
         IF( IINFO.NE.0 ) THEN
            INFO = N + 10
            RETURN
         END IF
      END IF
*
*     Permute the matrix to make it more nearly triangular
*     Workspace layout:  (8*N words -- "work" requires 6*N words)
*        left_permutation, right_permutation, work...
*
      ILEFT = 1
      IRIGHT = N + 1
      IWORK = IRIGHT + N
      CALL SGGBAL( 'B', N, A, LDA, B, LDB, ILO, IHI, WORK( ILEFT ),
     $             WORK( IRIGHT ), WORK( IWORK ), IINFO )
      IF( IINFO.NE.0 ) THEN
         INFO = N + 1
         GO TO 120
      END IF
*
*     Reduce B to triangular form, and initialize VL and/or VR
*     Workspace layout:  ("work..." must have at least N words)
*        left_permutation, right_permutation, tau, work...
*
      IROWS = IHI + 1 - ILO
      IF( ILV ) THEN
         ICOLS = N + 1 - ILO
      ELSE
         ICOLS = IROWS
      END IF
      ITAU = IWORK
      IWORK = ITAU + IROWS
      CALL SGEQRF( IROWS, ICOLS, B( ILO, ILO ), LDB, WORK( ITAU ),
     $             WORK( IWORK ), LWORK+1-IWORK, IINFO )
      IF( IINFO.GE.0 )
     $   LWKOPT = MAX( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 )
      IF( IINFO.NE.0 ) THEN
         INFO = N + 2
         GO TO 120
      END IF
*
      CALL SORMQR( 'L', 'T', IROWS, ICOLS, IROWS, B( ILO, ILO ), LDB,
     $             WORK( ITAU ), A( ILO, ILO ), LDA, WORK( IWORK ),
     $             LWORK+1-IWORK, IINFO )
      IF( IINFO.GE.0 )
     $   LWKOPT = MAX( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 )
      IF( IINFO.NE.0 ) THEN
         INFO = N + 3
         GO TO 120
      END IF
*
      IF( ILVL ) THEN
         CALL SLAZRO( N, N, ZERO, ONE, VL, LDVL )
         CALL SLACPY( 'L', IROWS-1, IROWS-1, B( ILO+1, ILO ), LDB,
     $                VL( ILO+1, ILO ), LDVL )
         CALL SORGQR( IROWS, IROWS, IROWS, VL( ILO, ILO ), LDVL,
     $                WORK( ITAU ), WORK( IWORK ), LWORK+1-IWORK,
     $                IINFO )
         IF( IINFO.GE.0 )
     $      LWKOPT = MAX( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 )
         IF( IINFO.NE.0 ) THEN
            INFO = N + 4
            GO TO 120
         END IF
      END IF
*
      IF( ILVR )
     $   CALL SLAZRO( N, N, ZERO, ONE, VR, LDVR )
*
*     Reduce to generalized Hessenberg form
*
      IF( ILV ) THEN
*
*        Eigenvectors requested -- work on whole matrix.
*
         CALL SGGHRD( JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB, VL,
     $                LDVL, VR, LDVR, IINFO )
      ELSE
         CALL SGGHRD( 'N', 'N', IROWS, 1, IROWS, A( ILO, ILO ), LDA,
     $                B( ILO, ILO ), LDB, VL, LDVL, VR, LDVR,
     $                IINFO )
      END IF
      IF( IINFO.NE.0 ) THEN
         INFO = N + 5
         GO TO 120
      END IF
*
*     Perform QZ algorithm
*     Workspace layout:  ("work..." must have at least 1 word)
*        left_permutation, right_permutation, work...
*
      IWORK = ITAU
      IF( ILV ) THEN
         CHTEMP = 'S'
      ELSE
         CHTEMP = 'E'
      END IF
      CALL SHGEQZ( CHTEMP, JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB,
     $             ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR,
     $             WORK( IWORK ), LWORK+1-IWORK, IINFO )
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
         GO TO 120
      END IF
*
      IF( ILV ) THEN
*
*        Compute Eigenvectors  (STGEVC requires 6*N words of workspace)
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
         CALL STGEVC( 'B', CHTEMP, LDUMMA, N, A, LDA, B, LDB, VL,
     $                LDVL, VR, LDVR, N, IN, WORK( IWORK ), IINFO )
         IF( IINFO.NE.0 ) THEN
            INFO = N + 7
            GO TO 120
         END IF
*
*        Undo balancing on VL and VR, rescale
*
         IF( ILVL ) THEN
            CALL SGGBAK( 'B', 'L', N, ILO, IHI, WORK( ILEFT ),
     $                   WORK( IRIGHT ), N, VL, LDVL, IINFO )
            IF( IINFO.NE.0 ) THEN
               INFO = N + 8
               GO TO 120
            END IF
            DO 50 JC = 1, N
               IF( ALPHAI( JC ).LT.ZERO )
     $            GO TO 50
               TEMP = ZERO
               IF( ALPHAI( JC ).EQ.ZERO ) THEN
                  DO 10 JR = 1, N
                     TEMP = MAX( TEMP, ABS( VL( JR, JC ) ) )
   10             CONTINUE
               ELSE
                  DO 20 JR = 1, N
                     TEMP = MAX( TEMP, ABS( VL( JR, JC ) )+
     $                      ABS( VL( JR, JC+1 ) ) )
   20             CONTINUE
               END IF
               IF( TEMP.LT.SAFMIN )
     $            GO TO 50
               TEMP = ONE / TEMP
               IF( ALPHAI( JC ).EQ.ZERO ) THEN
                  DO 30 JR = 1, N
                     VL( JR, JC ) = VL( JR, JC )*TEMP
   30             CONTINUE
               ELSE
                  DO 40 JR = 1, N
                     VL( JR, JC ) = VL( JR, JC )*TEMP
                     VL( JR, JC+1 ) = VL( JR, JC+1 )*TEMP
   40             CONTINUE
               END IF
   50       CONTINUE
         END IF
         IF( ILVR ) THEN
            CALL SGGBAK( 'B', 'R', N, ILO, IHI, WORK( ILEFT ),
     $                   WORK( IRIGHT ), N, VR, LDVR, IINFO )
            IF( IINFO.NE.0 ) THEN
               INFO = N + 9
               GO TO 120
            END IF
            DO 100 JC = 1, N
               IF( ALPHAI( JC ).LT.ZERO )
     $            GO TO 100
               TEMP = ZERO
               IF( ALPHAI( JC ).EQ.ZERO ) THEN
                  DO 60 JR = 1, N
                     TEMP = MAX( TEMP, ABS( VR( JR, JC ) ) )
   60             CONTINUE
               ELSE
                  DO 70 JR = 1, N
                     TEMP = MAX( TEMP, ABS( VR( JR, JC ) )+
     $                      ABS( VR( JR, JC+1 ) ) )
   70             CONTINUE
               END IF
               IF( TEMP.LT.SAFMIN )
     $            GO TO 100
               TEMP = ONE / TEMP
               IF( ALPHAI( JC ).EQ.ZERO ) THEN
                  DO 80 JR = 1, N
                     VR( JR, JC ) = VR( JR, JC )*TEMP
   80             CONTINUE
               ELSE
                  DO 90 JR = 1, N
                     VR( JR, JC ) = VR( JR, JC )*TEMP
                     VR( JR, JC+1 ) = VR( JR, JC+1 )*TEMP
   90             CONTINUE
               END IF
  100       CONTINUE
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
      DO 110 JC = 1, N
         ABSAR = ABS( ALPHAR( JC ) )
         ABSAI = ABS( ALPHAI( JC ) )
         ABSB = ABS( BETA( JC ) )
         SALFAR = ANRM*ALPHAR( JC )
         SALFAI = ANRM*ALPHAI( JC )
         SBETA = BNRM*BETA( JC )
         ILIMIT = .FALSE.
         SCALE = ONE
*
*        Check for significant underflow in ALPHAI
*
         IF( ABS( SALFAI ).LT.SAFMIN .AND. ABSAI.GE.
     $       MAX( SAFMIN, EPS*ABSAR, EPS*ABSB ) ) THEN
            ILIMIT = .TRUE.
            SCALE = ( ONEPLS*SAFMIN / ANRM1 ) /
     $              MAX( ONEPLS*SAFMIN, ANRM2*ABSAI )
*
         ELSE IF( SALFAI.EQ.ZERO ) THEN
*
*           If insignificant underflow in ALPHAI, then make the
*           conjugate eigenvalue real.
*
            IF( ALPHAI( JC ).LT.ZERO .AND. JC.GT.1 ) THEN
               ALPHAI( JC-1 ) = ZERO
            ELSE IF( ALPHAI( JC ).GT.ZERO .AND. JC.LT.N ) THEN
               ALPHAI( JC+1 ) = ZERO
            END IF
         END IF
*
*        Check for significant underflow in ALPHAR
*
         IF( ABS( SALFAR ).LT.SAFMIN .AND. ABSAR.GE.
     $       MAX( SAFMIN, EPS*ABSAI, EPS*ABSB ) ) THEN
            ILIMIT = .TRUE.
            SCALE = MAX( SCALE, ( ONEPLS*SAFMIN / ANRM1 ) /
     $              MAX( ONEPLS*SAFMIN, ANRM2*ABSAR ) )
         END IF
*
*        Check for significant underflow in BETA
*
         IF( ABS( SBETA ).LT.SAFMIN .AND. ABSB.GE.
     $       MAX( SAFMIN, EPS*ABSAR, EPS*ABSAI ) ) THEN
            ILIMIT = .TRUE.
            SCALE = MAX( SCALE, ( ONEPLS*SAFMIN / BNRM1 ) /
     $              MAX( ONEPLS*SAFMIN, BNRM2*ABSB ) )
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
*        Recompute un-scaled ALPHAR, ALPHAI, BETA if necessary.
*
         IF( ILIMIT ) THEN
            SALFAR = ( SCALE*ALPHAR( JC ) )*ANRM
            SALFAI = ( SCALE*ALPHAI( JC ) )*ANRM
            SBETA = ( SCALE*BETA( JC ) )*BNRM
         END IF
         ALPHAR( JC ) = SALFAR
         ALPHAI( JC ) = SALFAI
         BETA( JC ) = SBETA
  110 CONTINUE
*
  120 CONTINUE
      WORK( 1 ) = LWKOPT
*
      RETURN
*
*     End of SGEGV
*
      END
