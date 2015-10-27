      SUBROUTINE DCHKPB( DOTYPE, NN, NVAL, NNB, NBVAL, NRHS, THRESH,
     $                   TSTERR, NMAX, A, AFAC, AINV, B, X, XACT, WORK,
     $                   RWORK, IWORK, NOUT )
*
*  -- LAPACK test routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      LOGICAL            TSTERR
      INTEGER            NMAX, NN, NNB, NOUT, NRHS
      DOUBLE PRECISION   THRESH
*     ..
*     .. Array Arguments ..
      LOGICAL            DOTYPE( * )
      INTEGER            IWORK( * ), NBVAL( * ), NVAL( * )
      DOUBLE PRECISION   A( * ), AFAC( * ), AINV( * ), B( * ),
     $                   RWORK( * ), WORK( * ), X( * ), XACT( * )
*     ..
*
*  Purpose
*  =======
*
*  DCHKPB tests DPBTRF, -TRS, -RFS, and -CON.
*
*  Arguments
*  =========
*
*  DOTYPE  (input) LOGICAL array, dimension (NTYPES)
*          The matrix types to be used for testing.  Matrices of type j
*          (for 1 <= j <= NTYPES) are used for testing if DOTYPE(j) =
*          .TRUE.; if DOTYPE(j) = .FALSE., then type j is not used.
*
*  NN      (input) INTEGER
*          The number of values of N contained in the vector NVAL.
*
*  NVAL    (input) INTEGER array, dimension (NN)
*          The values of the matrix dimension N.
*
*  NNB     (input) INTEGER
*          The number of values of NB contained in the vector NBVAL.
*
*  NBVAL   (input) INTEGER array, dimension (NBVAL)
*          The values of the blocksize NB.
*
*  NRHS    (input) INTEGER
*          The number of right hand side vectors to be generated for
*          each linear system.
*
*  THRESH  (input) DOUBLE PRECISION
*          The threshold value for the test ratios.  A result is
*          included in the output file if RESULT >= THRESH.  To have
*          every test ratio printed, use THRESH = 0.
*
*  TSTERR  (input) LOGICAL
*          Flag that indicates whether error exits are to be tested.
*
*  NMAX    (input) INTEGER
*          The maximum value permitted for N, used in dimensioning the
*          work arrays.
*
*  A       (workspace) DOUBLE PRECISION array, dimension (NMAX*NMAX)
*
*  AFAC    (workspace) DOUBLE PRECISION array, dimension (NMAX*NMAX)
*
*  AINV    (workspace) DOUBLE PRECISION array, dimension (NMAX*NMAX)
*
*  B       (workspace) DOUBLE PRECISION array, dimension (NMAX*NRHS)
*
*  X       (workspace) DOUBLE PRECISION array, dimension (NMAX*NRHS)
*
*  XACT    (workspace) DOUBLE PRECISION array, dimension (NMAX*NRHS)
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension
*                      (NMAX*max(3,NRHS))
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension
*                      (max(NMAX,2*NRHS))
*
*  IWORK   (workspace) INTEGER array, dimension (NMAX)
*
*  NOUT    (input) INTEGER
*          The unit number for output.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      INTEGER            NTYPES, NTESTS
      PARAMETER          ( NTYPES = 8, NTESTS = 7 )
      INTEGER            NBW
      PARAMETER          ( NBW = 4 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ZEROT
      CHARACTER          DIST, PACKIT, TYPE, UPLO, XTYPE
      CHARACTER*3        PATH
      INTEGER            I, I1, I2, IKD, IMAT, IN, INB, INFO, IOFF,
     $                   IUPLO, IW, IZERO, K, KD, KL, KOFF, KU, LDA,
     $                   LDAB, MODE, N, NB, NERRS, NFAIL, NIMAT, NKD,
     $                   NRUN, NT
      DOUBLE PRECISION   AINVNM, ANORM, CNDNUM, RCOND, RCONDC
*     ..
*     .. Local Arrays ..
      INTEGER            ISEED( 4 ), ISEEDY( 4 ), KDVAL( NBW )
      DOUBLE PRECISION   RESULT( NTESTS )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DGET06, DLANGE, DLANSB
      EXTERNAL           DGET06, DLANGE, DLANSB
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAERH, ALAHD, ALASUM, DCOPY, DERRPO, DGET04,
     $                   DLACPY, DLARHS, DLATB4, DLATMS, DLAZRO, DPBCON,
     $                   DPBRFS, DPBT01, DPBT02, DPBT05, DPBTRF, DPBTRS,
     $                   DSWAP, XLAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Scalars in Common ..
      LOGICAL            LERR, OK
      CHARACTER*6        SRNAMT
      INTEGER            INFOT, NUNIT
*     ..
*     .. Common blocks ..
      COMMON             / INFOC / INFOT, NUNIT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Data statements ..
      DATA               ISEEDY / 1988, 1989, 1990, 1991 /
*     ..
*     .. Executable Statements ..
*
*     Initialize constants and the random number seed.
*
      PATH( 1: 1 ) = 'Double precision'
      PATH( 2: 3 ) = 'PB'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      DO 10 I = 1, 4
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE
*
*     Test the error exits
*
      IF( TSTERR )
     $   CALL DERRPO( PATH, NOUT )
      INFOT = 0
      CALL XLAENV( 2, 2 )
      KDVAL( 1 ) = 0
*
*     Do for each value of N in NVAL
*
      DO 90 IN = 1, NN
         N = NVAL( IN )
         LDA = MAX( N, 1 )
         XTYPE = 'N'
*
*        Set limits on the number of loop iterations.
*
         NKD = MAX( 1, MIN( N, 4 ) )
         NIMAT = NTYPES
         IF( N.EQ.0 )
     $      NIMAT = 1
*
         KDVAL( 2 ) = N + ( N+1 ) / 4
         KDVAL( 3 ) = ( 3*N-1 ) / 4
         KDVAL( 4 ) = ( N+1 ) / 4
*
         DO 80 IKD = 1, NKD
*
*           Do for KD = 0, (5*N+1)/4, (3N-1)/4, and (N+1)/4. This order
*           makes it easier to skip redundant values for small values
*           of N.
*
            KD = KDVAL( IKD )
            LDAB = KD + 1
*
*           Do first for UPLO = 'U', then for UPLO = 'L'
*
            DO 70 IUPLO = 1, 2
               KOFF = 1
               IF( IUPLO.EQ.1 ) THEN
                  UPLO = 'U'
                  KOFF = MAX( 1, KD+2-N )
                  PACKIT = 'Q'
               ELSE
                  UPLO = 'L'
                  PACKIT = 'B'
               END IF
*
               DO 60 IMAT = 1, NIMAT
*
*                 Do the tests only if DOTYPE( IMAT ) is true.
*
                  IF( .NOT.DOTYPE( IMAT ) )
     $               GO TO 60
*
*                 Skip types 2, 3, or 4 if the matrix size is too small.
*
                  ZEROT = IMAT.GE.2 .AND. IMAT.LE.4
                  IF( ZEROT .AND. N.LT.IMAT-1 )
     $               GO TO 60
*
                  IF( .NOT.ZEROT .OR. .NOT.DOTYPE( 1 ) ) THEN
*
*                    Set up parameters with DLATB4 and generate a test
*                    matrix with DLATMS.
*
                     CALL DLATB4( PATH, IMAT, N, N, TYPE, KL, KU, ANORM,
     $                            MODE, CNDNUM, DIST )
*
                     SRNAMT = 'DLATMS'
                     CALL DLATMS( N, N, DIST, ISEED, TYPE, RWORK, MODE,
     $                            CNDNUM, ANORM, KD, KD, PACKIT,
     $                            A( KOFF ), LDAB, WORK, INFO )
*
*                    Check error code from DLATMS.
*
                     IF( INFO.NE.0 ) THEN
                        CALL ALAERH( PATH, 'DLATMS', INFO, 0, UPLO, N,
     $                               N, KD, KD, -1, IMAT, NFAIL, NERRS,
     $                               NOUT )
                        GO TO 60
                     END IF
                  ELSE IF( IZERO.GT.0 ) THEN
*
*                    Use the same matrix for types 3 and 4 as for type
*                    2 by copying back the zeroed out column,
*
                     IW = 2*LDA + 1
                     IF( IUPLO.EQ.1 ) THEN
                        IOFF = ( IZERO-1 )*LDAB + KD + 1
                        CALL DCOPY( IZERO-I1, WORK( IW ), 1,
     $                              A( IOFF-IZERO+I1 ), 1 )
                        IW = IW + IZERO - I1
                        CALL DCOPY( I2-IZERO+1, WORK( IW ), 1,
     $                              A( IOFF ), MAX( LDAB-1, 1 ) )
                     ELSE
                        IOFF = ( I1-1 )*LDAB + 1
                        CALL DCOPY( IZERO-I1, WORK( IW ), 1,
     $                              A( IOFF+IZERO-I1 ),
     $                              MAX( LDAB-1, 1 ) )
                        IOFF = ( IZERO-1 )*LDAB + 1
                        IW = IW + IZERO - I1
                        CALL DCOPY( I2-IZERO+1, WORK( IW ), 1,
     $                              A( IOFF ), 1 )
                     END IF
                  END IF
*
*                 For types 2-4, zero one row and column of the matrix
*                 to test that INFO is returned correctly.
*
                  IZERO = 0
                  IF( ZEROT ) THEN
                     IF( IMAT.EQ.2 ) THEN
                        IZERO = 1
                     ELSE IF( IMAT.EQ.3 ) THEN
                        IZERO = N
                     ELSE
                        IZERO = N / 2 + 1
                     END IF
*
*                    Save the zeroed out row and column in WORK(*,3)
*
                     IW = 2*LDA
                     DO 20 I = 1, MIN( 2*KD+1, N )
                        WORK( IW+I ) = ZERO
   20                CONTINUE
                     IW = IW + 1
                     I1 = MAX( IZERO-KD, 1 )
                     I2 = MIN( IZERO+KD, N )
*
                     IF( IUPLO.EQ.1 ) THEN
                        IOFF = ( IZERO-1 )*LDAB + KD + 1
                        CALL DSWAP( IZERO-I1, A( IOFF-IZERO+I1 ), 1,
     $                              WORK( IW ), 1 )
                        IW = IW + IZERO - I1
                        CALL DSWAP( I2-IZERO+1, A( IOFF ),
     $                              MAX( LDAB-1, 1 ), WORK( IW ), 1 )
                     ELSE
                        IOFF = ( I1-1 )*LDAB + 1
                        CALL DSWAP( IZERO-I1, A( IOFF+IZERO-I1 ),
     $                              MAX( LDAB-1, 1 ), WORK( IW ), 1 )
                        IOFF = ( IZERO-1 )*LDAB + 1
                        IW = IW + IZERO - I1
                        CALL DSWAP( I2-IZERO+1, A( IOFF ), 1,
     $                              WORK( IW ), 1 )
                     END IF
                  END IF
*
*                 Do for each value of NB in NBVAL
*
                  DO 50 INB = 1, NNB
                     NB = NBVAL( INB )
                     CALL XLAENV( 1, NB )
*
*                    Compute the L*L' or U'*U factorization of the band
*                    matrix.
*
                     CALL DLACPY( 'Full', KD+1, N, A, LDAB, AFAC, LDAB )
                     SRNAMT = 'DPBTRF'
                     CALL DPBTRF( UPLO, N, KD, AFAC, LDAB, INFO )
*
*                    Check error code from DPBTRF.
*
                     IF( INFO.NE.IZERO ) THEN
                        CALL ALAERH( PATH, 'DPBTRF', INFO, IZERO, UPLO,
     $                               N, N, KD, KD, NB, IMAT, NFAIL,
     $                               NERRS, NOUT )
                        GO TO 50
                     END IF
*
*                    Skip the tests if INFO is not 0.
*
                     IF( INFO.NE.0 )
     $                  GO TO 50
*
*+    TEST 1
*                    Reconstruct matrix from factors and compute
*                    residual.
*
                     CALL DLACPY( 'Full', KD+1, N, AFAC, LDAB, AINV,
     $                            LDAB )
                     CALL DPBT01( UPLO, N, KD, A, LDAB, AINV, LDAB,
     $                            RWORK, RESULT( 1 ) )
                     NT = 1
*
*                    Only do other tests if this is the first blocksize.
*
                     IF( INB.GT.1 )
     $                  GO TO 30
*
*                    Form the inverse of A so we can get a good estimate
*                    of RCONDC = 1/(norm(A) * norm(inv(A))).
*
                     CALL DLAZRO( N, N, ZERO, ONE, AINV, LDA )
                     SRNAMT = 'DPBTRS'
                     CALL DPBTRS( UPLO, N, KD, N, AFAC, LDAB, AINV, LDA,
     $                            INFO )
*
*                    Compute RCONDC = 1/(norm(A) * norm(inv(A))).
*
                     ANORM = DLANSB( '1', UPLO, N, KD, A, LDAB, RWORK )
                     AINVNM = DLANGE( '1', N, N, AINV, LDA, RWORK )
                     IF( ANORM.LE.ZERO .OR. AINVNM.LE.ZERO ) THEN
                        RCONDC = ONE
                     ELSE
                        RCONDC = ( ONE / ANORM ) / AINVNM
                     END IF
*
*+    TEST 2
*                    Solve and compute residual for A * X = B.
*
                     SRNAMT = 'DLARHS'
                     CALL DLARHS( PATH, XTYPE, UPLO, ' ', N, N, KD, KD,
     $                            NRHS, A, LDAB, XACT, LDA, B, LDA,
     $                            ISEED, INFO )
                     XTYPE = 'C'
                     CALL DLACPY( 'Full', N, NRHS, B, LDA, X, LDA )
*
                     SRNAMT = 'DPBTRS'
                     CALL DPBTRS( UPLO, N, KD, NRHS, AFAC, LDAB, X, LDA,
     $                            INFO )
*
*                    Check error code from DPBTRS.
*
                     IF( INFO.NE.0 )
     $                  CALL ALAERH( PATH, 'DPBTRS', INFO, 0, UPLO, N,
     $                               N, KD, KD, NRHS, IMAT, NFAIL,
     $                               NERRS, NOUT )
*
                     CALL DLACPY( 'Full', N, NRHS, B, LDA, WORK, LDA )
                     CALL DPBT02( UPLO, N, KD, NRHS, A, LDAB, X, LDA,
     $                            WORK, LDA, RWORK, RESULT( 2 ) )
*
*+    TEST 3
*                    Check solution from generated exact solution.
*
                     CALL DGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC,
     $                            RESULT( 3 ) )
*
*+    TESTS 4, 5, and 6
*                    Use iterative refinement to improve the solution.
*
                     SRNAMT = 'DPBRFS'
                     CALL DPBRFS( UPLO, N, KD, NRHS, A, LDAB, AFAC,
     $                            LDAB, B, LDA, X, LDA, RWORK,
     $                            RWORK( NRHS+1 ), WORK, IWORK, INFO )
*
*                    Check error code from DPBRFS.
*
                     IF( INFO.NE.0 )
     $                  CALL ALAERH( PATH, 'DPBRFS', INFO, 0, UPLO, N,
     $                               N, KD, KD, NRHS, IMAT, NFAIL,
     $                               NERRS, NOUT )
*
                     CALL DGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC,
     $                            RESULT( 4 ) )
                     CALL DPBT05( UPLO, N, KD, NRHS, A, LDAB, B, LDA, X,
     $                            LDA, XACT, LDA, RWORK,
     $                            RWORK( NRHS+1 ), RESULT( 5 ) )
*
*+    TEST 7
*                    Get an estimate of RCOND = 1/CNDNUM.
*
                     SRNAMT = 'DPBCON'
                     CALL DPBCON( UPLO, N, KD, AFAC, LDAB, ANORM, RCOND,
     $                            WORK, IWORK, INFO )
*
*                    Check error code from DPBCON.
*
                     IF( INFO.NE.0 )
     $                  CALL ALAERH( PATH, 'DPBCON', INFO, 0, UPLO, N,
     $                               N, KD, KD, -1, IMAT, NFAIL, NERRS,
     $                               NOUT )
*
                     RESULT( 7 ) = DGET06( RCOND, RCONDC )
                     NT = 7
*
*                    Print information about the tests that did not pass
*                    the threshold.
*
   30                CONTINUE
                     DO 40 K = 1, NT
                        IF( RESULT( K ).GE.THRESH ) THEN
                           IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 )
     $                        CALL ALAHD( NOUT, PATH )
                           WRITE( NOUT, FMT = 9999 )UPLO, N, KD, NB,
     $                        IMAT, K, RESULT( K )
                           NFAIL = NFAIL + 1
                        END IF
   40                CONTINUE
                     NRUN = NRUN + NT
   50             CONTINUE
*
   60          CONTINUE
   70       CONTINUE
   80    CONTINUE
   90 CONTINUE
*
*     Print a summary of the results.
*
      CALL ALASUM( PATH, NOUT, NFAIL, NRUN, NERRS )
*
 9999 FORMAT( ' UPLO=''', A1, ''', N=', I5, ', KD=', I5, ', NB=', I4,
     $      ', type ', I2, ', test ', I2, ', ratio= ', G12.5 )
      RETURN
*
*     End of DCHKPB
*
      END
