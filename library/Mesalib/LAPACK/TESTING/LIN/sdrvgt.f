      SUBROUTINE SDRVGT( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, A, AF,
     $                   B, X, XACT, WORK, RWORK, IWORK, NOUT )
*
*  -- LAPACK test routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      LOGICAL            TSTERR
      INTEGER            NN, NOUT, NRHS
      REAL               THRESH
*     ..
*     .. Array Arguments ..
      LOGICAL            DOTYPE( * )
      INTEGER            IWORK( * ), NVAL( * )
      REAL               A( * ), AF( * ), B( * ), RWORK( * ), WORK( * ),
     $                   X( * ), XACT( * )
*     ..
*
*  Purpose
*  =======
*
*  SDRVGT tests SGTSV and -SVX.
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
*  THRESH  (input) REAL
*          The threshold value for the test ratios.  A result is
*          included in the output file if RESULT >= THRESH.  To have
*          every test ratio printed, use THRESH = 0.
*
*  TSTERR  (input) LOGICAL
*          Flag that indicates whether error exits are to be tested.
*
*  A       (workspace) REAL array, dimension (NMAX*4)
*
*  AF      (workspace) REAL array, dimension (NMAX*4)
*
*  B       (workspace) REAL array, dimension (NMAX*NRHS)
*
*  X       (workspace) REAL array, dimension (NMAX*NRHS)
*
*  XACT    (workspace) REAL array, dimension (NMAX*NRHS)
*
*  WORK    (workspace) REAL array, dimension
*                      (NMAX*max(3,NRHS))
*
*  RWORK   (workspace) REAL array, dimension
*                      (max(NMAX,2*NRHS))
*
*  IWORK   (workspace) INTEGER array, dimension (2*NMAX)
*
*  NOUT    (input) INTEGER
*          The unit number for output.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
      INTEGER            NTYPES
      PARAMETER          ( NTYPES = 12 )
      INTEGER            NTESTS
      PARAMETER          ( NTESTS = 6 )
*     ..
*     .. Local Scalars ..
      LOGICAL            TRFCON, ZEROT
      CHARACTER          DIST, FACT, TRANS, TYPE
      CHARACTER*3        PATH
      INTEGER            I, IFACT, IMAT, IN, INFO, ITRAN, IX, IZERO, J,
     $                   K, K1, KL, KOFF, KU, LDA, MODE, N, NERRS,
     $                   NFAIL, NIMAT, NRUN, NT
      REAL               AINVNM, ANORM, ANORMI, ANORMO, COND, RCOND,
     $                   RCONDC, RCONDI, RCONDO
*     ..
*     .. Local Arrays ..
      CHARACTER          TRANSS( 3 )
      INTEGER            ISEED( 4 ), ISEEDY( 4 )
      REAL               RESULT( NTESTS ), Z( 3 )
*     ..
*     .. External Functions ..
      REAL               SASUM, SGET06, SLANGT
      EXTERNAL           SASUM, SGET06, SLANGT
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALADHD, ALAERH, ALASVM, SCOPY, SERRVX, SGET04,
     $                   SGTSV, SGTSVX, SGTT01, SGTT02, SGTT05, SGTTRF,
     $                   SGTTRS, SLACPY, SLAGTM, SLARNV, SLASET, SLATB4,
     $                   SLATMS, SSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
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
      DATA               ISEEDY / 0, 0, 0, 1 / , TRANSS / 'N', 'T',
     $                   'C' /
*     ..
*     .. Executable Statements ..
*
      PATH( 1: 1 ) = 'Single precision'
      PATH( 2: 3 ) = 'GT'
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
     $   CALL SERRVX( PATH, NOUT )
      INFOT = 0
*
      DO 140 IN = 1, NN
*
*        Do for each value of N in NVAL.
*
         N = NVAL( IN )
         LDA = MAX( 1, N )
         NIMAT = NTYPES
         IF( N.LE.0 )
     $      NIMAT = 1
*
         DO 130 IMAT = 1, NIMAT
*
*           Do the tests only if DOTYPE( IMAT ) is true.
*
            IF( .NOT.DOTYPE( IMAT ) )
     $         GO TO 130
*
*           Set up parameters with SLATB4.
*
            CALL SLATB4( PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE,
     $                   COND, DIST )
*
            ZEROT = IMAT.GE.8 .AND. IMAT.LE.10
            IF( IMAT.LE.6 ) THEN
*
*              Types 1-6:  generate matrices of known condition number.
*
               KOFF = MAX( 2-KU, 3-MAX( 1, N ) )
               SRNAMT = 'SLATMS'
               CALL SLATMS( N, N, DIST, ISEED, TYPE, RWORK, MODE, COND,
     $                      ANORM, KL, KU, 'Z', AF( KOFF ), 3, WORK,
     $                      INFO )
*
*              Check the error code from SLATMS.
*
               IF( INFO.NE.0 ) THEN
                  CALL ALAERH( PATH, 'SLATMS', INFO, 0, ' ', N, N, KL,
     $                         KU, -1, IMAT, NFAIL, NERRS, NOUT )
                  GO TO 130
               END IF
               IZERO = 0
*
               IF( N.GT.1 ) THEN
                  CALL SCOPY( N-1, AF( 4 ), 3, A, 1 )
                  CALL SCOPY( N-1, AF( 3 ), 3, A( 2*N ), 1 )
               END IF
               CALL SCOPY( N, AF( 2 ), 3, A( N ), 1 )
            ELSE
*
*              Types 7-12:  generate tridiagonal matrices with
*              unknown condition numbers.
*
               IF( .NOT.ZEROT .OR. .NOT.DOTYPE( 7 ) ) THEN
*
*                 Generate a matrix with elements from [-1,1].
*
                  CALL SLARNV( 2, ISEED, 3*N-2, A )
                  IF( ANORM.NE.ONE )
     $               CALL SSCAL( 3*N-2, ANORM, A, 1 )
               ELSE IF( IZERO.GT.0 ) THEN
*
*                 Reuse the last matrix by copying back the zeroed out
*                 elements.
*
                  IF( IZERO.EQ.1 ) THEN
                     A( N ) = Z( 2 )
                     IF( N.GT.1 )
     $                  A( 1 ) = Z( 3 )
                  ELSE IF( IZERO.EQ.N ) THEN
                     A( 3*N-2 ) = Z( 1 )
                     A( 2*N-1 ) = Z( 2 )
                  ELSE
                     A( 2*N-2+IZERO ) = Z( 1 )
                     A( N-1+IZERO ) = Z( 2 )
                     A( IZERO ) = Z( 3 )
                  END IF
               END IF
*
*              If IMAT > 7, set one column of the matrix to 0.
*
               IF( .NOT.ZEROT ) THEN
                  IZERO = 0
               ELSE IF( IMAT.EQ.8 ) THEN
                  IZERO = 1
                  Z( 2 ) = A( N )
                  A( N ) = ZERO
                  IF( N.GT.1 ) THEN
                     Z( 3 ) = A( 1 )
                     A( 1 ) = ZERO
                  END IF
               ELSE IF( IMAT.EQ.9 ) THEN
                  IZERO = N
                  Z( 1 ) = A( 3*N-2 )
                  Z( 2 ) = A( 2*N-1 )
                  A( 3*N-2 ) = ZERO
                  A( 2*N-1 ) = ZERO
               ELSE
                  IZERO = ( N+1 ) / 2
                  DO 20 I = IZERO, N - 1
                     A( 2*N-2+I ) = ZERO
                     A( N-1+I ) = ZERO
                     A( I ) = ZERO
   20             CONTINUE
                  A( 3*N-2 ) = ZERO
                  A( 2*N-1 ) = ZERO
               END IF
            END IF
*
            DO 120 IFACT = 1, 2
               IF( IFACT.EQ.1 ) THEN
                  FACT = 'F'
               ELSE
                  FACT = 'N'
               END IF
*
*              Compute the condition number for comparison with
*              the value returned by SGTSVX.
*
               IF( ZEROT ) THEN
                  IF( IFACT.EQ.1 )
     $               GO TO 120
                  RCONDO = ZERO
                  RCONDI = ZERO
*
               ELSE IF( IFACT.EQ.1 ) THEN
                  IF( N.GT.0 )
     $               CALL SCOPY( 3*N-2, A, 1, AF, 1 )
*
*                 Compute the 1-norm and infinity-norm of A.
*
                  ANORMO = SLANGT( '1', N, A, A( N ), A( 2*N ) )
                  ANORMI = SLANGT( 'I', N, A, A( N ), A( 2*N ) )
*
*                 Factor the matrix A.
*
                  CALL SGTTRF( N, AF, AF( N ), AF( 2*N ), AF( 3*N-1 ),
     $                         IWORK, INFO )
*
*                 Use SGTTRS to solve for one column at a time of
*                 inv(A), computing the maximum column sum as we go.
*
                  AINVNM = ZERO
                  DO 40 I = 1, N
                     DO 30 J = 1, N
                        X( J ) = ZERO
   30                CONTINUE
                     X( I ) = ONE
                     CALL SGTTRS( 'No transpose', N, 1, AF, AF( N ),
     $                            AF( 2*N ), AF( 3*N-1 ), IWORK, X, LDA,
     $                            INFO )
                     AINVNM = MAX( AINVNM, SASUM( N, X, 1 ) )
   40             CONTINUE
*
*                 Compute the 1-norm condition number of A.
*
                  IF( ANORMO.LE.ZERO .OR. AINVNM.LE.ZERO ) THEN
                     RCONDO = ONE
                  ELSE
                     RCONDO = ( ONE / ANORMO ) / AINVNM
                  END IF
*
*                 Use SGTTRS to solve for one column at a time of
*                 inv(A'), computing the maximum column sum as we go.
*
                  AINVNM = ZERO
                  DO 60 I = 1, N
                     DO 50 J = 1, N
                        X( J ) = ZERO
   50                CONTINUE
                     X( I ) = ONE
                     CALL SGTTRS( 'Transpose', N, 1, AF, AF( N ),
     $                            AF( 2*N ), AF( 3*N-1 ), IWORK, X, LDA,
     $                            INFO )
                     AINVNM = MAX( AINVNM, SASUM( N, X, 1 ) )
   60             CONTINUE
*
*                 Compute the infinity-norm condition number of A.
*
                  IF( ANORMI.LE.ZERO .OR. AINVNM.LE.ZERO ) THEN
                     RCONDI = ONE
                  ELSE
                     RCONDI = ( ONE / ANORMI ) / AINVNM
                  END IF
               END IF
*
               DO 110 ITRAN = 1, 3
                  TRANS = TRANSS( ITRAN )
                  IF( ITRAN.EQ.1 ) THEN
                     RCONDC = RCONDO
                  ELSE
                     RCONDC = RCONDI
                  END IF
*
*                 Generate NRHS random solution vectors.
*
                  IX = 1
                  DO 70 J = 1, NRHS
                     CALL SLARNV( 2, ISEED, N, XACT( IX ) )
                     IX = IX + LDA
   70             CONTINUE
*
*                 Set the right hand side.
*
                  CALL SLAGTM( TRANS, N, NRHS, ONE, A, A( N ), A( 2*N ),
     $                         XACT, LDA, ZERO, B, LDA )
*
                  IF( IFACT.EQ.2 .AND. ITRAN.EQ.1 ) THEN
*
*                    --- Test SGTSV  ---
*
*                    Solve the system using Gaussian elimination with
*                    partial pivoting.
*
                     IF( N.GT.0 )
     $                  CALL SCOPY( 3*N-2, A, 1, AF, 1 )
                     CALL SLACPY( 'Full', N, NRHS, B, LDA, X, LDA )
*
                     SRNAMT = 'SGTSV '
                     CALL SGTSV( N, NRHS, AF, AF( N ), AF( 2*N ), X,
     $                           LDA, INFO )
*
*                    Check error code from SGTSV .
*
                     IF( INFO.NE.IZERO )
     $                  CALL ALAERH( PATH, 'SGTSV ', INFO, IZERO, ' ',
     $                               N, N, 1, 1, NRHS, IMAT, NFAIL,
     $                               NERRS, NOUT )
                     NT = 1
                     IF( IZERO.EQ.0 ) THEN
*
*                       Check residual of computed solution.
*
                        CALL SLACPY( 'Full', N, NRHS, B, LDA, WORK,
     $                               LDA )
                        CALL SGTT02( TRANS, N, NRHS, A, A( N ),
     $                               A( 2*N ), X, LDA, WORK, LDA, RWORK,
     $                               RESULT( 2 ) )
*
*                       Check solution from generated exact solution.
*
                        CALL SGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC,
     $                               RESULT( 3 ) )
                        NT = 3
                     END IF
*
*                    Print information about the tests that did not pass
*                    the threshold.
*
                     DO 80 K = 2, NT
                        IF( RESULT( K ).GE.THRESH ) THEN
                           IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 )
     $                        CALL ALADHD( NOUT, PATH )
                           WRITE( NOUT, FMT = 9999 )'SGTSV ', N, IMAT,
     $                        K, RESULT( K )
                           NFAIL = NFAIL + 1
                        END IF
   80                CONTINUE
                     NRUN = NRUN + NT - 1
                  END IF
*
*                 --- Test SGTSVX ---
*
                  IF( IFACT.GT.1 ) THEN
*
*                    Initialize AF to zero.
*
                     DO 90 I = 1, 3*N - 2
                        AF( I ) = ZERO
   90                CONTINUE
                  END IF
                  CALL SLASET( 'Full', N, NRHS, ZERO, ZERO, X, LDA )
*
*                 Solve the system and compute the condition number and
*                 error bounds using SGTSVX.
*
                  SRNAMT = 'SGTSVX'
                  CALL SGTSVX( FACT, TRANS, N, NRHS, A, A( N ),
     $                         A( 2*N ), AF, AF( N ), AF( 2*N ),
     $                         AF( 3*N-1 ), IWORK, B, LDA, X, LDA,
     $                         RCOND, RWORK, RWORK( NRHS+1 ), WORK,
     $                         IWORK( N+1 ), INFO )
*
*                 Check the error code from SGTSVX.
*
                  IF( INFO.NE.IZERO )
     $               CALL ALAERH( PATH, 'SGTSVX', INFO, IZERO,
     $                            FACT // TRANS, N, N, 1, 1, NRHS, IMAT,
     $                            NFAIL, NERRS, NOUT )
*
                  IF( IFACT.GE.2 ) THEN
*
*                    Reconstruct matrix from factors and compute
*                    residual.
*
                     CALL SGTT01( N, A, A( N ), A( 2*N ), AF, AF( N ),
     $                            AF( 2*N ), AF( 3*N-1 ), IWORK, WORK,
     $                            LDA, RWORK, RESULT( 1 ) )
                     K1 = 1
                  ELSE
                     K1 = 2
                  END IF
*
                  IF( INFO.EQ.0 ) THEN
                     TRFCON = .FALSE.
*
*                    Check residual of computed solution.
*
                     CALL SLACPY( 'Full', N, NRHS, B, LDA, WORK, LDA )
                     CALL SGTT02( TRANS, N, NRHS, A, A( N ), A( 2*N ),
     $                            X, LDA, WORK, LDA, RWORK,
     $                            RESULT( 2 ) )
*
*                    Check solution from generated exact solution.
*
                     CALL SGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC,
     $                            RESULT( 3 ) )
*
*                    Check the error bounds from iterative refinement.
*
                     CALL SGTT05( TRANS, N, NRHS, A, A( N ), A( 2*N ),
     $                            B, LDA, X, LDA, XACT, LDA, RWORK,
     $                            RWORK( NRHS+1 ), RESULT( 4 ) )
                     NT = 5
                  END IF
*
*                 Print information about the tests that did not pass
*                 the threshold.
*
                  DO 100 K = K1, NT
                     IF( RESULT( K ).GE.THRESH ) THEN
                        IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 )
     $                     CALL ALADHD( NOUT, PATH )
                        WRITE( NOUT, FMT = 9998 )'SGTSVX', FACT, TRANS,
     $                     N, IMAT, K, RESULT( K )
                        NFAIL = NFAIL + 1
                     END IF
  100             CONTINUE
*
*                 Check the reciprocal of the condition number.
*
                  RESULT( 6 ) = SGET06( RCOND, RCONDC )
                  IF( RESULT( 6 ).GE.THRESH ) THEN
                     IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 )
     $                  CALL ALADHD( NOUT, PATH )
                     WRITE( NOUT, FMT = 9998 )'SGTSVX', FACT, TRANS, N,
     $                  IMAT, K, RESULT( K )
                     NFAIL = NFAIL + 1
                  END IF
                  NRUN = NRUN + NT - K1 + 2
*
  110          CONTINUE
  120       CONTINUE
  130    CONTINUE
  140 CONTINUE
*
*     Print a summary of the results.
*
      CALL ALASVM( PATH, NOUT, NFAIL, NRUN, NERRS )
*
 9999 FORMAT( 1X, A6, ', N =', I5, ', type ', I2, ', test ', I2,
     $      ', ratio = ', G12.5 )
 9998 FORMAT( 1X, A6, ', FACT=''', A1, ''', TRANS=''', A1, ''', N =',
     $      I5, ', type ', I2, ', test ', I2, ', ratio = ', G12.5 )
      RETURN
*
*     End of SDRVGT
*
      END