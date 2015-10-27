      SUBROUTINE CDRVPP( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, NMAX,
     $                   A, AFAC, ASAV, B, BSAV, X, XACT, S, WORK,
     $                   RWORK, NOUT )
*
*  -- LAPACK test routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      LOGICAL            TSTERR
      INTEGER            NMAX, NN, NOUT, NRHS
      REAL               THRESH
*     ..
*     .. Array Arguments ..
      LOGICAL            DOTYPE( * )
      INTEGER            NVAL( * )
      REAL               RWORK( * ), S( * )
      COMPLEX            A( * ), AFAC( * ), ASAV( * ), B( * ),
     $                   BSAV( * ), WORK( * ), X( * ), XACT( * )
*     ..
*
*  Purpose
*  =======
*
*  CDRVPP tests the driver routines CPPSV and -SVX.
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
*  NRHS    (input) INTEGER
*          The number of right hand side vectors to be generated for
*          each linear system.
*
*  THRESH  (input) REAL
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
*  A       (workspace) COMPLEX array, dimension (NMAX*(NMAX+1)/2)
*
*  AFAC    (workspace) COMPLEX array, dimension (NMAX*(NMAX+1)/2)
*
*  ASAV    (workspace) COMPLEX array, dimension (NMAX*(NMAX+1)/2)
*
*  B       (workspace) COMPLEX array, dimension (NMAX*NRHS)
*
*  BSAV    (workspace) COMPLEX array, dimension (NMAX*NRHS)
*
*  X       (workspace) COMPLEX array, dimension (NMAX*NRHS)
*
*  XACT    (workspace) COMPLEX array, dimension (NMAX*NRHS)
*
*  S       (workspace) REAL array, dimension (NMAX)
*
*  WORK    (workspace) COMPLEX array, dimension
*                      (NMAX*max(3,NRHS))
*
*  RWORK   (workspace) REAL array, dimension (NMAX+2*NRHS)
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
      PARAMETER          ( NTYPES = 9 )
      INTEGER            NTESTS
      PARAMETER          ( NTESTS = 6 )
*     ..
*     .. Local Scalars ..
      LOGICAL            EQUIL, NOFACT, PREFAC, ZEROT
      CHARACTER          DIST, EQUED, FACT, PACKIT, TYPE, UPLO, XTYPE
      CHARACTER*3        PATH
      INTEGER            I, IEQUED, IFACT, IMAT, IN, INFO, IOFF, IUPLO,
     $                   IZERO, K, K1, KL, KU, LDA, MODE, N, NERRS,
     $                   NFACT, NFAIL, NIMAT, NPP, NRUN, NT
      REAL               AINVNM, AMAX, ANORM, CNDNUM, RCOND, RCONDC,
     $                   ROLDC, SCOND
*     ..
*     .. Local Arrays ..
      CHARACTER          EQUEDS( 2 ), FACTS( 3 ), PACKS( 2 ), UPLOS( 2 )
      INTEGER            ISEED( 4 ), ISEEDY( 4 )
      REAL               RESULT( NTESTS )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               CLANSP, SGET06
      EXTERNAL           LSAME, CLANSP, SGET06
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALADHD, ALAERH, ALASVM, CCOPY, CERRVX, CGET04,
     $                   CLACPY, CLAQSP, CLARHS, CLASET, CLATB4, CLATMS,
     $                   CPPEQU, CPPSV, CPPSVX, CPPT01, CPPT02, CPPT05,
     $                   CPPTRF, CPPTRI
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
*     .. Intrinsic Functions ..
      INTRINSIC          CMPLX, MAX
*     ..
*     .. Data statements ..
      DATA               ISEEDY / 1988, 1989, 1990, 1991 /
      DATA               UPLOS / 'U', 'L' / , FACTS / 'F', 'N', 'E' / ,
     $                   PACKS / 'C', 'R' / , EQUEDS / 'N', 'Y' /
*     ..
*     .. Executable Statements ..
*
*     Initialize constants and the random number seed.
*
      PATH( 1: 1 ) = 'Complex precision'
      PATH( 2: 3 ) = 'PP'
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
     $   CALL CERRVX( PATH, NOUT )
      INFOT = 0
*
*     Do for each value of N in NVAL
*
      DO 150 IN = 1, NN
         N = NVAL( IN )
         LDA = MAX( N, 1 )
         NPP = N*( N+1 ) / 2
         XTYPE = 'N'
         NIMAT = NTYPES
         IF( N.LE.0 )
     $      NIMAT = 1
*
         DO 140 IMAT = 1, NIMAT
*
*           Do the tests only if DOTYPE( IMAT ) is true.
*
            IF( .NOT.DOTYPE( IMAT ) )
     $         GO TO 140
*
*           Skip types 3, 4, or 5 if the matrix size is too small.
*
            ZEROT = IMAT.GE.3 .AND. IMAT.LE.5
            IF( ZEROT .AND. N.LT.IMAT-2 )
     $         GO TO 140
*
*           Do first for UPLO = 'U', then for UPLO = 'L'
*
            DO 130 IUPLO = 1, 2
               UPLO = UPLOS( IUPLO )
               PACKIT = PACKS( IUPLO )
*
*              Set up parameters with CLATB4 and generate a test matrix
*              with CLATMS.
*
               CALL CLATB4( PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE,
     $                      CNDNUM, DIST )
               RCONDC = ONE / CNDNUM
*
               SRNAMT = 'CLATMS'
               CALL CLATMS( N, N, DIST, ISEED, TYPE, RWORK, MODE,
     $                      CNDNUM, ANORM, KL, KU, PACKIT, A, LDA, WORK,
     $                      INFO )
*
*              Check error code from CLATMS.
*
               IF( INFO.NE.0 ) THEN
                  CALL ALAERH( PATH, 'CLATMS', INFO, 0, UPLO, N, N, -1,
     $                         -1, -1, IMAT, NFAIL, NERRS, NOUT )
                  GO TO 130
               END IF
*
*              For types 3-5, zero one row and column of the matrix to
*              test that INFO is returned correctly.
*
               IF( ZEROT ) THEN
                  IF( IMAT.EQ.3 ) THEN
                     IZERO = 1
                  ELSE IF( IMAT.EQ.4 ) THEN
                     IZERO = N
                  ELSE
                     IZERO = N / 2 + 1
                  END IF
*
*                 Set row and column IZERO of A to 0.
*
                  IF( IUPLO.EQ.1 ) THEN
                     IOFF = ( IZERO-1 )*IZERO / 2
                     DO 20 I = 1, IZERO - 1
                        A( IOFF+I ) = ZERO
   20                CONTINUE
                     IOFF = IOFF + IZERO
                     DO 30 I = IZERO, N
                        A( IOFF ) = ZERO
                        IOFF = IOFF + I
   30                CONTINUE
                  ELSE
                     IOFF = IZERO
                     DO 40 I = 1, IZERO - 1
                        A( IOFF ) = ZERO
                        IOFF = IOFF + N - I
   40                CONTINUE
                     IOFF = IOFF - IZERO
                     DO 50 I = IZERO, N
                        A( IOFF+I ) = ZERO
   50                CONTINUE
                  END IF
               ELSE
                  IZERO = 0
               END IF
*
*              Save a copy of the matrix A in ASAV.
*
               CALL CCOPY( NPP, A, 1, ASAV, 1 )
*
               DO 120 IEQUED = 1, 2
                  EQUED = EQUEDS( IEQUED )
                  IF( IEQUED.EQ.1 ) THEN
                     NFACT = 3
                  ELSE
                     NFACT = 1
                  END IF
*
                  DO 110 IFACT = 1, NFACT
                     FACT = FACTS( IFACT )
                     PREFAC = LSAME( FACT, 'F' )
                     NOFACT = LSAME( FACT, 'N' )
                     EQUIL = LSAME( FACT, 'E' )
*
                     IF( ZEROT ) THEN
                        IF( PREFAC )
     $                     GO TO 110
                        RCONDC = ZERO
*
                     ELSE IF( .NOT.LSAME( FACT, 'N' ) ) THEN
*
*                       Compute the condition number for comparison with
*                       the value returned by CPPSVX (FACT = 'N' reuses
*                       the condition number from the previous iteration
*                          with FACT = 'F').
*
                        CALL CCOPY( NPP, ASAV, 1, AFAC, 1 )
                        IF( EQUIL .OR. IEQUED.GT.1 ) THEN
*
*                          Compute row and column scale factors to
*                          equilibrate the matrix A.
*
                           CALL CPPEQU( UPLO, N, AFAC, S, SCOND, AMAX,
     $                                  INFO )
                           IF( INFO.EQ.0 .AND. N.GT.0 ) THEN
                              IF( IEQUED.GT.1 )
     $                           SCOND = ZERO
*
*                             Equilibrate the matrix.
*
                              CALL CLAQSP( UPLO, N, AFAC, S, SCOND,
     $                                     AMAX, EQUED )
                           END IF
                        END IF
*
*                       Save the condition number of the
*                       non-equilibrated system for use in CGET04.
*
                        IF( EQUIL )
     $                     ROLDC = RCONDC
*
*                       Compute the 1-norm of A.
*
                        ANORM = CLANSP( '1', UPLO, N, AFAC, RWORK )
*
*                       Factor the matrix A.
*
                        CALL CPPTRF( UPLO, N, AFAC, INFO )
*
*                       Form the inverse of A.
*
                        CALL CCOPY( NPP, AFAC, 1, A, 1 )
                        CALL CPPTRI( UPLO, N, A, INFO )
*
*                       Compute the 1-norm condition number of A.
*
                        AINVNM = CLANSP( '1', UPLO, N, A, RWORK )
                        IF( ANORM.LE.ZERO .OR. AINVNM.LE.ZERO ) THEN
                           RCONDC = ONE
                        ELSE
                           RCONDC = ( ONE / ANORM ) / AINVNM
                        END IF
                     END IF
   60                CONTINUE
*
*                    Restore the matrix A.
*
                     CALL CCOPY( NPP, ASAV, 1, A, 1 )
*
*                    Form an exact solution and set the right hand side.
*
                     SRNAMT = 'CLARHS'
                     CALL CLARHS( PATH, XTYPE, UPLO, ' ', N, N, KL, KU,
     $                            NRHS, A, LDA, XACT, LDA, B, LDA,
     $                            ISEED, INFO )
                     XTYPE = 'C'
                     CALL CLACPY( 'Full', N, NRHS, B, LDA, BSAV, LDA )
*
                     IF( NOFACT ) THEN
*
*                       --- Test CPPSV  ---
*
*                       Compute the L*L' or U'*U factorization of the
*                       matrix and solve the system.
*
                        CALL CCOPY( NPP, A, 1, AFAC, 1 )
                        CALL CLACPY( 'Full', N, NRHS, B, LDA, X, LDA )
*
                        SRNAMT = 'CPPSV '
                        CALL CPPSV( UPLO, N, NRHS, AFAC, X, LDA, INFO )
*
*                       Check error code from CPPSV .
*
                        IF( INFO.NE.IZERO ) THEN
                           CALL ALAERH( PATH, 'CPPSV ', INFO, IZERO,
     $                                  UPLO, N, N, -1, -1, NRHS, IMAT,
     $                                  NFAIL, NERRS, NOUT )
                           GO TO 80
                        ELSE IF( INFO.NE.0 ) THEN
                           GO TO 80
                        END IF
*
*                       Reconstruct matrix from factors and compute
*                       residual.
*
                        CALL CPPT01( UPLO, N, A, AFAC, RWORK,
     $                               RESULT( 1 ) )
*
*                       Compute residual of the computed solution.
*
                        CALL CLACPY( 'Full', N, NRHS, B, LDA, WORK,
     $                               LDA )
                        CALL CPPT02( UPLO, N, NRHS, A, X, LDA, WORK,
     $                               LDA, RWORK, RESULT( 2 ) )
*
*                       Check solution from generated exact solution.
*
                        CALL CGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC,
     $                               RESULT( 3 ) )
                        NT = 3
*
*                       Print information about the tests that did not
*                       pass the threshold.
*
                        DO 70 K = 1, NT
                           IF( RESULT( K ).GE.THRESH ) THEN
                              IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 )
     $                           CALL ALADHD( NOUT, PATH )
                              WRITE( NOUT, FMT = 9999 )'CPPSV ', UPLO,
     $                           N, IMAT, K, RESULT( K )
                              NFAIL = NFAIL + 1
                           END IF
   70                   CONTINUE
                        NRUN = NRUN + NT
   80                   CONTINUE
                     END IF
*
*                    --- Test CPPSVX ---
*
                     IF( .NOT.PREFAC .AND. NPP.GT.0 )
     $                  CALL CLASET( 'Full', NPP, 1, CMPLX( ZERO ),
     $                               CMPLX( ZERO ), AFAC, NPP )
                     CALL CLASET( 'Full', N, NRHS, CMPLX( ZERO ),
     $                            CMPLX( ZERO ), X, LDA )
                     IF( IEQUED.GT.1 .AND. N.GT.0 ) THEN
*
*                       Equilibrate the matrix if FACT='F' and
*                       EQUED='Y'.
*
                        CALL CLAQSP( UPLO, N, A, S, SCOND, AMAX, EQUED )
                     END IF
*
*                    Solve the system and compute the condition number
*                    and error bounds using CPPSVX.
*
                     SRNAMT = 'CPPSVX'
                     CALL CPPSVX( FACT, UPLO, N, NRHS, A, AFAC, EQUED,
     $                            S, B, LDA, X, LDA, RCOND, RWORK,
     $                            RWORK( NRHS+1 ), WORK,
     $                            RWORK( 2*NRHS+1 ), INFO )
*
*                    Check the error code from CPPSVX.
*
                     IF( INFO.NE.IZERO ) THEN
                        CALL ALAERH( PATH, 'CPPSVX', INFO, IZERO,
     $                               FACT // UPLO, N, N, -1, -1, NRHS,
     $                               IMAT, NFAIL, NERRS, NOUT )
                        GO TO 100
                     END IF
*
                     IF( INFO.EQ.0 ) THEN
                        IF( .NOT.PREFAC ) THEN
*
*                          Reconstruct matrix from factors and compute
*                          residual.
*
                           CALL CPPT01( UPLO, N, A, AFAC,
     $                                  RWORK( 2*NRHS+1 ), RESULT( 1 ) )
                           K1 = 1
                        ELSE
                           K1 = 2
                        END IF
*
*                       Compute residual of the computed solution.
*
                        CALL CLACPY( 'Full', N, NRHS, BSAV, LDA, WORK,
     $                               LDA )
                        CALL CPPT02( UPLO, N, NRHS, ASAV, X, LDA, WORK,
     $                               LDA, RWORK( 2*NRHS+1 ),
     $                               RESULT( 2 ) )
*
*                       Check solution from generated exact solution.
*
                        IF( NOFACT .OR. ( PREFAC .AND. LSAME( EQUED,
     $                      'N' ) ) ) THEN
                           CALL CGET04( N, NRHS, X, LDA, XACT, LDA,
     $                                  RCONDC, RESULT( 3 ) )
                        ELSE
                           CALL CGET04( N, NRHS, X, LDA, XACT, LDA,
     $                                  ROLDC, RESULT( 3 ) )
                        END IF
*
*                       Check the error bounds from iterative
*                       refinement.
*
                        CALL CPPT05( UPLO, N, NRHS, ASAV, B, LDA, X,
     $                               LDA, XACT, LDA, RWORK,
     $                               RWORK( NRHS+1 ), RESULT( 4 ) )
                     ELSE
                        K1 = 6
                     END IF
*
*                    Compare RCOND from CPPSVX with the computed value
*                    in RCONDC.
*
                     RESULT( 6 ) = SGET06( RCOND, RCONDC )
*
*                    Print information about the tests that did not pass
*                    the threshold.
*
                     DO 90 K = K1, 6
                        IF( RESULT( K ).GE.THRESH ) THEN
                           IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 )
     $                        CALL ALADHD( NOUT, PATH )
                           IF( PREFAC ) THEN
                              WRITE( NOUT, FMT = 9997 )'CPPSVX', FACT,
     $                           UPLO, N, EQUED, IMAT, K, RESULT( K )
                           ELSE
                              WRITE( NOUT, FMT = 9998 )'CPPSVX', FACT,
     $                           UPLO, N, IMAT, K, RESULT( K )
                           END IF
                           NFAIL = NFAIL + 1
                        END IF
   90                CONTINUE
                     NRUN = NRUN + 7 - K1
  100                CONTINUE
  110             CONTINUE
  120          CONTINUE
  130       CONTINUE
  140    CONTINUE
  150 CONTINUE
*
*     Print a summary of the results.
*
      CALL ALASVM( PATH, NOUT, NFAIL, NRUN, NERRS )
*
 9999 FORMAT( 1X, A6, ', UPLO=''', A1, ''', N =', I5, ', type ', I1,
     $      ', test(', I1, ')=', G12.5 )
 9998 FORMAT( 1X, A6, ', FACT=''', A1, ''', UPLO=''', A1, ''', N=', I5,
     $      ', type ', I1, ', test(', I1, ')=', G12.5 )
 9997 FORMAT( 1X, A6, ', FACT=''', A1, ''', UPLO=''', A1, ''', N=', I5,
     $      ', EQUED=''', A1, ''', type ', I1, ', test(', I1, ')=',
     $      G12.5 )
      RETURN
*
*     End of CDRVPP
*
      END
