      SUBROUTINE SERRST( PATH, NUNIT )
*
*  -- LAPACK test routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1992
*
*     .. Scalar Arguments ..
      CHARACTER*3        PATH
      INTEGER            NUNIT
*     ..
*
*  Purpose
*  =======
*
*  SERRST tests the error exits for SSYTRD, SORGTR, SORMTR, SSPTRD,
*  SOPGTR, SOPMTR, SSTEQR, SSTERF, SSTEBZ, SSTEIN, SPTEQR, and SSBTRD.
*
*  Arguments
*  =========
*
*  PATH    (input) CHARACTER*3
*          The LAPACK path name for the routines to be tested.
*
*  NUNIT   (input) INTEGER
*          The unit number for output.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NMAX, LIW, LW
      PARAMETER          ( NMAX = 3, LIW = 3*NMAX, LW = 4*NMAX )
*     ..
*     .. Local Scalars ..
      CHARACTER*2        C2
      INTEGER            I, IL, INFO, IU, J, M, NSPLIT, NT
      REAL               TOL, VL, VU
*     ..
*     .. Local Arrays ..
      INTEGER            I1( NMAX ), I2( NMAX ), I3( NMAX ), IW( LIW )
      REAL               A( NMAX, NMAX ), C( NMAX, NMAX ), D( NMAX ),
     $                   E( NMAX ), TAU( NMAX ), W( LW ), X( NMAX ),
     $                   Z( NMAX, NMAX )
*     ..
*     .. External Functions ..
      LOGICAL            LSAMEN
      EXTERNAL           LSAMEN
*     ..
*     .. External Subroutines ..
      EXTERNAL           CHKXER, SOPGTR, SOPMTR, SORGTR, SORMTR, SPTEQR,
     $                   SSBTRD, SSPTRD, SSTEBZ, SSTEIN, SSTEQR, SSTERF,
     $                   SSYTRD
*     ..
*     .. Scalars in Common ..
      LOGICAL            LERR, OK
      CHARACTER*6        SRNAMT
      INTEGER            INFOT, NOUT
*     ..
*     .. Common blocks ..
      COMMON             / INFOC / INFOT, NOUT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          REAL
*     ..
*     .. Executable Statements ..
*
      NOUT = NUNIT
      WRITE( NOUT, FMT = * )
      C2 = PATH( 2: 3 )
*
*     Set the variables to innocuous values.
*
      DO 20 J = 1, NMAX
         DO 10 I = 1, NMAX
            A( I, J ) = 1. / REAL( I+J )
   10    CONTINUE
   20 CONTINUE
      OK = .TRUE.
      NT = 0
*
*     Test error exits for the ST path.
*
      IF( LSAMEN( 2, C2, 'ST' ) ) THEN
*
*        SSYTRD
*
         SRNAMT = 'SSYTRD'
         INFOT = 1
         CALL SSYTRD( '/', 0, A, 1, D, E, TAU, W, 1, INFO )
         CALL CHKXER( 'SSYTRD', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SSYTRD( 'U', -1, A, 1, D, E, TAU, W, 1, INFO )
         CALL CHKXER( 'SSYTRD', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL SSYTRD( 'U', 2, A, 1, D, E, TAU, W, 1, INFO )
         CALL CHKXER( 'SSYTRD', INFOT, NOUT, LERR, OK )
         NT = NT + 3
*
*        SORGTR
*
         SRNAMT = 'SORGTR'
         INFOT = 1
         CALL SORGTR( '/', 0, A, 1, TAU, W, 1, INFO )
         CALL CHKXER( 'SORGTR', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SORGTR( 'U', -1, A, 1, TAU, W, 1, INFO )
         CALL CHKXER( 'SORGTR', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL SORGTR( 'U', 2, A, 1, TAU, W, 1, INFO )
         CALL CHKXER( 'SORGTR', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL SORGTR( 'U', 3, A, 3, TAU, W, 1, INFO )
         CALL CHKXER( 'SORGTR', INFOT, NOUT, LERR, OK )
         NT = NT + 4
*
*        SORMTR
*
         SRNAMT = 'SORMTR'
         INFOT = 1
         CALL SORMTR( '/', 'U', 'N', 0, 0, A, 1, TAU, C, 1, W, 1, INFO )
         CALL CHKXER( 'SORMTR', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SORMTR( 'L', '/', 'N', 0, 0, A, 1, TAU, C, 1, W, 1, INFO )
         CALL CHKXER( 'SORMTR', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL SORMTR( 'L', 'U', '/', 0, 0, A, 1, TAU, C, 1, W, 1, INFO )
         CALL CHKXER( 'SORMTR', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL SORMTR( 'L', 'U', 'N', -1, 0, A, 1, TAU, C, 1, W, 1,
     $                INFO )
         CALL CHKXER( 'SORMTR', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL SORMTR( 'L', 'U', 'N', 0, -1, A, 1, TAU, C, 1, W, 1,
     $                INFO )
         CALL CHKXER( 'SORMTR', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL SORMTR( 'L', 'U', 'N', 2, 0, A, 1, TAU, C, 2, W, 1, INFO )
         CALL CHKXER( 'SORMTR', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL SORMTR( 'R', 'U', 'N', 0, 2, A, 1, TAU, C, 1, W, 1, INFO )
         CALL CHKXER( 'SORMTR', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL SORMTR( 'L', 'U', 'N', 2, 0, A, 2, TAU, C, 1, W, 1, INFO )
         CALL CHKXER( 'SORMTR', INFOT, NOUT, LERR, OK )
         INFOT = 12
         CALL SORMTR( 'L', 'U', 'N', 0, 2, A, 1, TAU, C, 1, W, 1, INFO )
         CALL CHKXER( 'SORMTR', INFOT, NOUT, LERR, OK )
         INFOT = 12
         CALL SORMTR( 'R', 'U', 'N', 2, 0, A, 1, TAU, C, 2, W, 1, INFO )
         CALL CHKXER( 'SORMTR', INFOT, NOUT, LERR, OK )
         NT = NT + 10
*
*        SSPTRD
*
         SRNAMT = 'SSPTRD'
         INFOT = 1
         CALL SSPTRD( '/', 0, A, D, E, TAU, INFO )
         CALL CHKXER( 'SSPTRD', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SSPTRD( 'U', -1, A, D, E, TAU, INFO )
         CALL CHKXER( 'SSPTRD', INFOT, NOUT, LERR, OK )
         NT = NT + 2
*
*        SOPGTR
*
         SRNAMT = 'SOPGTR'
         INFOT = 1
         CALL SOPGTR( '/', 0, A, TAU, Z, 1, W, INFO )
         CALL CHKXER( 'SOPGTR', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SOPGTR( 'U', -1, A, TAU, Z, 1, W, INFO )
         CALL CHKXER( 'SOPGTR', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL SOPGTR( 'U', 2, A, TAU, Z, 1, W, INFO )
         CALL CHKXER( 'SOPGTR', INFOT, NOUT, LERR, OK )
         NT = NT + 3
*
*        SOPMTR
*
         SRNAMT = 'SOPMTR'
         INFOT = 1
         CALL SOPMTR( '/', 'U', 'N', 0, 0, A, TAU, C, 1, W, INFO )
         CALL CHKXER( 'SOPMTR', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SOPMTR( 'L', '/', 'N', 0, 0, A, TAU, C, 1, W, INFO )
         CALL CHKXER( 'SOPMTR', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL SOPMTR( 'L', 'U', '/', 0, 0, A, TAU, C, 1, W, INFO )
         CALL CHKXER( 'SOPMTR', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL SOPMTR( 'L', 'U', 'N', -1, 0, A, TAU, C, 1, W, INFO )
         CALL CHKXER( 'SOPMTR', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL SOPMTR( 'L', 'U', 'N', 0, -1, A, TAU, C, 1, W, INFO )
         CALL CHKXER( 'SOPMTR', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL SOPMTR( 'L', 'U', 'N', 2, 0, A, TAU, C, 1, W, INFO )
         CALL CHKXER( 'SOPMTR', INFOT, NOUT, LERR, OK )
         NT = NT + 6
*
*        SPTEQR
*
         SRNAMT = 'SPTEQR'
         INFOT = 1
         CALL SPTEQR( '/', 0, D, E, Z, 1, W, INFO )
         CALL CHKXER( 'SPTEQR', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SPTEQR( 'N', -1, D, E, Z, 1, W, INFO )
         CALL CHKXER( 'SPTEQR', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL SPTEQR( 'V', 2, D, E, Z, 1, W, INFO )
         CALL CHKXER( 'SPTEQR', INFOT, NOUT, LERR, OK )
         NT = NT + 3
*
*        SSTEBZ
*
         SRNAMT = 'SSTEBZ'
         INFOT = 1
         CALL SSTEBZ( '/', 'E', 0, VL, VU, IL, IU, TOL, D, E, M, NSPLIT,
     $                X, I1, I2, W, IW, INFO )
         CALL CHKXER( 'SSTEBZ', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SSTEBZ( 'A', '/', 0, VL, VU, IL, IU, TOL, D, E, M, NSPLIT,
     $                X, I1, I2, W, IW, INFO )
         CALL CHKXER( 'SSTEBZ', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL SSTEBZ( 'A', 'E', -1, VL, VU, IL, IU, TOL, D, E, M,
     $                NSPLIT, X, I1, I2, W, IW, INFO )
         CALL CHKXER( 'SSTEBZ', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL SSTEBZ( 'V', 'E', 0, 0.0, 0.0, IL, IU, TOL, D, E, M,
     $                NSPLIT, X, I1, I2, W, IW, INFO )
         CALL CHKXER( 'SSTEBZ', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL SSTEBZ( 'I', 'E', 0, VL, VU, 0, 0, TOL, D, E, M, NSPLIT,
     $                X, I1, I2, W, IW, INFO )
         CALL CHKXER( 'SSTEBZ', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL SSTEBZ( 'I', 'E', 1, VL, VU, 1, 0, TOL, D, E, M, NSPLIT,
     $                X, I1, I2, W, IW, INFO )
         CALL CHKXER( 'SSTEBZ', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL SSTEBZ( 'I', 'E', 0, VL, VU, 1, 1, TOL, D, E, M, NSPLIT,
     $                X, I1, I2, W, IW, INFO )
         CALL CHKXER( 'SSTEBZ', INFOT, NOUT, LERR, OK )
         NT = NT + 7
*
*        SSTEIN
*
         SRNAMT = 'SSTEIN'
         INFOT = 1
         CALL SSTEIN( -1, D, E, 0, X, I1, I2, Z, 1, W, IW, I3, INFO )
         CALL CHKXER( 'SSTEIN', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL SSTEIN( 0, D, E, -1, X, I1, I2, Z, 1, W, IW, I3, INFO )
         CALL CHKXER( 'SSTEIN', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL SSTEIN( 0, D, E, 1, X, I1, I2, Z, 1, W, IW, I3, INFO )
         CALL CHKXER( 'SSTEIN', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL SSTEIN( 2, D, E, 0, X, I1, I2, Z, 1, W, IW, I3, INFO )
         CALL CHKXER( 'SSTEIN', INFOT, NOUT, LERR, OK )
         NT = NT + 4
*
*        SSTEQR
*
         SRNAMT = 'SSTEQR'
         INFOT = 1
         CALL SSTEQR( '/', 0, D, E, Z, 1, W, INFO )
         CALL CHKXER( 'SSTEQR', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SSTEQR( 'N', -1, D, E, Z, 1, W, INFO )
         CALL CHKXER( 'SSTEQR', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL SSTEQR( 'V', 2, D, E, Z, 1, W, INFO )
         CALL CHKXER( 'SSTEQR', INFOT, NOUT, LERR, OK )
         NT = NT + 3
*
*        SSTERF
*
         SRNAMT = 'SSTERF'
         INFOT = 1
         CALL SSTERF( -1, D, E, INFO )
         CALL CHKXER( 'SSTERF', INFOT, NOUT, LERR, OK )
         NT = NT + 1
*
*     Test error exits for the SB path.
*
      ELSE IF( LSAMEN( 2, C2, 'SB' ) ) THEN
*
*        SSBTRD
*
         SRNAMT = 'SSBTRD'
         INFOT = 1
         CALL SSBTRD( '/', 'U', 0, 0, A, 1, D, E, Z, 1, W, INFO )
         CALL CHKXER( 'SSBTRD', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SSBTRD( 'N', '/', 0, 0, A, 1, D, E, Z, 1, W, INFO )
         CALL CHKXER( 'SSBTRD', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL SSBTRD( 'N', 'U', -1, 0, A, 1, D, E, Z, 1, W, INFO )
         CALL CHKXER( 'SSBTRD', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL SSBTRD( 'N', 'U', 0, -1, A, 1, D, E, Z, 1, W, INFO )
         CALL CHKXER( 'SSBTRD', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL SSBTRD( 'N', 'U', 1, 1, A, 1, D, E, Z, 1, W, INFO )
         CALL CHKXER( 'SSBTRD', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL SSBTRD( 'V', 'U', 2, 0, A, 1, D, E, Z, 1, W, INFO )
         CALL CHKXER( 'SSBTRD', INFOT, NOUT, LERR, OK )
         NT = NT + 6
      END IF
*
*     Print a summary line.
*
      IF( OK ) THEN
         WRITE( NOUT, FMT = 9999 )PATH, NT
      ELSE
         WRITE( NOUT, FMT = 9998 )PATH
      END IF
*
 9999 FORMAT( 1X, A3, ' routines passed the tests of the error exits',
     $      ' (', I3, ' tests done)' )
 9998 FORMAT( ' *** ', A3, ' routines failed the tests of the error ',
     $      'exits ***' )
*
      RETURN
*
*     End of SERRST
*
      END
