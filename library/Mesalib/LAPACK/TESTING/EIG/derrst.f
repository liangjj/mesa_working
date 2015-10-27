      SUBROUTINE DERRST( PATH, NUNIT )
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
*  DERRST tests the error exits for DSYTRD, DORGTR, DORMTR, DSPTRD,
*  DOPGTR, DOPMTR, DSTEQR, SSTERF, SSTEBZ, SSTEIN, DPTEQR, and DSBTRD.
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
      DOUBLE PRECISION   TOL, VL, VU
*     ..
*     .. Local Arrays ..
      INTEGER            I1( NMAX ), I2( NMAX ), I3( NMAX ), IW( LIW )
      DOUBLE PRECISION   A( NMAX, NMAX ), C( NMAX, NMAX ), D( NMAX ),
     $                   E( NMAX ), TAU( NMAX ), W( LW ), X( NMAX ),
     $                   Z( NMAX, NMAX )
*     ..
*     .. External Functions ..
      LOGICAL            LSAMEN
      EXTERNAL           LSAMEN
*     ..
*     .. External Subroutines ..
      EXTERNAL           CHKXER, DOPGTR, DOPMTR, DORGTR, DORMTR, DPTEQR,
     $                   DSBTRD, DSPTRD, DSTEBZ, DSTEIN, DSTEQR, DSTERF,
     $                   DSYTRD
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
      INTRINSIC          DBLE
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
            A( I, J ) = 1.D0 / DBLE( I+J )
   10    CONTINUE
   20 CONTINUE
      OK = .TRUE.
      NT = 0
*
*     Test error exits for the ST path.
*
      IF( LSAMEN( 2, C2, 'ST' ) ) THEN
*
*        DSYTRD
*
         SRNAMT = 'DSYTRD'
         INFOT = 1
         CALL DSYTRD( '/', 0, A, 1, D, E, TAU, W, 1, INFO )
         CALL CHKXER( 'DSYTRD', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DSYTRD( 'U', -1, A, 1, D, E, TAU, W, 1, INFO )
         CALL CHKXER( 'DSYTRD', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL DSYTRD( 'U', 2, A, 1, D, E, TAU, W, 1, INFO )
         CALL CHKXER( 'DSYTRD', INFOT, NOUT, LERR, OK )
         NT = NT + 3
*
*        DORGTR
*
         SRNAMT = 'DORGTR'
         INFOT = 1
         CALL DORGTR( '/', 0, A, 1, TAU, W, 1, INFO )
         CALL CHKXER( 'DORGTR', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DORGTR( 'U', -1, A, 1, TAU, W, 1, INFO )
         CALL CHKXER( 'DORGTR', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL DORGTR( 'U', 2, A, 1, TAU, W, 1, INFO )
         CALL CHKXER( 'DORGTR', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL DORGTR( 'U', 3, A, 3, TAU, W, 1, INFO )
         CALL CHKXER( 'DORGTR', INFOT, NOUT, LERR, OK )
         NT = NT + 4
*
*        DORMTR
*
         SRNAMT = 'DORMTR'
         INFOT = 1
         CALL DORMTR( '/', 'U', 'N', 0, 0, A, 1, TAU, C, 1, W, 1, INFO )
         CALL CHKXER( 'DORMTR', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DORMTR( 'L', '/', 'N', 0, 0, A, 1, TAU, C, 1, W, 1, INFO )
         CALL CHKXER( 'DORMTR', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DORMTR( 'L', 'U', '/', 0, 0, A, 1, TAU, C, 1, W, 1, INFO )
         CALL CHKXER( 'DORMTR', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL DORMTR( 'L', 'U', 'N', -1, 0, A, 1, TAU, C, 1, W, 1,
     $                INFO )
         CALL CHKXER( 'DORMTR', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL DORMTR( 'L', 'U', 'N', 0, -1, A, 1, TAU, C, 1, W, 1,
     $                INFO )
         CALL CHKXER( 'DORMTR', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL DORMTR( 'L', 'U', 'N', 2, 0, A, 1, TAU, C, 2, W, 1, INFO )
         CALL CHKXER( 'DORMTR', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL DORMTR( 'R', 'U', 'N', 0, 2, A, 1, TAU, C, 1, W, 1, INFO )
         CALL CHKXER( 'DORMTR', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL DORMTR( 'L', 'U', 'N', 2, 0, A, 2, TAU, C, 1, W, 1, INFO )
         CALL CHKXER( 'DORMTR', INFOT, NOUT, LERR, OK )
         INFOT = 12
         CALL DORMTR( 'L', 'U', 'N', 0, 2, A, 1, TAU, C, 1, W, 1, INFO )
         CALL CHKXER( 'DORMTR', INFOT, NOUT, LERR, OK )
         INFOT = 12
         CALL DORMTR( 'R', 'U', 'N', 2, 0, A, 1, TAU, C, 2, W, 1, INFO )
         CALL CHKXER( 'DORMTR', INFOT, NOUT, LERR, OK )
         NT = NT + 10
*
*        DSPTRD
*
         SRNAMT = 'DSPTRD'
         INFOT = 1
         CALL DSPTRD( '/', 0, A, D, E, TAU, INFO )
         CALL CHKXER( 'DSPTRD', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DSPTRD( 'U', -1, A, D, E, TAU, INFO )
         CALL CHKXER( 'DSPTRD', INFOT, NOUT, LERR, OK )
         NT = NT + 2
*
*        DOPGTR
*
         SRNAMT = 'DOPGTR'
         INFOT = 1
         CALL DOPGTR( '/', 0, A, TAU, Z, 1, W, INFO )
         CALL CHKXER( 'DOPGTR', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DOPGTR( 'U', -1, A, TAU, Z, 1, W, INFO )
         CALL CHKXER( 'DOPGTR', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL DOPGTR( 'U', 2, A, TAU, Z, 1, W, INFO )
         CALL CHKXER( 'DOPGTR', INFOT, NOUT, LERR, OK )
         NT = NT + 3
*
*        DOPMTR
*
         SRNAMT = 'DOPMTR'
         INFOT = 1
         CALL DOPMTR( '/', 'U', 'N', 0, 0, A, TAU, C, 1, W, INFO )
         CALL CHKXER( 'DOPMTR', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DOPMTR( 'L', '/', 'N', 0, 0, A, TAU, C, 1, W, INFO )
         CALL CHKXER( 'DOPMTR', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DOPMTR( 'L', 'U', '/', 0, 0, A, TAU, C, 1, W, INFO )
         CALL CHKXER( 'DOPMTR', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL DOPMTR( 'L', 'U', 'N', -1, 0, A, TAU, C, 1, W, INFO )
         CALL CHKXER( 'DOPMTR', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL DOPMTR( 'L', 'U', 'N', 0, -1, A, TAU, C, 1, W, INFO )
         CALL CHKXER( 'DOPMTR', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL DOPMTR( 'L', 'U', 'N', 2, 0, A, TAU, C, 1, W, INFO )
         CALL CHKXER( 'DOPMTR', INFOT, NOUT, LERR, OK )
         NT = NT + 6
*
*        DPTEQR
*
         SRNAMT = 'DPTEQR'
         INFOT = 1
         CALL DPTEQR( '/', 0, D, E, Z, 1, W, INFO )
         CALL CHKXER( 'DPTEQR', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DPTEQR( 'N', -1, D, E, Z, 1, W, INFO )
         CALL CHKXER( 'DPTEQR', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL DPTEQR( 'V', 2, D, E, Z, 1, W, INFO )
         CALL CHKXER( 'DPTEQR', INFOT, NOUT, LERR, OK )
         NT = NT + 3
*
*        DSTEBZ
*
         SRNAMT = 'DSTEBZ'
         INFOT = 1
         CALL DSTEBZ( '/', 'E', 0, VL, VU, IL, IU, TOL, D, E, M, NSPLIT,
     $                X, I1, I2, W, IW, INFO )
         CALL CHKXER( 'DSTEBZ', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DSTEBZ( 'A', '/', 0, VL, VU, IL, IU, TOL, D, E, M, NSPLIT,
     $                X, I1, I2, W, IW, INFO )
         CALL CHKXER( 'DSTEBZ', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DSTEBZ( 'A', 'E', -1, VL, VU, IL, IU, TOL, D, E, M,
     $                NSPLIT, X, I1, I2, W, IW, INFO )
         CALL CHKXER( 'DSTEBZ', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL DSTEBZ( 'V', 'E', 0, 0.0D0, 0.0D0, IL, IU, TOL, D, E, M,
     $                NSPLIT, X, I1, I2, W, IW, INFO )
         CALL CHKXER( 'DSTEBZ', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL DSTEBZ( 'I', 'E', 0, VL, VU, 0, 0, TOL, D, E, M, NSPLIT,
     $                X, I1, I2, W, IW, INFO )
         CALL CHKXER( 'DSTEBZ', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL DSTEBZ( 'I', 'E', 1, VL, VU, 1, 0, TOL, D, E, M, NSPLIT,
     $                X, I1, I2, W, IW, INFO )
         CALL CHKXER( 'DSTEBZ', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL DSTEBZ( 'I', 'E', 0, VL, VU, 1, 1, TOL, D, E, M, NSPLIT,
     $                X, I1, I2, W, IW, INFO )
         CALL CHKXER( 'DSTEBZ', INFOT, NOUT, LERR, OK )
         NT = NT + 7
*
*        DSTEIN
*
         SRNAMT = 'DSTEIN'
         INFOT = 1
         CALL DSTEIN( -1, D, E, 0, X, I1, I2, Z, 1, W, IW, I3, INFO )
         CALL CHKXER( 'DSTEIN', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL DSTEIN( 0, D, E, -1, X, I1, I2, Z, 1, W, IW, I3, INFO )
         CALL CHKXER( 'DSTEIN', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL DSTEIN( 0, D, E, 1, X, I1, I2, Z, 1, W, IW, I3, INFO )
         CALL CHKXER( 'DSTEIN', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL DSTEIN( 2, D, E, 0, X, I1, I2, Z, 1, W, IW, I3, INFO )
         CALL CHKXER( 'DSTEIN', INFOT, NOUT, LERR, OK )
         NT = NT + 4
*
*        DSTEQR
*
         SRNAMT = 'DSTEQR'
         INFOT = 1
         CALL DSTEQR( '/', 0, D, E, Z, 1, W, INFO )
         CALL CHKXER( 'DSTEQR', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DSTEQR( 'N', -1, D, E, Z, 1, W, INFO )
         CALL CHKXER( 'DSTEQR', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL DSTEQR( 'V', 2, D, E, Z, 1, W, INFO )
         CALL CHKXER( 'DSTEQR', INFOT, NOUT, LERR, OK )
         NT = NT + 3
*
*        DSTERF
*
         SRNAMT = 'DSTERF'
         INFOT = 1
         CALL DSTERF( -1, D, E, INFO )
         CALL CHKXER( 'DSTERF', INFOT, NOUT, LERR, OK )
         NT = NT + 1
*
*     Test error exits for the SB path.
*
      ELSE IF( LSAMEN( 2, C2, 'SB' ) ) THEN
*
*        DSBTRD
*
         SRNAMT = 'DSBTRD'
         INFOT = 1
         CALL DSBTRD( '/', 'U', 0, 0, A, 1, D, E, Z, 1, W, INFO )
         CALL CHKXER( 'DSBTRD', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DSBTRD( 'N', '/', 0, 0, A, 1, D, E, Z, 1, W, INFO )
         CALL CHKXER( 'DSBTRD', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DSBTRD( 'N', 'U', -1, 0, A, 1, D, E, Z, 1, W, INFO )
         CALL CHKXER( 'DSBTRD', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL DSBTRD( 'N', 'U', 0, -1, A, 1, D, E, Z, 1, W, INFO )
         CALL CHKXER( 'DSBTRD', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL DSBTRD( 'N', 'U', 1, 1, A, 1, D, E, Z, 1, W, INFO )
         CALL CHKXER( 'DSBTRD', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL DSBTRD( 'V', 'U', 2, 0, A, 1, D, E, Z, 1, W, INFO )
         CALL CHKXER( 'DSBTRD', INFOT, NOUT, LERR, OK )
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
*     End of DERRST
*
      END
