      SUBROUTINE CERRST( PATH, NUNIT )
*
*  -- LAPACK test routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER*3        PATH
      INTEGER            NUNIT
*     ..
*
*  Purpose
*  =======
*
*  CERRST tests the error exits for CHETRD, CUNGTR, CUNMTR, CHPTRD,
*  CUPGTR, CUPMTR, CSTEQR, CSTEIN, CPTEQR, and CHBTRD.
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
      INTEGER            NMAX, LW
      PARAMETER          ( NMAX = 3, LW = 5*NMAX )
*     ..
*     .. Local Scalars ..
      CHARACTER*2        C2
      INTEGER            I, INFO, J, NT
*     ..
*     .. Local Arrays ..
      INTEGER            I1( NMAX ), I2( NMAX ), I3( NMAX ), IW( LW )
      REAL               D( NMAX ), E( NMAX ), RW( LW ), X( NMAX )
      COMPLEX            A( NMAX, NMAX ), C( NMAX, NMAX ), TAU( NMAX ),
     $                   W( LW ), Z( NMAX, NMAX )
*     ..
*     .. External Functions ..
      LOGICAL            LSAMEN
      EXTERNAL           LSAMEN
*     ..
*     .. External Subroutines ..
      EXTERNAL           CHBTRD, CHETRD, CHKXER, CHPTRD, CPTEQR, CSTEIN,
     $                   CSTEQR, CUNGTR, CUNMTR, CUPGTR, CUPMTR
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
*        CHETRD
*
         SRNAMT = 'CHETRD'
         INFOT = 1
         CALL CHETRD( '/', 0, A, 1, D, E, TAU, W, 1, INFO )
         CALL CHKXER( 'CHETRD', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CHETRD( 'U', -1, A, 1, D, E, TAU, W, 1, INFO )
         CALL CHKXER( 'CHETRD', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL CHETRD( 'U', 2, A, 1, D, E, TAU, W, 1, INFO )
         CALL CHKXER( 'CHETRD', INFOT, NOUT, LERR, OK )
         NT = NT + 3
*
*        CUNGTR
*
         SRNAMT = 'CUNGTR'
         INFOT = 1
         CALL CUNGTR( '/', 0, A, 1, TAU, W, 1, INFO )
         CALL CHKXER( 'CUNGTR', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CUNGTR( 'U', -1, A, 1, TAU, W, 1, INFO )
         CALL CHKXER( 'CUNGTR', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL CUNGTR( 'U', 2, A, 1, TAU, W, 1, INFO )
         CALL CHKXER( 'CUNGTR', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL CUNGTR( 'U', 3, A, 3, TAU, W, 1, INFO )
         CALL CHKXER( 'CUNGTR', INFOT, NOUT, LERR, OK )
         NT = NT + 4
*
*        CUNMTR
*
         SRNAMT = 'CUNMTR'
         INFOT = 1
         CALL CUNMTR( '/', 'U', 'N', 0, 0, A, 1, TAU, C, 1, W, 1, INFO )
         CALL CHKXER( 'CUNMTR', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CUNMTR( 'L', '/', 'N', 0, 0, A, 1, TAU, C, 1, W, 1, INFO )
         CALL CHKXER( 'CUNMTR', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL CUNMTR( 'L', 'U', '/', 0, 0, A, 1, TAU, C, 1, W, 1, INFO )
         CALL CHKXER( 'CUNMTR', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL CUNMTR( 'L', 'U', 'N', -1, 0, A, 1, TAU, C, 1, W, 1,
     $                INFO )
         CALL CHKXER( 'CUNMTR', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL CUNMTR( 'L', 'U', 'N', 0, -1, A, 1, TAU, C, 1, W, 1,
     $                INFO )
         CALL CHKXER( 'CUNMTR', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL CUNMTR( 'L', 'U', 'N', 2, 0, A, 1, TAU, C, 2, W, 1, INFO )
         CALL CHKXER( 'CUNMTR', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL CUNMTR( 'R', 'U', 'N', 0, 2, A, 1, TAU, C, 1, W, 1, INFO )
         CALL CHKXER( 'CUNMTR', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL CUNMTR( 'L', 'U', 'N', 2, 0, A, 2, TAU, C, 1, W, 1, INFO )
         CALL CHKXER( 'CUNMTR', INFOT, NOUT, LERR, OK )
         INFOT = 12
         CALL CUNMTR( 'L', 'U', 'N', 0, 2, A, 1, TAU, C, 1, W, 1, INFO )
         CALL CHKXER( 'CUNMTR', INFOT, NOUT, LERR, OK )
         INFOT = 12
         CALL CUNMTR( 'R', 'U', 'N', 2, 0, A, 1, TAU, C, 2, W, 1, INFO )
         CALL CHKXER( 'CUNMTR', INFOT, NOUT, LERR, OK )
         NT = NT + 10
*
*        CHPTRD
*
         SRNAMT = 'CHPTRD'
         INFOT = 1
         CALL CHPTRD( '/', 0, A, D, E, TAU, INFO )
         CALL CHKXER( 'CHPTRD', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CHPTRD( 'U', -1, A, D, E, TAU, INFO )
         CALL CHKXER( 'CHPTRD', INFOT, NOUT, LERR, OK )
         NT = NT + 2
*
*        CUPGTR
*
         SRNAMT = 'CUPGTR'
         INFOT = 1
         CALL CUPGTR( '/', 0, A, TAU, Z, 1, W, INFO )
         CALL CHKXER( 'CUPGTR', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CUPGTR( 'U', -1, A, TAU, Z, 1, W, INFO )
         CALL CHKXER( 'CUPGTR', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL CUPGTR( 'U', 2, A, TAU, Z, 1, W, INFO )
         CALL CHKXER( 'CUPGTR', INFOT, NOUT, LERR, OK )
         NT = NT + 3
*
*        CUPMTR
*
         SRNAMT = 'CUPMTR'
         INFOT = 1
         CALL CUPMTR( '/', 'U', 'N', 0, 0, A, TAU, C, 1, W, INFO )
         CALL CHKXER( 'CUPMTR', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CUPMTR( 'L', '/', 'N', 0, 0, A, TAU, C, 1, W, INFO )
         CALL CHKXER( 'CUPMTR', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL CUPMTR( 'L', 'U', '/', 0, 0, A, TAU, C, 1, W, INFO )
         CALL CHKXER( 'CUPMTR', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL CUPMTR( 'L', 'U', 'N', -1, 0, A, TAU, C, 1, W, INFO )
         CALL CHKXER( 'CUPMTR', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL CUPMTR( 'L', 'U', 'N', 0, -1, A, TAU, C, 1, W, INFO )
         CALL CHKXER( 'CUPMTR', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL CUPMTR( 'L', 'U', 'N', 2, 0, A, TAU, C, 1, W, INFO )
         CALL CHKXER( 'CUPMTR', INFOT, NOUT, LERR, OK )
         NT = NT + 6
*
*        CPTEQR
*
         SRNAMT = 'CPTEQR'
         INFOT = 1
         CALL CPTEQR( '/', 0, D, E, Z, 1, RW, INFO )
         CALL CHKXER( 'CPTEQR', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CPTEQR( 'N', -1, D, E, Z, 1, RW, INFO )
         CALL CHKXER( 'CPTEQR', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL CPTEQR( 'V', 2, D, E, Z, 1, RW, INFO )
         CALL CHKXER( 'CPTEQR', INFOT, NOUT, LERR, OK )
         NT = NT + 3
*
*        CSTEIN
*
         SRNAMT = 'CSTEIN'
         INFOT = 1
         CALL CSTEIN( -1, D, E, 0, X, I1, I2, Z, 1, RW, IW, I3, INFO )
         CALL CHKXER( 'CSTEIN', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL CSTEIN( 0, D, E, -1, X, I1, I2, Z, 1, RW, IW, I3, INFO )
         CALL CHKXER( 'CSTEIN', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL CSTEIN( 0, D, E, 1, X, I1, I2, Z, 1, RW, IW, I3, INFO )
         CALL CHKXER( 'CSTEIN', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL CSTEIN( 2, D, E, 0, X, I1, I2, Z, 1, RW, IW, I3, INFO )
         CALL CHKXER( 'CSTEIN', INFOT, NOUT, LERR, OK )
         NT = NT + 4
*
*        CSTEQR
*
         SRNAMT = 'CSTEQR'
         INFOT = 1
         CALL CSTEQR( '/', 0, D, E, Z, 1, RW, INFO )
         CALL CHKXER( 'CSTEQR', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CSTEQR( 'N', -1, D, E, Z, 1, RW, INFO )
         CALL CHKXER( 'CSTEQR', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL CSTEQR( 'V', 2, D, E, Z, 1, RW, INFO )
         CALL CHKXER( 'CSTEQR', INFOT, NOUT, LERR, OK )
         NT = NT + 3
*
*     Test error exits for the HB path.
*
      ELSE IF( LSAMEN( 2, C2, 'HB' ) ) THEN
*
*        CHBTRD
*
         SRNAMT = 'CHBTRD'
         INFOT = 1
         CALL CHBTRD( '/', 'U', 0, 0, A, 1, D, E, Z, 1, W, INFO )
         CALL CHKXER( 'CHBTRD', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CHBTRD( 'N', '/', 0, 0, A, 1, D, E, Z, 1, W, INFO )
         CALL CHKXER( 'CHBTRD', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL CHBTRD( 'N', 'U', -1, 0, A, 1, D, E, Z, 1, W, INFO )
         CALL CHKXER( 'CHBTRD', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL CHBTRD( 'N', 'U', 0, -1, A, 1, D, E, Z, 1, W, INFO )
         CALL CHKXER( 'CHBTRD', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL CHBTRD( 'N', 'U', 1, 1, A, 1, D, E, Z, 1, W, INFO )
         CALL CHKXER( 'CHBTRD', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL CHBTRD( 'V', 'U', 2, 0, A, 1, D, E, Z, 1, W, INFO )
         CALL CHKXER( 'CHBTRD', INFOT, NOUT, LERR, OK )
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
*     End of CERRST
*
      END
