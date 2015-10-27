      SUBROUTINE ZERRGG( PATH, NUNIT )
*
*  -- LAPACK test routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      CHARACTER*3        PATH
      INTEGER            NUNIT
*     ..
*
*  Purpose
*  =======
*
*  ZERRGG tests the error exits for ZGGHRD, ZHGEQZ, and ZTGEVC.
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
      PARAMETER          ( NMAX = 3, LW = 6*NMAX )
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      CHARACTER*2        C2
      INTEGER            I, INFO, J, M, NT
*     ..
*     .. Local Arrays ..
      LOGICAL            SEL( NMAX )
      DOUBLE PRECISION   RW( LW )
      COMPLEX*16         A( NMAX, NMAX ), ALPHA( NMAX ),
     $                   B( NMAX, NMAX ), BETA( NMAX ), Q( NMAX, NMAX ),
     $                   W( LW ), Z( NMAX, NMAX )
*     ..
*     .. External Functions ..
      LOGICAL            LSAMEN
      EXTERNAL           LSAMEN
*     ..
*     .. External Subroutines ..
      EXTERNAL           CHKXER, ZGGHRD, ZHGEQZ, ZTGEVC
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
            A( I, J ) = ZERO
            B( I, J ) = ZERO
   10    CONTINUE
   20 CONTINUE
      DO 30 I = 1, NMAX
         A( I, I ) = ONE
         B( I, I ) = ONE
   30 CONTINUE
      OK = .TRUE.
      NT = 0
*
*     Test error exits for the GG path.
*
      IF( LSAMEN( 2, C2, 'GG' ) ) THEN
*
*        ZGGHRD
*
         SRNAMT = 'ZGGHRD'
         INFOT = 1
         CALL ZGGHRD( '/', 'N', 0, 1, 0, A, 1, B, 1, Q, 1, Z, 1, INFO )
         CALL CHKXER( 'ZGGHRD', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZGGHRD( 'N', '/', 0, 1, 0, A, 1, B, 1, Q, 1, Z, 1, INFO )
         CALL CHKXER( 'ZGGHRD', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZGGHRD( 'N', 'N', -1, 0, 0, A, 1, B, 1, Q, 1, Z, 1, INFO )
         CALL CHKXER( 'ZGGHRD', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL ZGGHRD( 'N', 'N', 0, 0, 0, A, 1, B, 1, Q, 1, Z, 1, INFO )
         CALL CHKXER( 'ZGGHRD', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL ZGGHRD( 'N', 'N', 0, 1, 1, A, 1, B, 1, Q, 1, Z, 1, INFO )
         CALL CHKXER( 'ZGGHRD', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL ZGGHRD( 'N', 'N', 2, 1, 1, A, 1, B, 2, Q, 1, Z, 1, INFO )
         CALL CHKXER( 'ZGGHRD', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL ZGGHRD( 'N', 'N', 2, 1, 1, A, 2, B, 1, Q, 1, Z, 1, INFO )
         CALL CHKXER( 'ZGGHRD', INFOT, NOUT, LERR, OK )
         INFOT = 11
         CALL ZGGHRD( 'V', 'N', 2, 1, 1, A, 2, B, 2, Q, 1, Z, 1, INFO )
         CALL CHKXER( 'ZGGHRD', INFOT, NOUT, LERR, OK )
         INFOT = 13
         CALL ZGGHRD( 'N', 'V', 2, 1, 1, A, 2, B, 2, Q, 1, Z, 1, INFO )
         CALL CHKXER( 'ZGGHRD', INFOT, NOUT, LERR, OK )
         NT = NT + 9
*
*        ZHGEQZ
*
         SRNAMT = 'ZHGEQZ'
         INFOT = 1
         CALL ZHGEQZ( '/', 'N', 'N', 0, 1, 0, A, 1, B, 1, ALPHA, BETA,
     $                Q, 1, Z, 1, W, 1, RW, INFO )
         CALL CHKXER( 'ZHGEQZ', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZHGEQZ( 'E', '/', 'N', 0, 1, 0, A, 1, B, 1, ALPHA, BETA,
     $                Q, 1, Z, 1, W, 1, RW, INFO )
         CALL CHKXER( 'ZHGEQZ', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZHGEQZ( 'E', 'N', '/', 0, 1, 0, A, 1, B, 1, ALPHA, BETA,
     $                Q, 1, Z, 1, W, 1, RW, INFO )
         CALL CHKXER( 'ZHGEQZ', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL ZHGEQZ( 'E', 'N', 'N', -1, 0, 0, A, 1, B, 1, ALPHA, BETA,
     $                Q, 1, Z, 1, W, 1, RW, INFO )
         CALL CHKXER( 'ZHGEQZ', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL ZHGEQZ( 'E', 'N', 'N', 0, 0, 0, A, 1, B, 1, ALPHA, BETA,
     $                Q, 1, Z, 1, W, 1, RW, INFO )
         CALL CHKXER( 'ZHGEQZ', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL ZHGEQZ( 'E', 'N', 'N', 0, 1, 1, A, 1, B, 1, ALPHA, BETA,
     $                Q, 1, Z, 1, W, 1, RW, INFO )
         CALL CHKXER( 'ZHGEQZ', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL ZHGEQZ( 'E', 'N', 'N', 2, 1, 1, A, 1, B, 2, ALPHA, BETA,
     $                Q, 1, Z, 1, W, 1, RW, INFO )
         CALL CHKXER( 'ZHGEQZ', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL ZHGEQZ( 'E', 'N', 'N', 2, 1, 1, A, 2, B, 1, ALPHA, BETA,
     $                Q, 1, Z, 1, W, 1, RW, INFO )
         CALL CHKXER( 'ZHGEQZ', INFOT, NOUT, LERR, OK )
         INFOT = 14
         CALL ZHGEQZ( 'E', 'V', 'N', 2, 1, 1, A, 2, B, 2, ALPHA, BETA,
     $                Q, 1, Z, 1, W, 1, RW, INFO )
         CALL CHKXER( 'ZHGEQZ', INFOT, NOUT, LERR, OK )
         INFOT = 16
         CALL ZHGEQZ( 'E', 'N', 'V', 2, 1, 1, A, 2, B, 2, ALPHA, BETA,
     $                Q, 1, Z, 1, W, 1, RW, INFO )
         CALL CHKXER( 'ZHGEQZ', INFOT, NOUT, LERR, OK )
         NT = NT + 10
*
*        ZTGEVC
*
         SRNAMT = 'ZTGEVC'
         INFOT = 1
         CALL ZTGEVC( '/', 'R', SEL, 0, A, 1, B, 1, Q, 1, Z, 1, 0, M, W,
     $                RW, INFO )
         CALL CHKXER( 'ZTGEVC', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZTGEVC( 'A', '/', SEL, 0, A, 1, B, 1, Q, 1, Z, 1, 0, M, W,
     $                RW, INFO )
         CALL CHKXER( 'ZTGEVC', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL ZTGEVC( 'A', 'R', SEL, -1, A, 1, B, 1, Q, 1, Z, 1, 0, M,
     $                W, RW, INFO )
         CALL CHKXER( 'ZTGEVC', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL ZTGEVC( 'A', 'R', SEL, 2, A, 1, B, 2, Q, 1, Z, 2, 0, M, W,
     $                RW, INFO )
         CALL CHKXER( 'ZTGEVC', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL ZTGEVC( 'A', 'R', SEL, 2, A, 2, B, 1, Q, 1, Z, 2, 0, M, W,
     $                RW, INFO )
         CALL CHKXER( 'ZTGEVC', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL ZTGEVC( 'A', 'L', SEL, 2, A, 2, B, 2, Q, 1, Z, 1, 0, M, W,
     $                RW, INFO )
         CALL CHKXER( 'ZTGEVC', INFOT, NOUT, LERR, OK )
         INFOT = 12
         CALL ZTGEVC( 'A', 'R', SEL, 2, A, 2, B, 2, Q, 1, Z, 1, 0, M, W,
     $                RW, INFO )
         CALL CHKXER( 'ZTGEVC', INFOT, NOUT, LERR, OK )
         INFOT = 13
         CALL ZTGEVC( 'A', 'R', SEL, 2, A, 2, B, 2, Q, 1, Z, 2, 1, M, W,
     $                RW, INFO )
         CALL CHKXER( 'ZTGEVC', INFOT, NOUT, LERR, OK )
         NT = NT + 8
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
 9999 FORMAT( 1X, A3, ' routines passed the tests of the error exits (',
     $      I3, ' tests done)' )
 9998 FORMAT( ' *** ', A3, ' routines failed the tests of the error ',
     $      'exits ***' )
*
      RETURN
*
*     End of ZERRGG
*
      END
