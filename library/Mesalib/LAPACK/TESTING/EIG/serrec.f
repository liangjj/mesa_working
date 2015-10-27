      SUBROUTINE SERREC( PATH, NUNIT )
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
*  SERREC tests the error exits for the routines for eigen- condition
*  estimation for REAL matrices:
*     STRSYL, STREXC, STRSNA and STRSEN.
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
      INTEGER            NMAX
      REAL               ONE, ZERO
      PARAMETER          ( NMAX = 4, ONE = 1.0E0, ZERO = 0.0E0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, INFO, J, M, NT
      REAL               SCALE
*     ..
*     .. Local Arrays ..
      LOGICAL            SEL( NMAX )
      INTEGER            IWORK( NMAX )
      REAL               A( NMAX, NMAX ), B( NMAX, NMAX ),
     $                   C( NMAX, NMAX ), S( NMAX ), SEP( NMAX ),
     $                   WI( NMAX ), WORK( NMAX ), WR( NMAX )
*     ..
*     .. External Subroutines ..
      EXTERNAL           CHKXER, STREXC, STRSEN, STRSNA, STRSYL
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
      OK = .TRUE.
      NT = 0
*
*     Initialize A, B and SEL
*
      DO 20 J = 1, NMAX
         DO 10 I = 1, NMAX
            A( I, J ) = ZERO
            B( I, J ) = ZERO
   10    CONTINUE
   20 CONTINUE
      DO 30 I = 1, NMAX
         A( I, I ) = ONE
         SEL( I ) = .TRUE.
   30 CONTINUE
*
*     Test STRSYL
*
      SRNAMT = 'STRSYL'
      INFOT = 1
      CALL STRSYL( 'X', 'N', 1, 0, 0, A, 1, B, 1, C, 1, SCALE, INFO )
      CALL CHKXER( 'STRSYL', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL STRSYL( 'N', 'X', 1, 0, 0, A, 1, B, 1, C, 1, SCALE, INFO )
      CALL CHKXER( 'STRSYL', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL STRSYL( 'N', 'N', 0, 0, 0, A, 1, B, 1, C, 1, SCALE, INFO )
      CALL CHKXER( 'STRSYL', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL STRSYL( 'N', 'N', 1, -1, 0, A, 1, B, 1, C, 1, SCALE, INFO )
      CALL CHKXER( 'STRSYL', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL STRSYL( 'N', 'N', 1, 0, -1, A, 1, B, 1, C, 1, SCALE, INFO )
      CALL CHKXER( 'STRSYL', INFOT, NOUT, LERR, OK )
      INFOT = 7
      CALL STRSYL( 'N', 'N', 1, 2, 0, A, 1, B, 1, C, 2, SCALE, INFO )
      CALL CHKXER( 'STRSYL', INFOT, NOUT, LERR, OK )
      INFOT = 9
      CALL STRSYL( 'N', 'N', 1, 0, 2, A, 1, B, 1, C, 1, SCALE, INFO )
      CALL CHKXER( 'STRSYL', INFOT, NOUT, LERR, OK )
      INFOT = 11
      CALL STRSYL( 'N', 'N', 1, 2, 0, A, 2, B, 1, C, 1, SCALE, INFO )
      CALL CHKXER( 'STRSYL', INFOT, NOUT, LERR, OK )
      NT = NT + 8
*
*     Test STREXC
*
      SRNAMT = 'STREXC'
      INFOT = 1
      CALL STREXC( 'X', 1, A, 1, B, 1, 1, 1, WORK, INFO )
      CALL CHKXER( 'STREXC', INFOT, NOUT, LERR, OK )
      INFOT = 7
      CALL STREXC( 'N', 0, A, 1, B, 1, 1, 1, WORK, INFO )
      CALL CHKXER( 'STREXC', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL STREXC( 'N', 2, A, 1, B, 1, 1, 2, WORK, INFO )
      CALL CHKXER( 'STREXC', INFOT, NOUT, LERR, OK )
      INFOT = 6
      CALL STREXC( 'V', 2, A, 2, B, 1, 1, 2, WORK, INFO )
      CALL CHKXER( 'STREXC', INFOT, NOUT, LERR, OK )
      INFOT = 7
      CALL STREXC( 'V', 1, A, 1, B, 1, 0, 1, WORK, INFO )
      CALL CHKXER( 'STREXC', INFOT, NOUT, LERR, OK )
      INFOT = 7
      CALL STREXC( 'V', 1, A, 1, B, 1, 2, 1, WORK, INFO )
      CALL CHKXER( 'STREXC', INFOT, NOUT, LERR, OK )
      INFOT = 8
      CALL STREXC( 'V', 1, A, 1, B, 1, 1, 0, WORK, INFO )
      CALL CHKXER( 'STREXC', INFOT, NOUT, LERR, OK )
      INFOT = 8
      CALL STREXC( 'V', 1, A, 1, B, 1, 1, 2, WORK, INFO )
      CALL CHKXER( 'STREXC', INFOT, NOUT, LERR, OK )
      NT = NT + 8
*
*     Test STRSNA
*
      SRNAMT = 'STRSNA'
      INFOT = 1
      CALL STRSNA( 'X', 'A', SEL, 0, A, 1, B, 1, C, 1, S, SEP, 1, M,
     $             WORK, 1, IWORK, INFO )
      CALL CHKXER( 'STRSNA', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL STRSNA( 'B', 'X', SEL, 0, A, 1, B, 1, C, 1, S, SEP, 1, M,
     $             WORK, 1, IWORK, INFO )
      CALL CHKXER( 'STRSNA', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL STRSNA( 'B', 'A', SEL, -1, A, 1, B, 1, C, 1, S, SEP, 1, M,
     $             WORK, 1, IWORK, INFO )
      CALL CHKXER( 'STRSNA', INFOT, NOUT, LERR, OK )
      INFOT = 6
      CALL STRSNA( 'V', 'A', SEL, 2, A, 1, B, 1, C, 1, S, SEP, 2, M,
     $             WORK, 2, IWORK, INFO )
      CALL CHKXER( 'STRSNA', INFOT, NOUT, LERR, OK )
      INFOT = 8
      CALL STRSNA( 'B', 'A', SEL, 2, A, 2, B, 1, C, 2, S, SEP, 2, M,
     $             WORK, 2, IWORK, INFO )
      CALL CHKXER( 'STRSNA', INFOT, NOUT, LERR, OK )
      INFOT = 10
      CALL STRSNA( 'B', 'A', SEL, 2, A, 2, B, 2, C, 1, S, SEP, 2, M,
     $             WORK, 2, IWORK, INFO )
      CALL CHKXER( 'STRSNA', INFOT, NOUT, LERR, OK )
      INFOT = 13
      CALL STRSNA( 'B', 'A', SEL, 1, A, 1, B, 1, C, 1, S, SEP, 0, M,
     $             WORK, 1, IWORK, INFO )
      CALL CHKXER( 'STRSNA', INFOT, NOUT, LERR, OK )
      INFOT = 13
      CALL STRSNA( 'B', 'S', SEL, 2, A, 2, B, 2, C, 2, S, SEP, 1, M,
     $             WORK, 2, IWORK, INFO )
      CALL CHKXER( 'STRSNA', INFOT, NOUT, LERR, OK )
      INFOT = 16
      CALL STRSNA( 'B', 'A', SEL, 2, A, 2, B, 2, C, 2, S, SEP, 2, M,
     $             WORK, 1, IWORK, INFO )
      CALL CHKXER( 'STRSNA', INFOT, NOUT, LERR, OK )
      NT = NT + 9
*
*     Test STRSEN
*
      SEL( 1 ) = .FALSE.
      SRNAMT = 'STRSEN'
      INFOT = 1
      CALL STRSEN( 'X', 'N', SEL, 0, A, 1, B, 1, WR, WI, M, S( 1 ),
     $             SEP( 1 ), WORK, 1, IWORK, 1, INFO )
      CALL CHKXER( 'STRSEN', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL STRSEN( 'N', 'X', SEL, 0, A, 1, B, 1, WR, WI, M, S( 1 ),
     $             SEP( 1 ), WORK, 1, IWORK, 1, INFO )
      CALL CHKXER( 'STRSEN', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL STRSEN( 'N', 'N', SEL, -1, A, 1, B, 1, WR, WI, M, S( 1 ),
     $             SEP( 1 ), WORK, 1, IWORK, 1, INFO )
      CALL CHKXER( 'STRSEN', INFOT, NOUT, LERR, OK )
      INFOT = 6
      CALL STRSEN( 'N', 'N', SEL, 2, A, 1, B, 1, WR, WI, M, S( 1 ),
     $             SEP( 1 ), WORK, 2, IWORK, 1, INFO )
      CALL CHKXER( 'STRSEN', INFOT, NOUT, LERR, OK )
      INFOT = 8
      CALL STRSEN( 'N', 'V', SEL, 2, A, 2, B, 1, WR, WI, M, S( 1 ),
     $             SEP( 1 ), WORK, 1, IWORK, 1, INFO )
      CALL CHKXER( 'STRSEN', INFOT, NOUT, LERR, OK )
      INFOT = 15
      CALL STRSEN( 'N', 'V', SEL, 2, A, 2, B, 2, WR, WI, M, S( 1 ),
     $             SEP( 1 ), WORK, 0, IWORK, 1, INFO )
      CALL CHKXER( 'STRSEN', INFOT, NOUT, LERR, OK )
      INFOT = 15
      CALL STRSEN( 'E', 'V', SEL, 3, A, 3, B, 3, WR, WI, M, S( 1 ),
     $             SEP( 1 ), WORK, 1, IWORK, 1, INFO )
      CALL CHKXER( 'STRSEN', INFOT, NOUT, LERR, OK )
      INFOT = 15
      CALL STRSEN( 'V', 'V', SEL, 3, A, 3, B, 3, WR, WI, M, S( 1 ),
     $             SEP( 1 ), WORK, 3, IWORK, 2, INFO )
      CALL CHKXER( 'STRSEN', INFOT, NOUT, LERR, OK )
      INFOT = 17
      CALL STRSEN( 'E', 'V', SEL, 2, A, 2, B, 2, WR, WI, M, S( 1 ),
     $             SEP( 1 ), WORK, 1, IWORK, 0, INFO )
      CALL CHKXER( 'STRSEN', INFOT, NOUT, LERR, OK )
      INFOT = 17
      CALL STRSEN( 'V', 'V', SEL, 3, A, 3, B, 3, WR, WI, M, S( 1 ),
     $             SEP( 1 ), WORK, 4, IWORK, 1, INFO )
      CALL CHKXER( 'STRSEN', INFOT, NOUT, LERR, OK )
      NT = NT + 10
*
*     Print a summary line.
*
      IF( OK ) THEN
         WRITE( NOUT, FMT = 9999 )PATH, NT
      ELSE
         WRITE( NOUT, FMT = 9998 )PATH
      END IF
*
      RETURN
 9999 FORMAT( 1X, A3, ' routines passed the tests of the error exits (',
     $      I3, ' tests done)' )
 9998 FORMAT( ' *** ', A3, ' routines failed the tests of the error ex',
     $      'its ***' )
*
*     End of SERREC
*
      END