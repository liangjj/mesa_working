      SUBROUTINE CTIM22( LINE, NSIZES, NN, NTYPES, DOTYPE, NPARMS, NNB,
     $                   LDAS, TIMMIN, NOUT, ISEED, A, D, E, E2, U, URE,
     $                   UIM, TAU, TAURE, Z, ZRE, ZIM, WORK, LWORK,
     $                   RWORK, LLWORK, IWORK, TIMES, LDT1, LDT2, LDT3,
     $                   OPCNTS, LDO1, LDO2, LDO3, INFO )
*
*  -- LAPACK timing routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993 
*
*     .. Scalar Arguments ..
      CHARACTER*80       LINE
      INTEGER            INFO, LDO1, LDO2, LDO3, LDT1, LDT2, LDT3,
     $                   LWORK, NOUT, NPARMS, NSIZES, NTYPES
      REAL               TIMMIN
*     ..
*     .. Array Arguments ..
      LOGICAL            DOTYPE( * ), LLWORK( * )
      INTEGER            ISEED( * ), IWORK( * ), LDAS( * ), NN( * ),
     $                   NNB( * )
      REAL               D( * ), E( * ), E2( * ),
     $                   OPCNTS( LDO1, LDO2, LDO3, * ), RWORK( * ),
     $                   TAURE( * ), TIMES( LDT1, LDT2, LDT3, * ),
     $                   UIM( * ), URE( * ), ZIM( * ), ZRE( * )
      COMPLEX            A( * ), TAU( * ), U( * ), WORK( * ), Z( * )
*     ..
*
*  Purpose
*  =======
*
*     CTIM22 times the LAPACK routines for the complex hermitian
*     eigenvalue problem.
*
*     For each N value in NN(1:NSIZES) and .TRUE. value in
*     DOTYPE(1:NTYPES), a matrix will be generated and used to test the
*     selected routines.  Thus, NSIZES*(number of .TRUE. values in
*     DOTYPE) matrices will be generated.
*
*  Arguments
*  =========
*
*  LINE    (input) CHARACTER*80
*          On entry, LINE contains the input line which requested
*          this routine.  This line may contain a subroutine name,
*          such as CHETRD, indicating that only routine CHETRD will
*          be timed, or it may contain a generic name, such as SST.
*          In this case, the rest of the line is scanned for the
*          first 11 non-blank characters, corresponding to the eleven
*          combinations of subroutine and options:
*          LAPACK:
*            1: CHETRD
*            2: CSTEQR(VECT='N')
*            3: CUNGTR+CSTEQR(VECT='V') (compare with IMTQL2+HTRIBK)
*            4: CPTEQR(VECT='N')
*            5: CUNGTR+CPTEQR(VECT='V')
*            6. SSTEBZ+CSTEIN+CUNMTR
*          EISPACK:
*            7: HTRIDI (compare with CHETRD)
*            8: IMTQL1 (compare w/ CSTEQR -- VECT='N')
*            9: IMTQL2+HTRIBK (compare w/ CUNGTR+CSTEQR(VECT='V') )
*          If a character is 'T' or 't', the corresponding routine in
*          this path is timed.  If the entire line is blank, all the
*          routines in the path are timed.
*
*  NSIZES  (input) INTEGER
*          The number of values of N contained in the vector NN.
*
*  NN      (input) INTEGER array, dimension( NSIZES )
*          The values of the matrix size N to be tested.  For each
*          N value in the array NN, and each .TRUE. value in DOTYPE,
*          a matrix A will be generated and used to test the routines.
*
*  NTYPES  (input) INTEGER
*          The number of types in DOTYPE.  Only the first MAXTYP
*          elements will be examined.  Exception: if NSIZES=1 and
*          NTYPES=MAXTYP+1, and DOTYPE=MAXTYP*f,t, then the input
*          value of A will be used.
*
*  DOTYPE  (input) LOGICAL
*          If DOTYPE(j) is .TRUE., then a matrix of type j will be
*          generated.  The matrix A has the form X**(-1) D X, where
*          X is unitary and D is diagonal with:
*          (j=1)  evenly spaced entries 1, ..., ULP with random signs.
*          (j=2)  geometrically spaced entries 1, ..., ULP with random
*                 signs.
*          (j=3)  "clustered" entries 1, ULP,..., ULP with random
*                 signs.
*          (j=4)  entries randomly chosen from ( ULP, 1 ).
*
*  NPARMS  (input) INTEGER
*          The number of values in each of the arrays NNB and LDAS.
*          For each matrix A generated according to NN and DOTYPE,
*          tests will be run with (NB,LDA)=
*          (NNB(1),LDAS(1)),...,(NNB(NPARMS), LDAS(NPARMS))
*
*  NNB     (input) INTEGER array, dimension( NPARMS )
*          The values of the blocksize ("NB") to be tested.
*
*  LDAS    (input) INTEGER array, dimension( NPARMS )
*          The values of LDA, the leading dimension of all matrices,
*          to be tested.
*
*  TIMMIN  (input) REAL
*          The minimum time a subroutine will be timed.
*
*  NOUT    (input) INTEGER
*          If NOUT > 0 then NOUT specifies the unit number
*          on which the output will be printed.  If NOUT <= 0, no
*          output is printed.
*
*  ISEED   (input/output) INTEGER array, dimension( 4 )
*          The random seed used by the random number generator, used
*          by the test matrix generator.  It is used and updated on
*          each call to CTIM22
*
*  A       (workspace) COMPLEX array, dimension( max(NN)*max(LDAS) )
*          The original matrix to be tested.
*
*  D       (workspace) REAL array, dimension( max(NN) )
*          The diagonal of the tridiagonal generated by CHETRD/HTRIDI.
*
*  E       (workspace) REAL array, dimension( max(NN) )
*          The off-diagonal of the tridiagonal generated by
*          CHETRD/HTRIDI.
*
*  E2      (workspace) REAL array, dimension( max(NN) )
*          The diagonal of a positive definite tridiagonal matrix
*          sent to CPTEQR.  The off-diagonal is in array E.
*
*  U       (workspace) COMPLEX array, dimension( max(NN)*max(LDAS) )
*          The array of Householder vectors output by CHETRD.  This
*          array is used only when URE and UIM are not; thus, on
*          nearly all computers, URE may be EQUIVALENCEd with the
*          first half of U in the main (calling) routine, and UIM with
*          the second half, although this is a violation of the
*          FORTRAN-77 standard.
*
*  URE     (workspace) REAL array,
*                      dimension( max(NN)*max(LDAS) )
*          The array of the real parts of Householder vectors output by
*          HTRIDI.  This array is used only when U is not -- see the
*          note description of U.
*
*  UIM     (workspace) REAL array,
*                      dimension( max(NN)*max(LDAS) )
*          The array of the imaginary parts of Householder vectors
*          output by HTRIDI.  This array is used only when U is not --
*          see the description of U.
*
*  TAU     (workspace) COMPLEX array, dimension( max(NN) )
*          The vector of coefficients for the Householder
*          transformations output by CHETRD.  This array is used only
*          when TAURE is not; thus, on nearly all computers, TAURE may
*          be EQUIVALENCEd with TAU in the main (calling) routine,
*          although this is a violation of the FORTRAN-77 standard.
*
*  TAURE   (workspace) REAL array, dimension( 2*max(NN) )
*          The vector of complex (modulus 1) factors output by HTRIDI.
*          This vector is used only when TAU is not -- see the
*          description of TAU.
*
*  Z       (workspace) COMPLEX array, dimension( max(NN)*max(LDAS) )
*          Various output arrays.  This array is used only when ZRE
*          and ZIM are not; thus, on nearly all computers, ZRE may be
*          EQUIVALENCEd with the first half of Z in the main (calling)
*          routine, and ZIM with the second half, although this is a
*          violation of the FORTRAN-77 standard.
*
*  ZRE     (workspace) REAL array,
*                      dimension( max(NN)*max(LDAS) )
*          Various output arrays (real parts).  This array is used
*          only when Z is not -- see the description of Z.
*
*  ZIM     (workspace) REAL array,
*                      dimension( max(NN)*max(LDAS) )
*          Various output arrays (imaginary parts).  This array is
*          used only when Z is not -- see the description of Z.
*
*  WORK    (workspace) COMPLEX array, dimension( LWORK )
*
*  LWORK   (input) INTEGER
*          Number of elements in WORK.  It must be at least
*          max( (NNB + 2 )*LDAS )
*
*  RWORK   (workspace) REAL array, dimension
*                   ( max( 4*max(LDAS), NSIZES*NTYPES*NPARMS ) )
*          This should *not* be equivalenced to other arrays.
*
*  LLWORK  (workspace) LOGICAL array, dimension( NPARMS )
*
*  IWORK   (workspace) INTEGER array, dimension ( 5*max(LDAS) )
*
*  TIMES   (workspace) REAL array,
*                      dimension (LDT1,LDT2,LDT3,NSUBS)
*          TIMES(i,j,k,l) will be set to the run time (in seconds) for
*          subroutine l, with N=NN(k), matrix type j, and LDA=LDAS(i),
*          NBLOCK=NNB(i).
*
*  LDT1    (input) INTEGER
*          The first dimension of TIMES.  LDT1 >= min( 1, NPARMS ).
*
*  LDT2    (input) INTEGER
*          The second dimension of TIMES.  LDT2 >= min( 1, NTYPES ).
*
*  LDT3    (input) INTEGER
*          The third dimension of TIMES.  LDT3 >= min( 1, NSIZES ).
*
*  OPCNTS  (output) REAL array,
*                   dimension (LDO1,LDO2,LDO3,NSUBS)
*          OPCNTS(i,j,k,l) will be set to the number of floating-point
*          operations executed by subroutine l, with N=NN(k), matrix
*          type j, and LDA=LDAS(i), NBLOCK=NNB(i).
*
*  LDO1    (input) INTEGER
*          The first dimension of OPCNTS.  LDO1 >= min( 1, NPARMS ).
*
*  LDO2    (input) INTEGER
*          The second dimension of OPCNTS.  LDO2 >= min( 1, NTYPES ).
*
*  LDO3    (input) INTEGER
*          The third dimension of OPCNTS.  LDO3 >= min( 1, NSIZES ).
*
*  INFO    (output) INTEGER
*          Error flag.  It will be set to zero if no error occurred.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            MAXTYP, NSUBS
      PARAMETER          ( MAXTYP = 4, NSUBS = 9 )
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            RUNHTR, RUNTRD
      CHARACTER*1        UPLO
      INTEGER            I, IC, IINFO, ILWORK, IMODE, IN, IPAR, ISUB,
     $                   ITYPE, J, J1, J2, J3, J4, LASTL, LDA, LDU,
     $                   MTYPES, N, NB, IL, IU, M, NSPLIT
      REAL               S1, S2, TIME, ULP, ULPINV, UNTIME, VL, VU,
     $                   ABSTOL
*     ..
*     .. Local Arrays ..
      LOGICAL            TIMSUB( NSUBS )
      CHARACTER*4        PNAMES( 4 )
      CHARACTER*20       SUBNAM( NSUBS )
      INTEGER            IDUMMA( 1 ), INPARM( NSUBS ), IOLDSD( 4 ),
     $                   KMODE( MAXTYP )
*     ..
*     .. External Functions ..
      REAL               SECOND, SLAMCH, SOPLA
      EXTERNAL           SECOND, SLAMCH, SOPLA
*     ..
*     .. External Subroutines ..
      EXTERNAL           ATIMIN, CHETRD, CLACPY, CLATMS, CPTEQR, CSTEIN,
     $                   CSTEQR, CUNGTR, CUNMTR, HTRIBK, HTRIDI, IMTQL1,
     $                   IMTQL2, SCOPY, SLAZRO, SPRTBE, SSTEBZ, XLAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, AIMAG, MAX, MIN, REAL
*     ..
*     .. Common blocks ..
      COMMON             / LATIME / OPS, ITCNT
*     ..
*     .. Scalars in Common ..
      REAL               ITCNT, OPS
*     ..
*     .. Data statements ..
      DATA               SUBNAM / 'CHETRD', 'CSTEQR(N)',
     $                   'CUNGTR+CSTEQR(V)', 'CPTEQR(N)',
     $                   'CUNGTR+CPTEQR(V)', 'SSTEBZ+CSTEIN+CUNMTR',
     $                   'HTRIDI', 'IMTQL1', 'IMTQL2+HTRIBK' /
      DATA               INPARM / 2, 1, 2, 1, 2, 2, 1, 1, 1 /
      DATA               PNAMES / 'LDA', 'NB', 'bad1', 'bad2' /
      DATA               KMODE / 4, 3, 1, 5 /
*     ..
*     .. Executable Statements ..
*
*
*     Extract the timing request from the input line.
*
      CALL ATIMIN( 'CST', LINE, NSUBS, SUBNAM, TIMSUB, NOUT, INFO )
      IF( INFO.NE.0 )
     $   RETURN
*
*     Check that N <= LDA for the input values.
*
      DO 20 J2 = 1, NSIZES
         DO 10 J1 = 1, NPARMS
            IF( NN( J2 ).GT.LDAS( J1 ) ) THEN
               INFO = -8
               WRITE( NOUT, FMT = 9999 )LINE( 1: 6 )
 9999          FORMAT( 1X, A, ' timing run not attempted -- N > LDA',
     $               / )
               RETURN
            END IF
   10    CONTINUE
   20 CONTINUE
*
*     Check LWORK
*
      ILWORK = 0
      DO 30 J1 = 1, NPARMS
         ILWORK = MAX( ILWORK, ( NNB( J1 )+2 )*LDAS( J1 ) )
   30 CONTINUE
      IF( ILWORK.GT.LWORK ) THEN
         INFO = -18
         WRITE( NOUT, FMT = 9998 )LINE( 1: 6 )
 9998    FORMAT( 1X, A, ' timing run not attempted -- LWORK too small.',
     $         / )
         RETURN
      END IF
*
*     Check to see whether CHETRD must be run.
*
*     RUNTRD -- if CHETRD must be run.
*
      RUNTRD = .FALSE.
      IF( TIMSUB( 2 ) .OR. TIMSUB( 3 ) .OR. TIMSUB( 4 ) .OR.
     $    TIMSUB( 5 ) )RUNTRD = .TRUE.
*
*     Check to see whether HTRIDI must be run.
*
*     RUNHTR -- if HTRIDI must be run.
*
      RUNHTR = .FALSE.
      IF( TIMSUB( 7 ) .OR. TIMSUB( 8 ) )
     $   RUNHTR = .TRUE.
*
*     Various Constants
*
      ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' )
      ULPINV = ONE / ULP
*
*     Zero out OPCNTS, TIMES
*
      DO 70 J4 = 1, NSUBS
         DO 60 J3 = 1, NSIZES
            DO 50 J2 = 1, NTYPES
               DO 40 J1 = 1, NPARMS
                  OPCNTS( J1, J2, J3, J4 ) = ZERO
                  TIMES( J1, J2, J3, J4 ) = ZERO
   40          CONTINUE
   50       CONTINUE
   60    CONTINUE
   70 CONTINUE
*
*     Do for each value of N:
*
      DO 540 IN = 1, NSIZES
*
         N = NN( IN )
*
*        Do for each .TRUE. value in DOTYPE:
*
         MTYPES = MIN( MAXTYP, NTYPES )
         IF( NTYPES.EQ.MAXTYP+1 .AND. NSIZES.EQ.1 )
     $      MTYPES = NTYPES
         DO 530 ITYPE = 1, MTYPES
            IF( .NOT.DOTYPE( ITYPE ) )
     $         GO TO 530
*
*           Save random number seed for error messages
*
            DO 80 J = 1, 4
               IOLDSD( J ) = ISEED( J )
   80       CONTINUE
*
*-----------------------------------------------------------------------
*
*           Time the LAPACK Routines
*
*           Generate A
*
            UPLO = 'L'
            IF( ITYPE.LE.MAXTYP ) THEN
               IMODE = KMODE( ITYPE )
               CALL CLATMS( N, N, 'S', ISEED, 'S', RWORK, IMODE, ULPINV,
     $                      ONE, N, N, UPLO, A, N, WORK, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUT, FMT = 9997 )'CLATMS', IINFO, N, ITYPE,
     $               0, IOLDSD
                  INFO = ABS( IINFO )
                  GO TO 530
               END IF
            END IF
*
*           Time CHETRD for each pair NNB(j), LDAS(j)
*
            IF( TIMSUB( 1 ) ) THEN
               DO 110 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
                  NB = MIN( N, NNB( IPAR ) )
                  CALL XLAENV( 1, NB )
                  CALL XLAENV( 2, 2 )
                  CALL XLAENV( 3, NB )
*
*                 Time CHETRD
*
                  IC = 0
                  OPS = ZERO
                  S1 = SECOND( )
   90             CONTINUE
                  CALL CLACPY( UPLO, N, N, A, N, U, LDA )
                  CALL CHETRD( UPLO, N, U, LDA, D, E, TAU, WORK, LWORK,
     $                         IINFO )
                  IF( IINFO.NE.0 ) THEN
                     WRITE( NOUT, FMT = 9997 )SUBNAM( 1 ), IINFO, N,
     $                  ITYPE, IPAR, IOLDSD
                     INFO = ABS( IINFO )
                     GO TO 190
                  END IF
*
                  S2 = SECOND( )
                  TIME = S2 - S1
                  IC = IC + 1
                  IF( TIME.LT.TIMMIN )
     $               GO TO 90
*
*                 Subtract the time used in CLACPY.
*
                  S1 = SECOND( )
                  DO 100 J = 1, IC
                     CALL CLACPY( UPLO, N, N, A, N, Z, LDA )
  100             CONTINUE
                  S2 = SECOND( )
                  UNTIME = S2 - S1
*
                  TIMES( IPAR, ITYPE, IN, 1 ) = MAX( TIME-UNTIME,
     $               ZERO ) / REAL( IC )
                  OPCNTS( IPAR, ITYPE, IN, 1 ) = SOPLA( 'CHETRD', N,
     $               0, 0, 0, NB )
                  LDU = LDA
  110          CONTINUE
            ELSE
               IF( RUNTRD ) THEN
                  CALL CLACPY( UPLO, N, N, A, N, U, N )
                  CALL CHETRD( UPLO, N, U, N, D, E, TAU, WORK, LWORK,
     $                         IINFO )
                  IF( IINFO.NE.0 ) THEN
                     WRITE( NOUT, FMT = 9997 )SUBNAM( 1 ), IINFO, N,
     $                  ITYPE, 0, IOLDSD
                     INFO = ABS( IINFO )
                     GO TO 190
                  END IF
                  LDU = N
               END IF
            END IF
*
*           Time CSTEQR for each distinct LDA=LDAS(j)
*
            IF( TIMSUB( 2 ) ) THEN
               DO 150 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
*
*                 If this value of LDA has come up before, just use
*                 the value previously computed.
*
                  LASTL = 0
                  DO 120 J = 1, IPAR - 1
                     IF( LDA.EQ.LDAS( J ) )
     $                  LASTL = J
  120             CONTINUE
                  IF( LASTL.EQ.0 ) THEN
*
*                    Time CSTEQR with VECT='N'
*
                     IC = 0
                     OPS = ZERO
                     S1 = SECOND( )
  130                CONTINUE
                     CALL SCOPY( N, D, 1, RWORK, 1 )
                     CALL SCOPY( N-1, E, 1, RWORK( LDA+1 ), 1 )
                     CALL CSTEQR( 'N', N, RWORK, RWORK( LDA+1 ), Z, LDA,
     $                            RWORK( 2*LDA+1 ), IINFO )
                     IF( IINFO.NE.0 ) THEN
                        WRITE( NOUT, FMT = 9997 )SUBNAM( 2 ), IINFO, N,
     $                     ITYPE, IPAR, IOLDSD
                        INFO = ABS( IINFO )
                        GO TO 150
                     END IF
                     S2 = SECOND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN )
     $                  GO TO 130
*
*                    Subtract the time used in SCOPY.
*
                     S1 = SECOND( )
                     DO 140 J = 1, IC
                        CALL SCOPY( N, D, 1, RWORK, 1 )
                        CALL SCOPY( N-1, E, 1, RWORK( LDA+1 ), 1 )
  140                CONTINUE
                     S2 = SECOND( )
                     UNTIME = S2 - S1
*
                     TIMES( IPAR, ITYPE, IN, 2 ) = MAX( TIME-UNTIME,
     $                  ZERO ) / REAL( IC )
                     OPCNTS( IPAR, ITYPE, IN, 2 ) = OPS / REAL( IC )
                  ELSE
                     OPCNTS( IPAR, ITYPE, IN, 2 ) = OPCNTS( LASTL,
     $                  ITYPE, IN, 2 )
                     TIMES( IPAR, ITYPE, IN, 2 ) = TIMES( LASTL, ITYPE,
     $                  IN, 2 )
                  END IF
  150          CONTINUE
            END IF
*
*           Time CUNGTR + CSTEQR(VECT='V') for each pair NNB(j), LDAS(j)
*
            IF( TIMSUB( 3 ) ) THEN
               DO 180 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
                  NB = MIN( N, NNB( IPAR ) )
                  CALL XLAENV( 1, NB )
                  CALL XLAENV( 2, 2 )
                  CALL XLAENV( 3, NB )
*
*                 Time CUNGTR + CSTEQR
*
                  IC = 0
                  OPS = ZERO
                  S1 = SECOND( )
  160             CONTINUE
                  CALL CLACPY( 'L', N, N, A, N, Z, LDA )
                  CALL CUNGTR( 'L', N, Z, LDA, TAU, WORK, LWORK, IINFO )
                  IF( IINFO.NE.0 ) THEN
                     WRITE( NOUT, FMT = 9997 )'CUNGTR', IINFO, N,
     $                  ITYPE, IPAR, IOLDSD
                     INFO = ABS( IINFO )
                     GO TO 180
                  END IF
                  CALL SCOPY( N, D, 1, RWORK, 1 )
                  CALL SCOPY( N-1, E, 1, RWORK( LDA+1 ), 1 )
                  CALL CSTEQR( 'V', N, RWORK, RWORK( LDA+1 ), Z, LDA,
     $                         RWORK( 2*LDA+1 ), IINFO )
                  IF( IINFO.NE.0 ) THEN
                     WRITE( NOUT, FMT = 9997 )SUBNAM( 3 ), IINFO, N,
     $                  ITYPE, IPAR, IOLDSD
                     INFO = ABS( IINFO )
                     GO TO 180
                  END IF
*
                  S2 = SECOND( )
                  TIME = S2 - S1
                  IC = IC + 1
                  IF( TIME.LT.TIMMIN )
     $               GO TO 160
*
*                 Subtract the time used in CLACPY.
*
                  S1 = SECOND( )
                  DO 170 J = 1, IC
                     CALL SCOPY( N, D, 1, RWORK, 1 )
                     CALL SCOPY( N-1, E, 1, RWORK( LDA+1 ), 1 )
                     CALL CLACPY( 'L', N, N, A, N, Z, LDA )
  170             CONTINUE
                  S2 = SECOND( )
                  UNTIME = S2 - S1
*
                  TIMES( IPAR, ITYPE, IN, 3 ) = MAX( TIME-UNTIME,
     $               ZERO ) / REAL( IC )
                  OPCNTS( IPAR, ITYPE, IN, 3 ) = OPS / REAL( IC )
                  LDU = LDA
  180          CONTINUE
            END IF
*
  190       CONTINUE
*
*           Time CPTEQR for each distinct LDA=LDAS(j)
*
            IF( TIMSUB( 4 ) ) THEN
               DO 240 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
*
*                 If this value of LDA has come up before, just use
*                 the value previously computed.
*
                  LASTL = 0
                  DO 200 J = 1, IPAR - 1
                     IF( LDA.EQ.LDAS( J ) )
     $                  LASTL = J
  200             CONTINUE
                  IF( LASTL.EQ.0 ) THEN
*
*                    Time CPTEQR with VECT='N'
*
*
*                    Modify the tridiagonal matrix to make it
*                    positive definite.
                     E2( 1 ) = ABS( D( 1 ) ) + ABS( E( 1 ) )
                     DO 210 I = 2, N - 1
                        E2( I ) = ABS( D( I ) ) + ABS( E( I ) ) +
     $                            ABS( E( I-1 ) )
  210                CONTINUE
                     E2( N ) = ABS( D( N ) ) + ABS( E( N-1 ) )
                     IC = 0
                     OPS = ZERO
                     S1 = SECOND( )
  220                CONTINUE
                     CALL SCOPY( N, E2, 1, RWORK( 1 ), 1 )
                     CALL SCOPY( N-1, E, 1, RWORK( LDA+1 ), 1 )
                     CALL CPTEQR( 'N', N, RWORK, RWORK( LDA+1 ), Z, LDA,
     $                            RWORK( 2*LDA+1 ), IINFO )
                     IF( IINFO.NE.0 ) THEN
                        WRITE( NOUT, FMT = 9997 )SUBNAM( 4 ), IINFO, N,
     $                     ITYPE, IPAR, IOLDSD
                        INFO = ABS( IINFO )
                        GO TO 240
                     END IF
                     S2 = SECOND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN )
     $                  GO TO 220
*
*                    Subtract the time used in SCOPY.
*
                     S1 = SECOND( )
                     DO 230 J = 1, IC
                        CALL SCOPY( N, E2, 1, RWORK, 1 )
                        CALL SCOPY( N-1, E, 1, RWORK( LDA+1 ), 1 )
  230                CONTINUE
                     S2 = SECOND( )
                     UNTIME = S2 - S1
*
                     TIMES( IPAR, ITYPE, IN, 4 ) = MAX( TIME-UNTIME,
     $                  ZERO ) / REAL( IC )
                     OPCNTS( IPAR, ITYPE, IN, 4 ) = OPS / REAL( IC )
                  ELSE
                     OPCNTS( IPAR, ITYPE, IN, 4 ) = OPCNTS( LASTL,
     $                  ITYPE, IN, 4 )
                     TIMES( IPAR, ITYPE, IN, 4 ) = TIMES( LASTL, ITYPE,
     $                  IN, 4 )
                  END IF
  240          CONTINUE
            END IF
*
*           Time CUNGTR + CPTEQR(VECT='V') for each pair NNB(j), LDAS(j)
*
            IF( TIMSUB( 5 ) ) THEN
               DO 290 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
                  NB = MIN( N, NNB( IPAR ) )
                  CALL XLAENV( 1, NB )
                  CALL XLAENV( 2, 2 )
                  CALL XLAENV( 3, NB )
*
*                 Time CUNGTR + CPTEQR
*
                  IC = 0
                  OPS = ZERO
                  S1 = SECOND( )
  250             CONTINUE
                  CALL CLACPY( 'L', N, N, A, N, Z, LDA )
                  CALL CUNGTR( 'L', N, Z, LDA, TAU, WORK, LWORK, IINFO )
                  IF( IINFO.NE.0 ) THEN
                     WRITE( NOUT, FMT = 9997 )'CUNGTR', IINFO, N,
     $                  ITYPE, IPAR, IOLDSD
                     INFO = ABS( IINFO )
                     GO TO 290
                  END IF
*
*                 Modify the tridiagonal matrix to make it
*                 positive definite.
                  E2( 1 ) = ABS( D( 1 ) ) + ABS( E( 1 ) )
                  DO 260 I = 2, N - 1
                     E2( I ) = ABS( D( I ) ) + ABS( E( I ) ) +
     $                         ABS( E( I-1 ) )
  260             CONTINUE
                  E2( N ) = ABS( D( N ) ) + ABS( E( N-1 ) )
*
                  CALL SCOPY( N, E2, 1, RWORK, 1 )
                  CALL SCOPY( N-1, E, 1, RWORK( LDA+1 ), 1 )
                  CALL CPTEQR( 'V', N, RWORK, RWORK( LDA+1 ), Z, LDA,
     $                         RWORK( 2*LDA+1 ), IINFO )
                  IF( IINFO.NE.0 ) THEN
                     WRITE( NOUT, FMT = 9997 )SUBNAM( 5 ), IINFO, N,
     $                  ITYPE, IPAR, IOLDSD
                     INFO = ABS( IINFO )
                     GO TO 290
                  END IF
*
                  S2 = SECOND( )
                  TIME = S2 - S1
                  IC = IC + 1
                  IF( TIME.LT.TIMMIN )
     $               GO TO 250
*
*                 Subtract the time used in CLACPY.
*
                  S1 = SECOND( )
                  DO 280 J = 1, IC
                     E2( 1 ) = ABS( D( 1 ) ) + ABS( E( 1 ) )
                     DO 270 I = 2, N - 1
                        E2( I ) = ABS( D( I ) ) + ABS( E( I ) ) +
     $                            ABS( E( I-1 ) )
  270                CONTINUE
                     E2( N ) = ABS( D( N ) ) + ABS( E( N-1 ) )
*
                     CALL SCOPY( N, E2, 1, RWORK, 1 )
                     CALL SCOPY( N-1, E, 1, RWORK( LDA+1 ), 1 )
                     CALL CLACPY( 'L', N, N, A, N, Z, LDA )
  280             CONTINUE
                  S2 = SECOND( )
                  UNTIME = S2 - S1
*
                  TIMES( IPAR, ITYPE, IN, 5 ) = MAX( TIME-UNTIME,
     $               ZERO ) / REAL( IC )
                  OPCNTS( IPAR, ITYPE, IN, 5 ) = OPS / REAL( IC )
                  LDU = LDA
  290          CONTINUE
            END IF
*
*           Time SSTEBZ+CSTEIN+CUNMTR for each pair NNB(j), LDAS(j)
*
            IF( TIMSUB( 6 ) ) THEN
               VL = ZERO
               VU = ZERO
               IL = 1
               IU = N
               ABSTOL = ZERO
               DO 181 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
                  NB = MIN( N, NNB( IPAR ) )
                  CALL XLAENV( 1, NB )
                  CALL XLAENV( 2, 2 )
*
*                 Time SSTEBZ + CSTEIN + CUNMTR
*
                  IC = 0
                  OPS = ZERO
                  S1 = SECOND( )
  161             CONTINUE
*
                  CALL SSTEBZ( 'A', 'B', N, VL, VU, IL, IU, ABSTOL, D,
     $                         E, M, NSPLIT, RWORK( 1 ), IWORK( 1 ),
     $                         IWORK( N+1 ), RWORK( 2*N+1 ),
     $                         IWORK( 2*N+1 ), IINFO )
                  IF( IINFO.NE.0 ) THEN
                     WRITE( NOUT, FMT = 9997 )'SSTEBZ', IINFO, N,
     $                  ITYPE, IPAR, IOLDSD
                     INFO = ABS( IINFO )
                     GO TO 181
                  END IF
*
                  CALL CSTEIN( N, D, E, N, RWORK( 1 ), IWORK( 1 ),
     $                         IWORK( N+1 ), Z, LDA, RWORK( N+1 ),
     $                         IWORK( 2*N+1 ), IWORK( 3*N+1 ), IINFO )
                  IF( IINFO.NE.0 ) THEN
                     WRITE( NOUT, FMT = 9997 )SUBNAM( 6 ), IINFO, N,
     $                  ITYPE, IPAR, IOLDSD
                     INFO = ABS( IINFO )
                     GO TO 181
                  END IF
*
                  CALL CUNMTR( 'L', 'L', 'N', N, N, U, LDU, TAU, Z,
     $                         LDA, WORK, LWORK, IINFO )
                  IF( IINFO.NE.0 ) THEN
                     WRITE( NOUT, FMT = 9997 )'CUNMTR', IINFO, N,
     $                  ITYPE, IPAR, IOLDSD
                     INFO = ABS( IINFO )
                     GO TO 181
                  END IF
*
                  S2 = SECOND( )
                  TIME = S2 - S1
                  IC = IC + 1
                  IF( TIME.LT.TIMMIN )
     $               GO TO 161
                  UNTIME = ZERO
*
                  TIMES( IPAR, ITYPE, IN, 6 ) = MAX( TIME-UNTIME,
     $               ZERO ) / REAL( IC )
                  OPCNTS( IPAR, ITYPE, IN, 6 ) = OPS / REAL( IC )
                  LDU = LDA
  181          CONTINUE
            END IF
*
*-----------------------------------------------------------------------
*
*           Time the EISPACK Routines
*
*           Skip routines if N <= 0 (EISPACK requirement)
*
            IF( N.LE.0 )
     $         GO TO 530
*
*           Time HTRIDI for each LDAS(j)
*
            IF( TIMSUB( 7 ) ) THEN
               DO 370 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
*
*                 If this value of LDA has come up before, just use
*                 the value previously computed.
*
                  LASTL = 0
                  DO 300 J = 1, IPAR - 1
                     IF( LDA.EQ.LDAS( J ) )
     $                  LASTL = J
  300             CONTINUE
*
                  IF( LASTL.EQ.0 ) THEN
*
*                    Time HTRIDI
*
                     IC = 0
                     OPS = ZERO
                     S1 = SECOND( )
  310                CONTINUE
                     DO 330 J2 = 0, N - 1
                        DO 320 J1 = 1, N
                           URE( J1+LDA*J2 ) = REAL( A( J1+N*J2 ) )
                           UIM( J1+LDA*J2 ) = AIMAG( A( J1+N*J2 ) )
  320                   CONTINUE
  330                CONTINUE
                     CALL HTRIDI( LDA, N, URE, UIM, D, E, RWORK, TAURE )
                     S2 = SECOND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN )
     $                  GO TO 310
*
*                    Subtract the time used in copying A.
*
                     S1 = SECOND( )
                     DO 360 J = 1, IC
                        DO 350 J2 = 0, N - 1
                           DO 340 J1 = 1, N
                              ZRE( J1+LDA*J2 ) = REAL( A( J1+N*J2 ) )
                              ZIM( J1+LDA*J2 ) = AIMAG( A( J1+N*J2 ) )
  340                      CONTINUE
  350                   CONTINUE
  360                CONTINUE
                     S2 = SECOND( )
                     UNTIME = S2 - S1
*
                     TIMES( IPAR, ITYPE, IN, 7 ) = MAX( TIME-UNTIME,
     $                  ZERO ) / REAL( IC )
                     OPCNTS( IPAR, ITYPE, IN, 7 ) = OPS / REAL( IC )
                     LDU = LDA
                  ELSE
                     OPCNTS( IPAR, ITYPE, IN, 7 ) = OPCNTS( LASTL,
     $                  ITYPE, IN, 7 )
                     TIMES( IPAR, ITYPE, IN, 7 ) = TIMES( LASTL, ITYPE,
     $                  IN, 7 )
                  END IF
  370          CONTINUE
            ELSE
               IF( RUNHTR ) THEN
                  DO 390 J2 = 0, N - 1
                     DO 380 J1 = 1, N
                        URE( J1+N*J2 ) = REAL( A( J1+N*J2 ) )
                        UIM( J1+N*J2 ) = AIMAG( A( J1+N*J2 ) )
  380                CONTINUE
  390             CONTINUE
                  CALL HTRIDI( N, N, URE, UIM, D, E, RWORK, TAURE )
                  LDU = N
               END IF
            END IF
*
*           Time IMTQL1 for each LDAS(j)
*
            IF( TIMSUB( 8 ) ) THEN
               DO 430 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
*
*                 If this value of LDA has come up before, just use
*                 the value previously computed.
*
                  LASTL = 0
                  DO 400 J = 1, IPAR - 1
                     IF( LDA.EQ.LDAS( J ) )
     $                  LASTL = J
  400             CONTINUE
*
                  IF( LASTL.EQ.0 ) THEN
*
*                    Time IMTQL1
*
                     IC = 0
                     OPS = ZERO
                     S1 = SECOND( )
  410                CONTINUE
                     CALL SCOPY( N, D, 1, RWORK, 1 )
                     CALL SCOPY( N-1, E, 1, RWORK( LDA+1 ), 1 )
                     CALL IMTQL1( N, RWORK, RWORK( LDA+1 ), IINFO )
                     IF( IINFO.NE.0 ) THEN
                        WRITE( NOUT, FMT = 9997 )SUBNAM( 8 ), IINFO, N,
     $                     ITYPE, IPAR, IOLDSD
                        INFO = ABS( IINFO )
                        GO TO 440
                     END IF
                     S2 = SECOND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN )
     $                  GO TO 410
*
*                    Subtract the time used in SCOPY
*
                     S1 = SECOND( )
                     DO 420 J = 1, IC
                        CALL SCOPY( N, D, 1, RWORK, 1 )
                        CALL SCOPY( N-1, E, 1, RWORK( LDA+1 ), 1 )
  420                CONTINUE
                     S2 = SECOND( )
                     UNTIME = S2 - S1
*
                     TIMES( IPAR, ITYPE, IN, 8 ) = MAX( TIME-UNTIME,
     $                  ZERO ) / REAL( IC )
                     OPCNTS( IPAR, ITYPE, IN, 8 ) = OPS / REAL( IC )
                  ELSE
                     OPCNTS( IPAR, ITYPE, IN, 8 ) = OPCNTS( LASTL,
     $                  ITYPE, IN, 8 )
                     TIMES( IPAR, ITYPE, IN, 8 ) = TIMES( LASTL, ITYPE,
     $                  IN, 8 )
                  END IF
  430          CONTINUE
            END IF
  440       CONTINUE
*
*           Time IMTQL2 + HTRIBK for each LDAS(j)
*
            IF( TIMSUB( 9 ) ) THEN
               DO 520 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
*
*                 If this value of LDA has come up before, just use
*                 the value previously computed.
*
                  LASTL = 0
                  DO 450 J = 1, IPAR - 1
                     IF( LDA.EQ.LDAS( J ) )
     $                  LASTL = J
  450             CONTINUE
*
                  IF( LASTL.EQ.0 ) THEN
*
*                    Change leading dimension of U
*
                     IF( LDA.GT.LDU ) THEN
                        DO 470 J2 = N - 1, 1, -1
                           DO 460 J1 = N, 1, -1
                              URE( J1+LDA*J2 ) = URE( J1+LDU*J2 )
                              UIM( J1+LDA*J2 ) = UIM( J1+LDU*J2 )
  460                      CONTINUE
  470                   CONTINUE
                        LDU = LDA
                     ELSE IF( LDA.LT.LDU ) THEN
                        DO 490 J2 = 1, N - 1
                           DO 480 J1 = 1, N
                              URE( J1+LDA*J2 ) = URE( J1+LDU*J2 )
                              UIM( J1+LDA*J2 ) = UIM( J1+LDU*J2 )
  480                      CONTINUE
  490                   CONTINUE
                        LDU = LDA
                     END IF
*
*                    Time IMTQL2 + HTRIBK
*
                     IC = 0
                     OPS = ZERO
                     S1 = SECOND( )
  500                CONTINUE
                     CALL SCOPY( N, D, 1, RWORK, 1 )
                     CALL SCOPY( N-1, E, 1, RWORK( LDA+1 ), 1 )
                     CALL SLAZRO( N, N, ZERO, ONE, ZRE, LDA )
                     CALL IMTQL2( LDA, N, RWORK, RWORK( LDA+1 ), ZRE,
     $                            IINFO )
                     IF( IINFO.NE.0 ) THEN
                        WRITE( NOUT, FMT = 9997 )SUBNAM( 9 ), IINFO, N,
     $                     ITYPE, IPAR, IOLDSD
                        INFO = ABS( IINFO )
                        GO TO 530
                     END IF
                     CALL HTRIBK( LDA, N, URE, UIM, TAURE, N, ZRE, ZIM )
                     S2 = SECOND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN )
     $                  GO TO 500
*
*                    Subtract the time used in copying
*
                     S1 = SECOND( )
                     DO 510 J = 1, IC
                        CALL SCOPY( N, D, 1, RWORK, 1 )
                        CALL SCOPY( N-1, E, 1, RWORK( LDA+1 ), 1 )
                        CALL SLAZRO( N, N, ZERO, ONE, ZRE, LDA )
  510                CONTINUE
                     S2 = SECOND( )
                     UNTIME = S2 - S1
*
                     TIMES( IPAR, ITYPE, IN, 9 ) = MAX( TIME-UNTIME,
     $                  ZERO ) / REAL( IC )
                     OPCNTS( IPAR, ITYPE, IN, 9 ) = OPS / REAL( IC )
                  ELSE
                     OPCNTS( IPAR, ITYPE, IN, 9 ) = OPCNTS( LASTL,
     $                  ITYPE, IN, 9 )
                     TIMES( IPAR, ITYPE, IN, 9 ) = TIMES( LASTL, ITYPE,
     $                  IN, 9 )
                  END IF
  520          CONTINUE
            END IF
*
  530    CONTINUE
  540 CONTINUE
*
*-----------------------------------------------------------------------
*
*     Print a table of results for each timed routine.
*
      DO 550 ISUB = 1, NSUBS
         IF( TIMSUB( ISUB ) ) THEN
            CALL SPRTBE( SUBNAM( ISUB ), MTYPES, DOTYPE, NSIZES, NN,
     $                   INPARM( ISUB ), PNAMES, NPARMS, LDAS, NNB,
     $                   IDUMMA, IDUMMA, OPCNTS( 1, 1, 1, ISUB ), LDO1,
     $                   LDO2, TIMES( 1, 1, 1, ISUB ), LDT1, LDT2,
     $                   RWORK, LLWORK, NOUT )
         END IF
  550 CONTINUE
*
 9997 FORMAT( ' CTIM22: ', A, ' returned INFO=', I6, '.', / 9X, 'N=',
     $      I6, ', ITYPE=', I6, ', IPAR=', I6, ', ISEED=(',
     $      3( I5, ',' ), I5, ')' )
*
      RETURN
*
*     End of CTIM22
*
      END
