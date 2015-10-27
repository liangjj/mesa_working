      SUBROUTINE DDRVSG( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH,
     $                   NOUNIT, A, LDA, B, LDB, SD, SE, D1, D2, D3, D4,
     $                   D5, WA1, WA2, WA3, WR, U, LDU, BB, V, TAU, Z,
     $                   UZ, WORK, NWORK, IWORK, RESULT, INFO )
*
*  -- LAPACK test routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDB, LDU, NOUNIT, NSIZES, NTYPES,
     $                   NWORK
      DOUBLE PRECISION   THRESH
*     ..
*     .. Array Arguments ..
      LOGICAL            DOTYPE( * )
      INTEGER            ISEED( 4 ), IWORK( * ), NN( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), BB( LDB, * ),
     $                   D1( * ), D2( * ), D3( * ), D4( * ), D5( * ),
     $                   RESULT( * ), SD( * ), SE( * ), TAU( * ),
     $                   U( * ), UZ( * ), V( LDU, * ), WA1( * ),
     $                   WA2( * ), WA3( * ), WORK( * ), WR( * ),
     $                   Z( LDU, * )
*     ..
*
*  Purpose
*  =======
*
*  DDRVSG checks the symmetric generalized eigenvalue problem routines.
*
*     DSYGST reduces a DOUBLE PRECISION generalized symmetric
*     eigenproblem, of the form A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or
*     B*A*x=(lambda)*x, to the standard symmetric eigenvalue
*     problem. Here A and B are assumed to be symmetric and B
*     is also positive definite.
*
*     DSYTRD factors A as  U S U' , where ' means transpose,
*     S is symmetric tridiagonal, and U is orthogonal.
*     U is represented as a product of Householder
*     transformations, whose vectors are stored in the first
*     n-1 columns of V, and whose scale factors are in TAU.
*
*     DORGTR constructs the matrix U from the contents of V and TAU.
*
*     DSTEQR factors S as  Z D1 Z' , where Z is the orthogonal
*     matrix of eigenvectors and D1 is a diagonal matrix with
*     the eigenvalues on the diagonal.  D2 is the matrix of
*     eigenvalues computed when Z is not computed.
*
*     DSTERF computes D3, the matrix of eigenvalues, by the
*     PWK method, which does not yield eigenvectors.
*
*       When DDRVSG is called, a number of matrix "sizes" ("n's") and a
*       number of matrix "types" are specified.  For each size ("n")
*       and each type of matrix, one matrix will be generated and used
*       to test the symmetric generalized eigenroutines.
*       For each matrix, one test will be performed:
*
*       (1)     | A - B V S V' | / ( |A| n ulp )  (i.e., using the
*                                                vectors in V)
*
*       The "sizes" are specified by an array NN(1:NSIZES); the value of
*       each element NN(j) specifies one size.
*       The "types" are specified by a logical array DOTYPE( 1:NTYPES );
*       if DOTYPE(j) is .TRUE., then matrix type "j" will be generated.
*       Currently, the list of possible types is:
*
*       (1)  The zero matrix.
*       (2)  The identity matrix.
*
*       (3)  A diagonal matrix with evenly spaced entries
*            1, ..., ULP  and random signs.
*            (ULP = (first number larger than 1) - 1 )
*       (4)  A diagonal matrix with geometrically spaced entries
*            1, ..., ULP  and random signs.
*       (5)  A diagonal matrix with "clustered" entries 1, ULP, ..., ULP
*            and random signs.
*
*       (6)  Same as (4), but multiplied by SQRT( overflow threshold )
*       (7)  Same as (4), but multiplied by SQRT( underflow threshold )
*
*       (8)  A matrix of the form  U' D U, where U is orthogonal and
*            D has evenly spaced entries 1, ..., ULP with random signs
*            on the diagonal.
*
*       (9)  A matrix of the form  U' D U, where U is orthogonal and
*            D has geometrically spaced entries 1, ..., ULP with random
*            signs on the diagonal.
*
*       (10) A matrix of the form  U' D U, where U is orthogonal and
*            D has "clustered" entries 1, ULP,..., ULP with random
*            signs on the diagonal.
*
*       (11) Same as (8), but multiplied by SQRT( overflow threshold )
*       (12) Same as (8), but multiplied by SQRT( underflow threshold )
*
*       (13) Symmetric matrix with random entries chosen from (-1,1).
*       (14) Same as (13), but multiplied by SQRT( overflow threshold )
*       (15) Same as (13), but multiplied by SQRT( underflow threshold )
*       (16) Same as (8), but diagonal elements are all positive.
*       (17) Same as (9), but diagonal elements are all positive.
*       (18) Same as (10), but diagonal elements are all positive.
*       (19) Same as (16), but multiplied by SQRT( overflow threshold )
*       (20) Same as (16), but multiplied by SQRT( underflow threshold )
*       (21) A diagonally dominant tridiagonal matrix with geometrically
*            spaced diagonal entries 1, ..., ULP.
*
*  Arguments
*  =========
*
*  NSIZES  INTEGER
*          The number of sizes of matrices to use.  If it is zero,
*          DDRVSG does nothing.  It must be at least zero.
*          Not modified.
*
*  NN      INTEGER array, dimension (NSIZES)
*          An array containing the sizes to be used for the matrices.
*          Zero values will be skipped.  The values must be at least
*          zero.
*          Not modified.
*
*  NTYPES  INTEGER
*          The number of elements in DOTYPE.   If it is zero, DDRVSG
*          does nothing.  It must be at least zero.  If it is MAXTYP+1
*          and NSIZES is 1, then an additional type, MAXTYP+1 is
*          defined, which is to use whatever matrix is in A.  This
*          is only useful if DOTYPE(1:MAXTYP) is .FALSE. and
*          DOTYPE(MAXTYP+1) is .TRUE. .
*          Not modified.
*
*  DOTYPE  LOGICAL array, dimension (NTYPES)
*          If DOTYPE(j) is .TRUE., then for each size in NN a
*          matrix of that size and of type j will be generated.
*          If NTYPES is smaller than the maximum number of types
*          defined (PARAMETER MAXTYP), then types NTYPES+1 through
*          MAXTYP will not be generated.  If NTYPES is larger
*          than MAXTYP, DOTYPE(MAXTYP+1) through DOTYPE(NTYPES)
*          will be ignored.
*          Not modified.
*
*  ISEED   INTEGER array, dimension (4)
*          On entry ISEED specifies the seed of the random number
*          generator. The array elements should be between 0 and 4095;
*          if not they will be reduced mod 4096.  Also, ISEED(4) must
*          be odd.  The random number generator uses a linear
*          congruential sequence limited to small integers, and so
*          should produce machine independent random numbers. The
*          values of ISEED are changed on exit, and can be used in the
*          next call to DDRVSG to continue the same random number
*          sequence.
*          Modified.
*
*  THRESH  DOUBLE PRECISION
*          A test will count as "failed" if the "error", computed as
*          described above, exceeds THRESH.  Note that the error
*          is scaled to be O(1), so THRESH should be a reasonably
*          small multiple of 1, e.g., 10 or 100.  In particular,
*          it should not depend on the precision (single vs. double)
*          or the size of the matrix.  It must be at least zero.
*          Not modified.
*
*  NOUNIT  INTEGER
*          The FORTRAN unit number for printing out error messages
*          (e.g., if a routine returns IINFO not equal to 0.)
*          Not modified.
*
*  A       DOUBLE PRECISION array, dimension (LDA , max(NN))
*          Used to hold the matrix whose eigenvalues are to be
*          computed.  On exit, A contains the last matrix actually
*          used.
*          Modified.
*
*  LDA     INTEGER
*          The leading dimension of A.  It must be at
*          least 1 and at least max( NN ).
*          Not modified.
*
*  B       DOUBLE PRECISION array, dimension (LDB , max(NN))
*          Used to hold the symmetric positive definite matrix for
*          the generailzed problem.
*          On exit, B contains the last matrix actually
*          used.
*          Modified.
*
*  LDB     INTEGER
*          The leading dimension of B.  It must be at
*          least 1 and at least max( NN ).
*          Not modified.
*
*  SD      DOUBLE PRECISION array, dimension (max(NN))
*          The diagonal of the tridiagonal matrix computed by DSYTRD.
*          On exit, SD and SE contain the tridiagonal form of the
*          matrix in A.
*          Modified.
*
*  SE      DOUBLE PRECISION array, dimension (max(NN))
*          The off-diagonal of the tridiagonal matrix computed by
*          DSYTRD.  On exit, SD and SE contain the tridiagonal form of
*          the matrix in A.
*          Modified.
*
*  D1      DOUBLE PRECISION array, dimension (max(NN))
*          The eigenvalues of A, as computed by DSTEQR simlutaneously
*          with Z.  On exit, the eigenvalues in D1 correspond with the
*          matrix in A.
*          Modified.
*
*  D2      DOUBLE PRECISION array, dimension (max(NN))
*          The eigenvalues of A, as computed by DSTEQR if Z is not
*          computed.  On exit, the eigenvalues in D2 correspond with
*          the matrix in A.
*          Modified.
*
*  D3      DOUBLE PRECISION array, dimension (max(NN))
*          The eigenvalues of A, as computed by DSTERF.  On exit, the
*          eigenvalues in D3 correspond with the matrix in A.
*          Modified.
*
*  U       DOUBLE PRECISION array, dimension (LDU, max(NN))
*          The orthogonal matrix computed by DSYTRD + DORGTR.
*          Modified.
*
*  LDU     INTEGER
*          The leading dimension of U, Z, V, and UZ.  It must be at
*          least 1 and at least max( NN ).
*          Not modified.
*
*  V       DOUBLE PRECISION array, dimension (LDU, max(NN))
*          The Housholder vectors computed by DSYTRD in reducing A to
*          tridiagonal form.
*          Modified.
*
*  TAU     DOUBLE PRECISION array, dimension (max(NN))
*          The Householder factors computed by DSYTRD in reducing A
*          to tridiagonal form.
*          Modified.
*
*  Z       DOUBLE PRECISION array, dimension (LDU, max(NN))
*          The orthogonal matrix of eigenvectors computed by DSTEQR,
*          DPTEQR, and DSTEIN.
*          Modified.
*
*  UZ      DOUBLE PRECISION array, dimension (LDU, max(NN))
*          The product of U times Z.
*          Modified.
*
*  WORK    DOUBLE PRECISION array, dimension (NWORK)
*          Workspace.
*          Modified.
*
*  NWORK   INTEGER
*          The number of entries in WORK.  This must be at least
*          2*max( NN(j), 2 )**2.
*          Not modified.
*
*  IWORK   INTEGER array, dimension (3*max(NN))
*          Workspace.
*          Modified.
*
*  RESULT  DOUBLE PRECISION array, dimension (17)
*          The values computed by the tests described above.
*          The values are currently limited to 1/ulp, to avoid
*          overflow.
*          Modified.
*
*  INFO    INTEGER
*          If 0, then everything ran OK.
*           -1: NSIZES < 0
*           -2: Some NN(j) < 0
*           -3: NTYPES < 0
*           -5: THRESH < 0
*           -9: LDA < 1 or LDA < NMAX, where NMAX is max( NN(j) ).
*          -16: LDU < 1 or LDU < NMAX.
*          -21: NWORK too small.
*          If  DLATMR, SLATMS, DSYTRD, DORGTR, DSTEQR, SSTERF,
*              or DORMTR returns an error code, the
*              absolute value of it is returned.
*          Modified.
*
*-----------------------------------------------------------------------
*
*       Some Local Variables and Parameters:
*       ---- ----- --------- --- ----------
*       ZERO, ONE       Real 0 and 1.
*       MAXTYP          The number of types defined.
*       NTEST           The number of tests performed, or which can
*                       be performed so far, for the current matrix.
*       NTESTT          The total number of tests performed so far.
*       NMAX            Largest value in NN.
*       NMATS           The number of matrices generated so far.
*       NERRS           The number of tests which have exceeded THRESH
*                       so far (computed by DLAFTS).
*       COND, IMODE     Values to be passed to the matrix generators.
*       ANORM           Norm of A; passed to matrix generators.
*
*       OVFL, UNFL      Overflow and underflow thresholds.
*       ULP, ULPINV     Finest relative precision and its inverse.
*       RTOVFL, RTUNFL  Square roots of the previous 2 values.
*               The following four arrays decode JTYPE:
*       KTYPE(j)        The general type (1-10) for type "j".
*       KMODE(j)        The MODE value to be passed to the matrix
*                       generator for type "j".
*       KMAGN(j)        The order of magnitude ( O(1),
*                       O(overflow^(1/2) ), O(underflow^(1/2) )
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, TEN
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,
     $                   TEN = 10.0D0 )
      DOUBLE PRECISION   HALF
      PARAMETER          ( HALF = ONE / TWO )
      INTEGER            MAXTYP
      PARAMETER          ( MAXTYP = 21 )
*     ..
*     .. Local Scalars ..
      LOGICAL            BADNN
      CHARACTER          UPLO
      INTEGER            I, IBTYPE, IBUPLO, IINFO, IJ, IMODE, ITYPE, J,
     $                   JCOL, JSIZE, JTYPE, MTYPES, N, NERRS, NERRS2,
     $                   NMATS, NMAX, NTEST, NTEST2, NTESTT
      DOUBLE PRECISION   ANINV, ANORM, COND, OVFL, RTOVFL, RTUNFL,
     $                   TEMP1, ULP, ULPINV, UNFL
*     ..
*     .. Local Arrays ..
      INTEGER            IDUMMA( 1 ), IOLDSD( 4 ), KMAGN( MAXTYP ),
     $                   KMODE( MAXTYP ), KTYPE( MAXTYP )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLABAD, DLACPY, DLAFTS, DLASUM, DLATMR, DLATMS,
     $                   DLAZRO, DSGT01, DSPGV, DSYGV, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN, SQRT
*     ..
*     .. Data statements ..
      DATA               KTYPE / 1, 2, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 8,
     $                   8, 8, 5, 5, 5, 5, 5, 4 /
      DATA               KMAGN / 1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1,
     $                   2, 3, 1, 1, 1, 2, 3, 1 /
      DATA               KMODE / 0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0,
     $                   0, 0, 4, 3, 1, 4, 4, 3 /
*     ..
*     .. Executable Statements ..
*
*     1)      Check for errors
*
      NTESTT = 0
      NTEST2 = 0
      INFO = 0
*
      BADNN = .FALSE.
      NMAX = 1
      DO 10 J = 1, NSIZES
         NMAX = MAX( NMAX, NN( J ) )
         IF( NN( J ).LT.0 )
     $      BADNN = .TRUE.
   10 CONTINUE
*
*     Check for errors
*
      IF( NSIZES.LT.0 ) THEN
         INFO = -1
      ELSE IF( BADNN ) THEN
         INFO = -2
      ELSE IF( NTYPES.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.NMAX ) THEN
         INFO = -9
      ELSE IF( LDU.LT.NMAX ) THEN
         INFO = -16
      ELSE IF( 2*MAX( 2, NMAX )**2.GT.NWORK ) THEN
         INFO = -21
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DDRVSG', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( NSIZES.EQ.0 .OR. NTYPES.EQ.0 )
     $   RETURN
*
*     More Important constants
*
      UNFL = DLAMCH( 'Safe minimum' )
      OVFL = DLAMCH( 'Overflow' )
      CALL DLABAD( UNFL, OVFL )
      ULP = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' )
      ULPINV = ONE / ULP
      RTUNFL = SQRT( UNFL )
      RTOVFL = SQRT( OVFL )
*
*     Loop over sizes, types
*
      NERRS = 0
      NERRS2 = 0
      NMATS = 0
*
      DO 170 JSIZE = 1, NSIZES
         N = NN( JSIZE )
         ANINV = ONE / DBLE( MAX( 1, N ) )
*
         IF( NSIZES.NE.1 ) THEN
            MTYPES = MIN( MAXTYP, NTYPES )
         ELSE
            MTYPES = MIN( MAXTYP+1, NTYPES )
         END IF
*
         DO 160 JTYPE = 1, MTYPES
            IF( .NOT.DOTYPE( JTYPE ) )
     $         GO TO 160
*
            NMATS = NMATS + 1
            NTEST = 0
*
            DO 20 J = 1, 4
               IOLDSD( J ) = ISEED( J )
   20       CONTINUE
*
*           2)      Compute "A"
*
*                   Control parameters:
*
*               KMAGN  KMODE        KTYPE
*           =1  O(1)   clustered 1  zero
*           =2  large  clustered 2  identity
*           =3  small  exponential  (none)
*           =4         arithmetic   diagonal, (w/ eigenvalues)
*           =5         random log   symmetric, w/ eigenvalues
*           =6         random       (none)
*           =7                      random diagonal
*           =8                      random symmetric
*
            IF( MTYPES.GT.MAXTYP )
     $         GO TO 90
*
            ITYPE = KTYPE( JTYPE )
            IMODE = KMODE( JTYPE )
*
*           Compute norm
*
            GO TO ( 30, 40, 50 )KMAGN( JTYPE )
*
   30       CONTINUE
            ANORM = ONE
            GO TO 60
*
   40       CONTINUE
            ANORM = ( RTOVFL*ULP )*ANINV
            GO TO 60
*
   50       CONTINUE
            ANORM = RTUNFL*N*ULPINV
            GO TO 60
*
   60       CONTINUE
*
            CALL DLAZRO( LDA, N, ZERO, ZERO, A, LDA )
            IINFO = 0
            COND = ULPINV
*
*           Special Matrices -- Identity & Jordan block
*
            IF( ITYPE.EQ.1 ) THEN
*
*              Zero
*
               IINFO = 0
*
            ELSE IF( ITYPE.EQ.2 ) THEN
*
*              Identity
*
               DO 70 JCOL = 1, N
                  A( JCOL, JCOL ) = ANORM
   70          CONTINUE
*
            ELSE IF( ITYPE.EQ.4 ) THEN
*
*              Diagonal Matrix, [Eigen]values Specified
*
               IF( JTYPE.NE.21 ) THEN
                  CALL DLATMS( N, N, 'S', ISEED, 'S', WORK, IMODE, COND,
     $                         ANORM, 0, 0, 'N', A, LDA, WORK( N+1 ),
     $                         IINFO )
               ELSE
                  COND = ULPINV / ( TEN*N )
                  CALL DLATMS( N, N, 'S', ISEED, 'P', WORK, IMODE, COND,
     $                         ANORM, 1, 1, 'N', A, LDA, WORK( N+1 ),
     $                         IINFO )
                  DO 80 I = 2, N
                     TEMP1 = ABS( A( I-1, I ) ) /
     $                       SQRT( ABS( A( I-1, I-1 )*A( I, I ) ) )
                     IF( TEMP1.GT.HALF ) THEN
                        A( I-1, I ) = HALF*SQRT( ABS( A( I-1,
     $                                I-1 )*A( I, I ) ) )
                        A( I, I-1 ) = A( I-1, I )
                     END IF
   80             CONTINUE
               END IF
*
            ELSE IF( ITYPE.EQ.5 ) THEN
*
*              Symmetric, eigenvalues specified
*
               IF( JTYPE.LE.15 ) THEN
                  CALL DLATMS( N, N, 'S', ISEED, 'S', WORK, IMODE, COND,
     $                         ANORM, N, N, 'N', A, LDA, WORK( N+1 ),
     $                         IINFO )
               ELSE
                  COND = ULPINV / ( TEN*N )
                  CALL DLATMS( N, N, 'S', ISEED, 'P', WORK, IMODE, COND,
     $                         ANORM, N, N, 'N', A, LDA, WORK( N+1 ),
     $                         IINFO )
               END IF
*
            ELSE IF( ITYPE.EQ.7 ) THEN
*
*              Diagonal, random eigenvalues
*
               CALL DLATMR( N, N, 'S', ISEED, 'S', WORK, 6, ONE, ONE,
     $                      'T', 'N', WORK( N+1 ), 1, ONE,
     $                      WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, 0, 0,
     $                      ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO )
*
            ELSE IF( ITYPE.EQ.8 ) THEN
*
*              Symmetric, random eigenvalues
*
               CALL DLATMR( N, N, 'S', ISEED, 'S', WORK, 6, ONE, ONE,
     $                      'T', 'N', WORK( N+1 ), 1, ONE,
     $                      WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, N, N,
     $                      ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO )
*
            ELSE
*
               IINFO = 1
            END IF
*
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'Generator', IINFO, N, JTYPE,
     $            IOLDSD
               INFO = ABS( IINFO )
               RETURN
            END IF
*
   90       CONTINUE
*
*           3)  Call DSYGV  or DSPGV  to compute S and U, do tests.
*
*               loop over the three generalized problems
*                 IBTYPE = 1: A*x = (lambda)*B*x
*                 IBTYPE = 2: A*B*x = (lambda)*x
*                 IBTYPE = 3: B*A*x = (lambda)*x
*
            DO 150 IBTYPE = 1, 3
*
*              loop over the setting UPLO
*
               DO 140 IBUPLO = 1, 2
                  IF( IBUPLO.EQ.1 )
     $               UPLO = 'Upper'
                  IF( IBUPLO.EQ.2 )
     $               UPLO = 'Lower'
*
                  CALL DLACPY( ' ', N, N, A, LDA, V, LDU )
*
                  NTEST = 1
*
                  CALL DLATMS( N, N, 'U', ISEED, 'P', WORK, 5, TEN, ONE,
     $                         N, N, UPLO, B, LDB, WORK( N+1 ), IINFO )
                  CALL DLACPY( UPLO, N, N, B, LDB, BB, LDB )
                  CALL DSYGV( IBTYPE, 'V', UPLO, N, V, LDU, B, LDB, D1,
     $                        WORK, NWORK, IINFO )
*
*                 Do Test
*
                  NTEST = 1
                  CALL DSGT01( IBTYPE, UPLO, N, A, LDA, BB, LDB, V, LDU,
     $                         D1, UZ, RESULT( 1 ) )
                  NTESTT = NTESTT + NTEST
                  CALL DLAFTS( 'DSG', N, N, JTYPE, NTEST, RESULT,
     $                         IOLDSD, THRESH, NOUNIT, NERRS )
*
*                 Copy the matrices into packed storage.
*
                  IF( LSAME( UPLO, 'U' ) ) THEN
                     IJ = 1
                     DO 110 J = 1, N
                        DO 100 I = 1, J
                           UZ( IJ ) = A( I, J )
                           U( IJ ) = BB( I, J )
                           IJ = IJ + 1
  100                   CONTINUE
  110                CONTINUE
                  ELSE
                     IJ = 1
                     DO 130 J = 1, N
                        DO 120 I = J, N
                           UZ( IJ ) = A( I, J )
                           U( IJ ) = BB( I, J )
                           IJ = IJ + 1
  120                   CONTINUE
  130                CONTINUE
                  END IF
*
                  CALL DSPGV( IBTYPE, 'V', UPLO, N, UZ, U, D1, V, LDU,
     $                        WORK, IINFO )
*
*                 Do Test
*
                  NTEST = 1
                  CALL DSGT01( IBTYPE, UPLO, N, A, LDA, BB, LDB, V, LDU,
     $                         D1, UZ, RESULT( 1 ) )
*
*                 End of Loop -- Check for RESULT(j) > THRESH
*
                  NTEST2 = NTEST2 + NTEST
                  CALL DLAFTS( 'DSG', N, N, JTYPE, NTEST, RESULT,
     $                         IOLDSD, THRESH, NOUNIT, NERRS2 )
*
  140          CONTINUE
  150       CONTINUE
  160    CONTINUE
  170 CONTINUE
*
*     Summary
*
      CALL DLASUM( 'DSG', NOUNIT, NERRS+NERRS2, NTESTT+NTEST2 )
*
 9999 FORMAT( ' DDRVSG: ', A, ' returned INFO=', I6, '.', / 9X, 'N=',
     $      I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' )
*
      RETURN
*
*     End of DDRVSG
*
      END
