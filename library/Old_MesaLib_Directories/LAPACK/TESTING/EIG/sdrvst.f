      SUBROUTINE SDRVST( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH,
     $                   NOUNIT, A, LDA, SD, SE, D1, D2, D3, D4, D5,
     $                   WA1, WA2, WA3, WR, U, LDU, V, TAU, Z, UZ, WORK,
     $                   NWORK, IWORK, RESULT, INFO )
*
*  -- LAPACK test routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDU, NOUNIT, NSIZES, NTYPES, NWORK
      REAL               THRESH
*     ..
*     .. Array Arguments ..
      LOGICAL            DOTYPE( * )
      INTEGER            ISEED( 4 ), IWORK( * ), NN( * )
      REAL               A( LDA, * ), D1( * ), D2( * ), D3( * ),
     $                   D4( * ), D5( * ), RESULT( * ), SD( * ),
     $                   SE( * ), TAU( * ), U( LDU, * ), UZ( LDU, * ),
     $                   V( LDU, * ), WA1( * ), WA2( * ), WA3( * ),
     $                   WORK( * ), WR( * ), Z( LDU, * )
*     ..
*
*  Purpose
*  =======
*
*       SDRVST  checks the symmetric eigenvalue problem drivers.
*
*               SSTEV computes all eigenvalues and, optionally,
*               eigenvectors of a real symmetric tridiagonal matrix.
*
*               SSTEVX computes selected eigenvalues and, optionally,
*               eigenvectors of a real symmetric tridiagonal matrix.
*
*               SSYEV computes all eigenvalues and, optionally,
*               eigenvectors of a real symmetric matrix.
*
*               SSYEVX computes selected eigenvalues and, optionally,
*               eigenvectors of a real symmetric matrix.
*
*               SSPEV computes all eigenvalues and, optionally,
*               eigenvectors of a real symmetric matrix in packed
*               storage.
*
*               SSPEVX computes selected eigenvalues and, optionally,
*               eigenvectors of a real symmetric matrix in packed
*               storage.
*
*               SSBEV computes all eigenvalues and, optionally,
*               eigenvectors of a real symmetric band matrix.
*
*               SSBEVX computes selected eigenvalues and, optionally,
*               eigenvectors of a real symmetric band matrix.
*
*       When SDRVST is called, a number of matrix "sizes" ("n's") and a
*       number of matrix "types" are specified.  For each size ("n")
*       and each type of matrix, one matrix will be generated and used
*       to test the appropriate drivers.  For each matrix and each
*       driver routine called, the following tests will be performed:
*
*       (1)     | A - Z D Z' | / ( |A| n ulp )
*
*       (2)     | I - Z Z' | / ( n ulp )
*
*       (3)     | D1 - D2 | / ( |D1| ulp )
*
*       where Z is the matrix of eigenvectors returned when the
*       eigenvector option is given and D1 and D2 are the eigenvalues
*       returned with and without the eigenvector option.
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
*       (16) A band matrix with half bandwidth randomly chosen between
*            0 and N-1, with evenly spaced eigenvalues 1, ..., ULP
*            with random signs.
*       (17) Same as (16), but multiplied by SQRT( overflow threshold )
*       (18) Same as (16), but multiplied by SQRT( underflow threshold )
*
*  Arguments
*  =========
*
*  NSIZES  INTEGER
*          The number of sizes of matrices to use.  If it is zero,
*          SDRVST does nothing.  It must be at least zero.
*          Not modified.
*
*  NN      INTEGER array, dimension (NSIZES)
*          An array containing the sizes to be used for the matrices.
*          Zero values will be skipped.  The values must be at least
*          zero.
*          Not modified.
*
*  NTYPES  INTEGER
*          The number of elements in DOTYPE.   If it is zero, SDRVST
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
*          next call to SDRVST to continue the same random number
*          sequence.
*          Modified.
*
*  THRESH  REAL
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
*  A       REAL array, dimension (LDA , max(NN))
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
*  SD      REAL array, dimension (max(NN))
*          The diagonal of the tridiagonal matrix computed by SSYTRD.
*          On exit, SD and SE contain the tridiagonal form of the
*          matrix in A.
*          Modified.
*
*  SE      REAL array, dimension (max(NN))
*          The off-diagonal of the tridiagonal matrix computed by
*          SSYTRD.  On exit, SD and SE contain the tridiagonal form of
*          the matrix in A.
*          Modified.
*
*  D1      REAL array, dimension (max(NN))
*          The eigenvalues of A, as computed by SSTEQR simlutaneously
*          with Z.  On exit, the eigenvalues in D1 correspond with the
*          matrix in A.
*          Modified.
*
*  D2      REAL array, dimension (max(NN))
*          The eigenvalues of A, as computed by SSTEQR if Z is not
*          computed.  On exit, the eigenvalues in D2 correspond with
*          the matrix in A.
*          Modified.
*
*  D3      REAL array, dimension (max(NN))
*          The eigenvalues of A, as computed by SSTERF.  On exit, the
*          eigenvalues in D3 correspond with the matrix in A.
*          Modified.
*
*  U       REAL array, dimension (LDU, max(NN))
*          The orthogonal matrix computed by SSYTRD + SORGTR.
*          Modified.
*
*  LDU     INTEGER
*          The leading dimension of U, Z, V, and UZ.  It must be at
*          least 1 and at least max( NN ).
*          Not modified.
*
*  V       REAL array, dimension (LDU, max(NN))
*          The Housholder vectors computed by SSYTRD in reducing A to
*          tridiagonal form.
*          Modified.
*
*  TAU     REAL array, dimension (max(NN))
*          The Householder factors computed by SSYTRD in reducing A
*          to tridiagonal form.
*          Modified.
*
*  Z       REAL array, dimension (LDU, max(NN))
*          The orthogonal matrix of eigenvectors computed by SSTEQR,
*          SPTEQR, and SSTEIN.
*          Modified.
*
*  UZ      REAL array, dimension (LDU, max(NN))
*          The product of U times Z.
*          Modified.
*
*  WORK    REAL array, dimension (NWORK)
*          Workspace.
*          Modified.
*
*  NWORK   INTEGER
*          The number of entries in WORK.  This must be at least
*          2*max( NN(j), 2 )**2.
*          Not modified.
*
*  IWORK   INTEGER array, dimension (6*max(NN))
*          Workspace.
*          Modified.
*
*  RESULT  REAL array, dimension (84)
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
*          If  SLATMR, SLATMS, SSYTRD, SORGTR, SSTEQR, SSTERF,
*              or SORMTR returns an error code, the
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
*                       so far (computed by SLAFTS).
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
      REAL               ZERO, ONE, TEN
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0, TEN = 10.0E0 )
      REAL               HALF
      PARAMETER          ( HALF = 0.5E0 )
      INTEGER            MAXTYP
      PARAMETER          ( MAXTYP = 18 )
*     ..
*     .. Local Scalars ..
      LOGICAL            BADNN
      CHARACTER          UPLO
      INTEGER            I, IDIAG, IHBW, IINFO, IL, IMODE, INDX, IROW,
     $                   ITEMP, ITYPE, IU, IUPLO, J, J1, J2, JCOL,
     $                   JSIZE, JTYPE, KD, M, M2, M3, MTYPES, N, NERRS,
     $                   NMATS, NMAX, NTEST, NTESTT
      REAL               ABSTOL, ANINV, ANORM, COND, OVFL, RTOVFL,
     $                   RTUNFL, TEMP1, TEMP2, TEMP3, ULP, ULPINV, UNFL,
     $                   VL, VU
*     ..
*     .. Local Arrays ..
      INTEGER            IDUMMA( 1 ), IOLDSD( 4 ), ISEED2( 4 ),
     $                   ISEED3( 4 ), KMAGN( MAXTYP ), KMODE( MAXTYP ),
     $                   KTYPE( MAXTYP )
*     ..
*     .. External Functions ..
      REAL               SLAMCH, SLARND, SSXT1
      EXTERNAL           SLAMCH, SLARND, SSXT1
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLABAD, SLACPY, SLAFTS, ALASVM, SLATMR, SLATMS,
     $                   SLAZRO, SSBEV, SSBEVX, SSPEV, SSPEVX, SSTEV,
     $                   SSTEVX, SSTT21, SSTT22, SSYEV, SSYEVX, SSYT21,
     $                   SSYT22, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, REAL, SQRT
*     ..
*     .. Data statements ..
      DATA               KTYPE / 1, 2, 5*4, 5*5, 3*8, 3*9 /
      DATA               KMAGN / 2*1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1,
     $                   2, 3, 1, 2, 3 /
      DATA               KMODE / 2*0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0,
     $                   0, 0, 4, 4, 4 /
*     ..
*     .. Executable Statements ..
*
*     1)      Check for errors
*
      NTESTT = 0
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
         CALL XERBLA( 'SDRVST', -INFO )
         RETURN
      END IF
*
*     Quick return if nothing to do
*
      IF( NSIZES.EQ.0 .OR. NTYPES.EQ.0 )
     $   RETURN
*
*     More Important constants
*
      UNFL = SLAMCH( 'Safe minimum' )
      OVFL = SLAMCH( 'Overflow' )
      CALL SLABAD( UNFL, OVFL )
      ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' )
      ULPINV = ONE / ULP
      RTUNFL = SQRT( UNFL )
      RTOVFL = SQRT( OVFL )
*
*     Loop over sizes, types
*
      DO 20 I = 1, 4
         ISEED2( I ) = ISEED( I )
         ISEED3( I ) = ISEED( I )
   20 CONTINUE
*
      NERRS = 0
      NMATS = 0
*
      DO 1230 JSIZE = 1, NSIZES
         N = NN( JSIZE )
         ANINV = ONE / REAL( MAX( 1, N ) )
*
         IF( NSIZES.NE.1 ) THEN
            MTYPES = MIN( MAXTYP, NTYPES )
         ELSE
            MTYPES = MIN( MAXTYP+1, NTYPES )
         END IF
*
         DO 1220 JTYPE = 1, MTYPES
            IF( .NOT.DOTYPE( JTYPE ) )
     $         GO TO 1220
            NMATS = NMATS + 1
            NTEST = 0
*
            DO 30 J = 1, 4
               IOLDSD( J ) = ISEED( J )
   30       CONTINUE
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
*           =9                      band symmetric, w/ eigenvalues
*
            IF( MTYPES.GT.MAXTYP )
     $         GO TO 110
*
            ITYPE = KTYPE( JTYPE )
            IMODE = KMODE( JTYPE )
*
*           Compute norm
*
            GO TO ( 40, 50, 60 )KMAGN( JTYPE )
*
   40       CONTINUE
            ANORM = ONE
            GO TO 70
*
   50       CONTINUE
            ANORM = ( RTOVFL*ULP )*ANINV
            GO TO 70
*
   60       CONTINUE
            ANORM = RTUNFL*N*ULPINV
            GO TO 70
*
   70       CONTINUE
*
            CALL SLAZRO( LDA, N, ZERO, ZERO, A, LDA )
            IINFO = 0
            COND = ULPINV
*
*           Special Matrices -- Identity & Jordan block
*
*                   Zero
*
            IF( ITYPE.EQ.1 ) THEN
               IINFO = 0
*
            ELSE IF( ITYPE.EQ.2 ) THEN
*
*              Identity
*
               DO 80 JCOL = 1, N
                  A( JCOL, JCOL ) = ANORM
   80          CONTINUE
*
            ELSE IF( ITYPE.EQ.4 ) THEN
*
*              Diagonal Matrix, [Eigen]values Specified
*
               CALL SLATMS( N, N, 'S', ISEED, 'S', WORK, IMODE, COND,
     $                      ANORM, 0, 0, 'N', A, LDA, WORK( N+1 ),
     $                      IINFO )
*
            ELSE IF( ITYPE.EQ.5 ) THEN
*
*              Symmetric, eigenvalues specified
*
               CALL SLATMS( N, N, 'S', ISEED, 'S', WORK, IMODE, COND,
     $                      ANORM, N, N, 'N', A, LDA, WORK( N+1 ),
     $                      IINFO )
*
            ELSE IF( ITYPE.EQ.7 ) THEN
*
*              Diagonal, random eigenvalues
*
               CALL SLATMR( N, N, 'S', ISEED, 'S', WORK, 6, ONE, ONE,
     $                      'T', 'N', WORK( N+1 ), 1, ONE,
     $                      WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, 0, 0,
     $                      ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO )
*
            ELSE IF( ITYPE.EQ.8 ) THEN
*
*              Symmetric, random eigenvalues
*
               CALL SLATMR( N, N, 'S', ISEED, 'S', WORK, 6, ONE, ONE,
     $                      'T', 'N', WORK( N+1 ), 1, ONE,
     $                      WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, N, N,
     $                      ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO )
*
            ELSE IF( ITYPE.EQ.9 ) THEN
*
*              Symmetric banded, eigenvalues specified
*
               IHBW = ( N-1 )*SLARND( 1, ISEED3 )
               CALL SLATMS( N, N, 'S', ISEED, 'S', WORK, IMODE, COND,
     $                      ANORM, IHBW, IHBW, 'Z', U, LDU, WORK( N+1 ),
     $                      IINFO )
*
*              Store as dense matrix for most routines.
*
               CALL SLAZRO( LDA, N, ZERO, ZERO, A, LDA )
               DO 100 IDIAG = -IHBW, IHBW
                  IROW = IHBW - IDIAG + 1
                  J1 = MAX( 1, IDIAG )
                  J2 = MIN( N, N+IDIAG )
                  DO 90 J = J1, J2
                     I = J - IDIAG
                     A( I, J ) = U( IROW, J )
   90             CONTINUE
  100          CONTINUE
            ELSE
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
  110       CONTINUE
*
            ABSTOL = ZERO
            IF( N.LE.1 ) THEN
               IL = 1
               IU = N
            ELSE
               IL = 1 + ( N-1 )*SLARND( 1, ISEED2 )
               IU = 1 + ( N-1 )*SLARND( 1, ISEED2 )
               IF( IL.GT.IU ) THEN
                  ITEMP = IL
                  IL = IU
                  IU = ITEMP
               END IF
            END IF
*
*           3)      If matrix is tridiagonal, call SSTEV and SSTEVX.
*
            IF( JTYPE.LE.7 ) THEN
               NTEST = 1
               DO 120 I = 1, N
                  D1( I ) = A( I, I )
  120          CONTINUE
               DO 130 I = 1, N - 1
                  D2( I ) = A( I+1, I )
  130          CONTINUE
               CALL SSTEV( 'V', N, D1, D2, Z, LDU, WORK, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'SSTEV(V)', IINFO, N,
     $               JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( 1 ) = ULPINV
                     RESULT( 2 ) = ULPINV
                     RESULT( 3 ) = ULPINV
                     GO TO 180
                  END IF
               END IF
*
*              Do tests 1 and 2.
*
               DO 140 I = 1, N
                  D3( I ) = A( I, I )
  140          CONTINUE
               DO 150 I = 1, N-1
                  D4( I ) = A( I+1, I )
  150          CONTINUE
               CALL SSTT21( N, 0, D3, D4, D1, D2, Z, LDU, WORK,
     $                      RESULT( 1 ) )
*
               NTEST = 3
               DO 160 I = 1, N - 1
                  D4( I ) = A( I+1, I )
  160          CONTINUE
               CALL SSTEV( 'N', N, D3, D4, Z, LDU, WORK, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'SSTEV(N)', IINFO, N,
     $               JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( 3 ) = ULPINV
                     GO TO 180
                  END IF
               END IF
*
*              Do test 3.
*
               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 170 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( D1( J ) ), ABS( D3( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( D1( J )-D3( J ) ) )
  170          CONTINUE
               RESULT( 3 ) = TEMP2 / MAX( UNFL,
     $                       ULP*MAX( TEMP1, TEMP2 ) )
*
  180          CONTINUE
*
               NTEST = 4
               DO 190 I = 1, N
                  D1( I ) = A( I, I )
  190          CONTINUE
               DO 200 I = 1, N - 1
                  D2( I ) = A( I+1, I )
  200          CONTINUE
               CALL SSTEVX( 'V', 'A', N, D1, D2, VL, VU, IL, IU, ABSTOL,
     $                      M, WA1, Z, LDU, WORK, IWORK, IWORK( 5*N+1 ),
     $                      IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'SSTEVX(V,A)', IINFO, N,
     $               JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( 4 ) = ULPINV
                     RESULT( 5 ) = ULPINV
                     RESULT( 6 ) = ULPINV
                     GO TO 250
                  END IF
               END IF
               IF( N.GT.0 ) THEN
                  TEMP3 = MAX( ABS( WA1( 1 ) ), ABS( WA1( N ) ) )
               ELSE
                  TEMP3 = ZERO
               END IF
*
*              Do tests 4 and 5.
*
               DO 210 I = 1, N
                  D3( I ) = A( I, I )
  210          CONTINUE
               DO 220 I = 1, N-1
                  D4( I ) = A( I+1, I )
  220          CONTINUE
               CALL SSTT21( N, 0, D3, D4, WA1, D2, Z, LDU, WORK,
     $                      RESULT( 4 ) )
*
               NTEST = 6
               DO 230 I = 1, N - 1
                  D4( I ) = A( I+1, I )
  230          CONTINUE
               CALL SSTEVX( 'N', 'A', N, D3, D4, VL, VU, IL, IU, ABSTOL,
     $                      M2, WA2, Z, LDU, WORK, IWORK,
     $                      IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'SSTEVX(N,A)', IINFO, N,
     $               JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( 6 ) = ULPINV
                     GO TO 250
                  END IF
               END IF
*
*              Do test 6.
*
               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 240 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( WA1( J ) ), ABS( WA2( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( WA1( J )-WA2( J ) ) )
  240          CONTINUE
               RESULT( 6 ) = TEMP2 / MAX( UNFL,
     $                       ULP*MAX( TEMP1, TEMP2 ) )
*
  250          CONTINUE
*
               NTEST = 7
               DO 260 I = 1, N
                  D1( I ) = A( I, I )
  260          CONTINUE
               DO 270 I = 1, N - 1
                  D2( I ) = A( I+1, I )
  270          CONTINUE
               CALL SSTEVX( 'V', 'I', N, D1, D2, VL, VU, IL, IU, ABSTOL,
     $                      M2, WA2, Z, LDU, WORK, IWORK,
     $                      IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'SSTEVX(V,I)', IINFO, N,
     $               JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( 7 ) = ULPINV
                     RESULT( 8 ) = ULPINV
                     RESULT( 9 ) = ULPINV
                     GO TO 310
                  END IF
               END IF
*
*              Do tests 7 and 8.
*
               DO 280 I = 1, N
                  D3( I ) = A( I, I )
  280          CONTINUE
               DO 290 I = 1, N-1
                  D4( I ) = A( I+1, I )
  290          CONTINUE
               CALL SSTT22( N, M2, 0, D3, D4, WA2, D2, Z, LDU, WORK,
     $                      MAX( 1, M2 ), RESULT( 7 ) )
*
               NTEST = 9
               DO 300 I = 1, N - 1
                  D4( I ) = A( I+1, I )
  300          CONTINUE
               CALL SSTEVX( 'N', 'I', N, D3, D4, VL, VU, IL, IU, ABSTOL,
     $                      M3, WA3, Z, LDU, WORK, IWORK,
     $                      IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'SSTEVX(N,I)', IINFO, N,
     $               JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( 9 ) = ULPINV
                     GO TO 310
                  END IF
               END IF
*
*              Do test 9.
*
               TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
               TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
               RESULT( 9 ) = ( TEMP1+TEMP2 ) / MAX( UNFL, ULP*TEMP3 )
*
  310          CONTINUE
*
               NTEST = 10
               IF( N.GT.0 ) THEN
               IF( IL.NE.1 ) THEN
                  VL = WA1( IL ) - MAX( HALF*( WA1( IL )-WA1( IL-1 ) ),
     $                 TEN*ULP*TEMP3, TEN*RTUNFL )
               ELSE
                  VL = WA1( 1 ) - MAX( HALF*( WA1( N )-WA1( 1 ) ),
     $                 TEN*ULP*TEMP3, TEN*RTUNFL )
               END IF
               IF( IU.NE.N ) THEN
                  VU = WA1( IU ) + MAX( HALF*( WA1( IU+1 )-WA1( IU ) ),
     $                 TEN*ULP*TEMP3, TEN*RTUNFL )
               ELSE
                  VU = WA1( N ) + MAX( HALF*( WA1( N )-WA1( 1 ) ),
     $                 TEN*ULP*TEMP3, TEN*RTUNFL )
               END IF
               ELSE
                  VL = ZERO
                  VU = ONE
               END IF
*
               DO 320 I = 1, N
                  D1( I ) = A( I, I )
  320          CONTINUE
               DO 330 I = 1, N - 1
                  D2( I ) = A( I+1, I )
  330          CONTINUE
*
               CALL SSTEVX( 'V', 'V', N, D1, D2, VL, VU, IL, IU, ABSTOL,
     $                      M2, WA2, Z, LDU, WORK, IWORK,
     $                      IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'SSTEVX(V,V)', IINFO, N,
     $               JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( 10 ) = ULPINV
                     RESULT( 11 ) = ULPINV
                     RESULT( 12 ) = ULPINV
                     GO TO 370
                  END IF
               END IF
*
               IF( M2.EQ.0 .AND. N.GT.0 ) THEN
                  RESULT( 10 ) = ULPINV
                  RESULT( 11 ) = ULPINV
                  RESULT( 12 ) = ULPINV
                  GO TO 370
               END IF
*
*              Do tests 10 and 11.
*
               DO 340 I = 1, N
                  D3( I ) = A( I, I )
  340          CONTINUE
               DO 350 I = 1, N-1
                  D4( I ) = A( I+1, I )
  350          CONTINUE
               CALL SSTT22( N, M2, 0, D3, D4, WA2, D2, Z, LDU, WORK,
     $                      MAX( 1, M2 ), RESULT( 10 ) )
*
               NTEST = 12
               DO 360 I = 1, N - 1
                  D4( I ) = A( I+1, I )
  360          CONTINUE
               CALL SSTEVX( 'N', 'V', N, D3, D4, VL, VU, IL, IU, ABSTOL,
     $                      M3, WA3, Z, LDU, WORK, IWORK,
     $                      IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'SSTEVX(N,V)', IINFO, N,
     $               JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( 12 ) = ULPINV
                     GO TO 370
                  END IF
               END IF
*
*              Do test 12.
*
               TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
               TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
               RESULT( 12 ) = ( TEMP1+TEMP2 ) / MAX( UNFL, TEMP3*ULP )
*
  370          CONTINUE
*
            ELSE
*
               DO 380 I = 1, 12
                  RESULT( I ) = ZERO
  380          CONTINUE
               NTEST = 12
            END IF
*
*           Perform remaining tests storing upper or lower triangular
*           part of matrix.
*
            DO 1210 IUPLO = 0, 1
               IF( IUPLO.EQ.0 ) THEN
                  UPLO = 'L'
               ELSE
                  UPLO = 'U'
               END IF
*
*              4)      Call SSYEV and SSYEVX.
*
               CALL SLACPY( ' ', N, N, A, LDA, V, LDU )
*
               NTEST = NTEST + 1
               CALL SSYEV( 'V', UPLO, N, A, LDU, D1, WORK, NWORK,
     $                     IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'SSYEV(V,' // UPLO // ')',
     $               IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 400
                  END IF
               END IF
*
*              Do tests 13 and 14.
*
               CALL SSYT21( 1, UPLO, N, 0, V, LDU, D1, D2, A, LDU, Z,
     $                      LDU, TAU, WORK, RESULT( NTEST ) )
*
               CALL SLACPY( ' ', N, N, V, LDU, A, LDA )
*
               NTEST = NTEST + 2
               CALL SSYEV( 'N', UPLO, N, A, LDU, D3, WORK, NWORK,
     $                     IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'SSYEV(N,' // UPLO // ')',
     $               IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 400
                  END IF
               END IF
*
*              Do test 15.
*
               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 390 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( D1( J ) ), ABS( D3( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( D1( J )-D3( J ) ) )
  390          CONTINUE
               RESULT( NTEST ) = TEMP2 / MAX( UNFL,
     $                           ULP*MAX( TEMP1, TEMP2 ) )
*
  400          CONTINUE
               CALL SLACPY( ' ', N, N, V, LDU, A, LDA )
*
               NTEST = NTEST + 1
*
               IF( N.GT.0 ) THEN
                  TEMP3 = MAX( ABS( D1( 1 ) ), ABS( D1( N ) ) )
               IF( IL.NE.1 ) THEN
                  VL = D1( IL ) - MAX( HALF*( D1( IL )-D1( IL-1 ) ),
     $                 TEN*ULP*TEMP3, TEN*RTUNFL )
               ELSE IF( N.GT.0 ) THEN
                  VL = D1( 1 ) - MAX( HALF*( D1( N )-D1( 1 ) ),
     $                 TEN*ULP*TEMP3, TEN*RTUNFL )
               END IF
               IF( IU.NE.N ) THEN
                  VU = D1( IU ) + MAX( HALF*( D1( IU+1 )-D1( IU ) ),
     $                 TEN*ULP*TEMP3, TEN*RTUNFL )
               ELSE IF( N.GT.0 ) THEN
                  VU = D1( N ) + MAX( HALF*( D1( N )-D1( 1 ) ),
     $                 TEN*ULP*TEMP3, TEN*RTUNFL )
               END IF
               ELSE
                  TEMP3 = ZERO
                  VL = ZERO
                  VU = ONE
               END IF
*
               CALL SSYEVX( 'V', 'A', UPLO, N, A, LDU, VL, VU, IL, IU,
     $                      ABSTOL, M, WA1, Z, LDU, WORK, NWORK, IWORK,
     $                      IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'SSYEVX(V,A,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 420
                  END IF
               END IF
*
*              Do tests 16 and 17.
*
               CALL SLACPY( ' ', N, N, V, LDU, A, LDA )
*
               CALL SSYT21( 1, UPLO, N, 0, A, LDU, D1, D2, Z, LDU, V,
     $                      LDU, TAU, WORK, RESULT( NTEST ) )
*
               NTEST = NTEST + 2
               CALL SSYEVX( 'N', 'A', UPLO, N, A, LDU, VL, VU, IL, IU,
     $                      ABSTOL, M2, WA2, Z, LDU, WORK, NWORK, IWORK,
     $                      IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'SSYEVX(N,A,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 420
                  END IF
               END IF
*
*              Do test 18.
*
               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 410 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( WA1( J ) ), ABS( WA2( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( WA1( J )-WA2( J ) ) )
  410          CONTINUE
               RESULT( NTEST ) = TEMP2 / MAX( UNFL,
     $                           ULP*MAX( TEMP1, TEMP2 ) )
*
  420          CONTINUE
*
               NTEST = NTEST + 1
               CALL SLACPY( ' ', N, N, V, LDU, A, LDA )
               CALL SSYEVX( 'V', 'I', UPLO, N, A, LDU, VL, VU, IL, IU,
     $                      ABSTOL, M2, WA2, Z, LDU, WORK, NWORK, IWORK,
     $                      IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'SSYEVX(V,I,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 430
                  END IF
               END IF
*
*              Do tests 19 and 20.
*
               CALL SLACPY( ' ', N, N, V, LDU, A, LDA )
*
               CALL SSYT22( 1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU,
     $                      V, LDU, TAU, WORK, RESULT( NTEST ) )
*
               NTEST = NTEST + 2
               CALL SLACPY( ' ', N, N, V, LDU, A, LDA )
               CALL SSYEVX( 'N', 'I', UPLO, N, A, LDU, VL, VU, IL, IU,
     $                      ABSTOL, M3, WA3, Z, LDU, WORK, NWORK, IWORK,
     $                      IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'SSYEVX(N,I,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 430
                  END IF
               END IF
*
*              Do test 21.
*
               TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
               TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
               RESULT( NTEST ) = ( TEMP1+TEMP2 ) /
     $                           MAX( UNFL, ULP*TEMP3 )
  430          CONTINUE
*
               NTEST = NTEST + 1
               CALL SLACPY( ' ', N, N, V, LDU, A, LDA )
               CALL SSYEVX( 'V', 'V', UPLO, N, A, LDU, VL, VU, IL, IU,
     $                      ABSTOL, M2, WA2, Z, LDU, WORK, NWORK, IWORK,
     $                      IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'SSYEVX(V,V,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 440
                  END IF
               END IF
*
*              Do tests 22 and 23.
*
               CALL SLACPY( ' ', N, N, V, LDU, A, LDA )
*
               CALL SSYT22( 1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU,
     $                      V, LDU, TAU, WORK, RESULT( NTEST ) )
*
               NTEST = NTEST + 2
               CALL SLACPY( ' ', N, N, V, LDU, A, LDA )
               CALL SSYEVX( 'N', 'V', UPLO, N, A, LDU, VL, VU, IL, IU,
     $                      ABSTOL, M3, WA3, Z, LDU, WORK, NWORK, IWORK,
     $                      IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'SSYEVX(N,V,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 440
                  END IF
               END IF
*
               IF( M3.EQ.0 .AND. N.GT.0 ) THEN
                  RESULT( NTEST ) = ULPINV
                  GO TO 440
               END IF
*
*              Do test 24.
*
               TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
               TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
               IF( N.GT.0 ) THEN
                  TEMP3 = MAX( ABS( WA1( 1 ) ), ABS( WA1( N ) ) )
               ELSE
                  TEMP3 = ZERO
               END IF
               RESULT( NTEST ) = ( TEMP1+TEMP2 ) /
     $                           MAX( UNFL, TEMP3*ULP )
*
  440          CONTINUE
*
*              5)      Call SSPEV and SSPEVX.
*
               CALL SLACPY( ' ', N, N, V, LDU, A, LDA )
*
*              Load array WORK with the upper or lower triangular
*              part of the matrix in packed form.
*
               IF( IUPLO.EQ.1 ) THEN
                  INDX = 1
                  DO 460 J = 1, N
                     DO 450 I = 1, J
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  450                CONTINUE
  460             CONTINUE
               ELSE
                  INDX = 1
                  DO 480 J = 1, N
                     DO 470 I = J, N
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  470                CONTINUE
  480             CONTINUE
               END IF
*
               NTEST = NTEST + 1
               CALL SSPEV( 'V', UPLO, N, WORK, D1, Z, LDU, V, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'SSPEV(V,' // UPLO // ')',
     $               IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 540
                  END IF
               END IF
*
*              Do tests 25 and 26.
*
               CALL SSYT21( 1, UPLO, N, 0, A, LDA, D1, D2, Z, LDU, V,
     $                      LDU, TAU, WORK, RESULT( NTEST ) )
*
               IF( IUPLO.EQ.1 ) THEN
                  INDX = 1
                  DO 500 J = 1, N
                     DO 490 I = 1, J
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  490                CONTINUE
  500             CONTINUE
               ELSE
                  INDX = 1
                  DO 520 J = 1, N
                     DO 510 I = J, N
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  510                CONTINUE
  520             CONTINUE
               END IF
*
               NTEST = NTEST + 2
               CALL SSPEV( 'N', UPLO, N, WORK, D3, Z, LDU, V, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'SSPEV(N,' // UPLO // ')',
     $               IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 540
                  END IF
               END IF
*
*              Do test 27.
*
               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 530 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( D1( J ) ), ABS( D3( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( D1( J )-D3( J ) ) )
  530          CONTINUE
               RESULT( NTEST ) = TEMP2 / MAX( UNFL,
     $                           ULP*MAX( TEMP1, TEMP2 ) )
*
*              Load array WORK with the upper or lower triangular part
*              of the matrix in packed form.
*
  540          CONTINUE
               IF( IUPLO.EQ.1 ) THEN
                  INDX = 1
                  DO 560 J = 1, N
                     DO 550 I = 1, J
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  550                CONTINUE
  560             CONTINUE
               ELSE
                  INDX = 1
                  DO 580 J = 1, N
                     DO 570 I = J, N
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  570                CONTINUE
  580             CONTINUE
               END IF
*
               NTEST = NTEST + 1
*
               IF( N.GT.0 ) THEN
                  TEMP3 = MAX( ABS( D1( 1 ) ), ABS( D1( N ) ) )
               IF( IL.NE.1 ) THEN
                  VL = D1( IL ) - MAX( HALF*( D1( IL )-D1( IL-1 ) ),
     $                 TEN*ULP*TEMP3, TEN*RTUNFL )
               ELSE IF( N.GT.0 ) THEN
                  VL = D1( 1 ) - MAX( HALF*( D1( N )-D1( 1 ) ),
     $                 TEN*ULP*TEMP3, TEN*RTUNFL )
               END IF
               IF( IU.NE.N ) THEN
                  VU = D1( IU ) + MAX( HALF*( D1( IU+1 )-D1( IU ) ),
     $                 TEN*ULP*TEMP3, TEN*RTUNFL )
               ELSE IF( N.GT.0 ) THEN
                  VU = D1( N ) + MAX( HALF*( D1( N )-D1( 1 ) ),
     $                 TEN*ULP*TEMP3, TEN*RTUNFL )
               END IF
               ELSE
                  TEMP3 = ZERO
                  VL = ZERO
                  VU = ONE
               END IF
*
               CALL SSPEVX( 'V', 'A', UPLO, N, WORK, VL, VU, IL, IU,
     $                      ABSTOL, M, WA1, Z, LDU, V, IWORK,
     $                      IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'SSPEVX(V,A,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 640
                  END IF
               END IF
*
*              Do tests 28 and 29.
*
               CALL SSYT21( 1, UPLO, N, 0, A, LDU, WA1, D2, Z, LDU, V,
     $                      LDU, TAU, WORK, RESULT( NTEST ) )
*
               NTEST = NTEST + 2
*
               IF( IUPLO.EQ.1 ) THEN
                  INDX = 1
                  DO 600 J = 1, N
                     DO 590 I = 1, J
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  590                CONTINUE
  600             CONTINUE
               ELSE
                  INDX = 1
                  DO 620 J = 1, N
                     DO 610 I = J, N
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  610                CONTINUE
  620             CONTINUE
               END IF
*
               CALL SSPEVX( 'N', 'A', UPLO, N, WORK, VL, VU, IL, IU,
     $                      ABSTOL, M2, WA2, Z, LDU, V, IWORK,
     $                      IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'SSPEVX(N,A,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 640
                  END IF
               END IF
*
*              Do test 30.
*
               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 630 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( WA1( J ) ), ABS( WA2( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( WA1( J )-WA2( J ) ) )
  630          CONTINUE
               RESULT( NTEST ) = TEMP2 / MAX( UNFL,
     $                           ULP*MAX( TEMP1, TEMP2 ) )
*
  640          CONTINUE
               IF( IUPLO.EQ.1 ) THEN
                  INDX = 1
                  DO 660 J = 1, N
                     DO 650 I = 1, J
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  650                CONTINUE
  660             CONTINUE
               ELSE
                  INDX = 1
                  DO 680 J = 1, N
                     DO 670 I = J, N
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  670                CONTINUE
  680             CONTINUE
               END IF
*
               NTEST = NTEST + 1
*
               CALL SSPEVX( 'V', 'I', UPLO, N, WORK, VL, VU, IL, IU,
     $                      ABSTOL, M2, WA2, Z, LDU, V, IWORK,
     $                      IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'SSPEVX(V,I,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 730
                  END IF
               END IF
*
*              Do tests 31 and 32.
*
               CALL SSYT22( 1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU,
     $                      V, LDU, TAU, WORK, RESULT( NTEST ) )
*
               NTEST = NTEST + 2
*
               IF( IUPLO.EQ.1 ) THEN
                  INDX = 1
                  DO 700 J = 1, N
                     DO 690 I = 1, J
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  690                CONTINUE
  700             CONTINUE
               ELSE
                  INDX = 1
                  DO 720 J = 1, N
                     DO 710 I = J, N
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  710                CONTINUE
  720             CONTINUE
               END IF
*
               CALL SSPEVX( 'N', 'I', UPLO, N, WORK, VL, VU, IL, IU,
     $                      ABSTOL, M3, WA3, Z, LDU, V, IWORK,
     $                      IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'SSPEVX(N,I,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 730
                  END IF
               END IF
*
               IF( M3.EQ.0 .AND. N.GT.0 ) THEN
                  RESULT( NTEST ) = ULPINV
                  GO TO 730
               END IF
*
*              Do test 33.
*
               TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
               TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
               IF( N.GT.0 ) THEN
                  TEMP3 = MAX( ABS( WA1( 1 ) ), ABS( WA1( N ) ) )
               ELSE
                  TEMP3 = ZERO
               END IF
               RESULT( NTEST ) = ( TEMP1+TEMP2 ) /
     $                           MAX( UNFL, TEMP3*ULP )
*
  730          CONTINUE
               IF( IUPLO.EQ.1 ) THEN
                  INDX = 1
                  DO 750 J = 1, N
                     DO 740 I = 1, J
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  740                CONTINUE
  750             CONTINUE
               ELSE
                  INDX = 1
                  DO 770 J = 1, N
                     DO 760 I = J, N
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  760                CONTINUE
  770             CONTINUE
               END IF
*
               NTEST = NTEST + 1
*
               CALL SSPEVX( 'V', 'V', UPLO, N, WORK, VL, VU, IL, IU,
     $                      ABSTOL, M2, WA2, Z, LDU, V, IWORK,
     $                      IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'SSPEVX(V,V,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 820
                  END IF
               END IF
*
*              Do tests 34 and 35.
*
               CALL SSYT22( 1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU,
     $                      V, LDU, TAU, WORK, RESULT( NTEST ) )
*
               NTEST = NTEST + 2
*
               IF( IUPLO.EQ.1 ) THEN
                  INDX = 1
                  DO 790 J = 1, N
                     DO 780 I = 1, J
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  780                CONTINUE
  790             CONTINUE
               ELSE
                  INDX = 1
                  DO 810 J = 1, N
                     DO 800 I = J, N
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  800                CONTINUE
  810             CONTINUE
               END IF
*
               CALL SSPEVX( 'N', 'V', UPLO, N, WORK, VL, VU, IL, IU,
     $                      ABSTOL, M3, WA3, Z, LDU, V, IWORK,
     $                      IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'SSPEVX(N,V,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 820
                  END IF
               END IF
*
               IF( M3.EQ.0 .AND. N.GT.0 ) THEN
                  RESULT( NTEST ) = ULPINV
                  GO TO 820
               END IF
*
*              Do test 36.
*
               TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
               TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
               IF( N.GT.0 ) THEN
                  TEMP3 = MAX( ABS( WA1( 1 ) ), ABS( WA1( N ) ) )
               ELSE
                  TEMP3 = ZERO
               END IF
               RESULT( NTEST ) = ( TEMP1+TEMP2 ) /
     $                           MAX( UNFL, TEMP3*ULP )
*
  820          CONTINUE
*
*              6)      Call SSBEV and SSBEVX.
*
               IF( JTYPE.LE.7 ) THEN
                  KD = 0
               ELSE IF( JTYPE.GE.8 .AND. JTYPE.LE.15 ) THEN
                  KD = MAX( N-1, 0 )
               ELSE
                  KD = IHBW
               END IF
*
*              Load array V with the upper or lower triangular part
*              of the matrix in band form.
*
               IF( IUPLO.EQ.1 ) THEN
                  DO 840 J = 1, N
                     DO 830 I = MAX( 1, J-KD ), J
                        V( KD+1+I-J, J ) = A( I, J )
  830                CONTINUE
  840             CONTINUE
               ELSE
                  DO 860 J = 1, N
                     DO 850 I = J, MIN( N, J+KD )
                        V( 1+I-J, J ) = A( I, J )
  850                CONTINUE
  860             CONTINUE
               END IF
*
               NTEST = NTEST + 1
               CALL SSBEV( 'V', UPLO, N, KD, V, LDU, D1, Z, LDU, WORK,
     $                     IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'SSBEV(V,' // UPLO // ')',
     $               IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 920
                  END IF
               END IF
*
*              Do tests 37 and 38.
*
               CALL SSYT21( 1, UPLO, N, 0, A, LDA, D1, D2, Z, LDU, V,
     $                      LDU, TAU, WORK, RESULT( NTEST ) )
*
               IF( IUPLO.EQ.1 ) THEN
                  DO 880 J = 1, N
                     DO 870 I = MAX( 1, J-KD ), J
                        V( KD+1+I-J, J ) = A( I, J )
  870                CONTINUE
  880             CONTINUE
               ELSE
                  DO 900 J = 1, N
                     DO 890 I = J, MIN( N, J+KD )
                        V( 1+I-J, J ) = A( I, J )
  890                CONTINUE
  900             CONTINUE
               END IF
*
               NTEST = NTEST + 2
               CALL SSBEV( 'N', UPLO, N, KD, V, LDU, D3, Z, LDU, WORK,
     $                     IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'SSBEV(N,' // UPLO // ')',
     $               IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 920
                  END IF
               END IF
*
*              Do test 39.
*
               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 910 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( D1( J ) ), ABS( D3( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( D1( J )-D3( J ) ) )
  910          CONTINUE
               RESULT( NTEST ) = TEMP2 / MAX( UNFL,
     $                           ULP*MAX( TEMP1, TEMP2 ) )
*
*              Load array V with the upper or lower triangular part
*              of the matrix in band form.
*
  920          CONTINUE
               IF( IUPLO.EQ.1 ) THEN
                  DO 940 J = 1, N
                     DO 930 I = MAX( 1, J-KD ), J
                        V( KD+1+I-J, J ) = A( I, J )
  930                CONTINUE
  940             CONTINUE
               ELSE
                  DO 960 J = 1, N
                     DO 950 I = J, MIN( N, J+KD )
                        V( 1+I-J, J ) = A( I, J )
  950                CONTINUE
  960             CONTINUE
               END IF
*
               NTEST = NTEST + 1
               CALL SSBEVX( 'V', 'A', UPLO, N, KD, V, LDU, U, LDU, VL,
     $                      VU, IL, IU, ABSTOL, M, WA2, Z, LDU, WORK,
     $                      IWORK, IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'SSBEVX(V,A,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 1020
                  END IF
               END IF
*
*              Do tests 40 and 41.
*
               CALL SSYT21( 1, UPLO, N, 0, A, LDU, WA2, D2, Z, LDU, V,
     $                      LDU, TAU, WORK, RESULT( NTEST ) )
*
               NTEST = NTEST + 2
*
               IF( IUPLO.EQ.1 ) THEN
                  DO 980 J = 1, N
                     DO 970 I = MAX( 1, J-KD ), J
                        V( KD+1+I-J, J ) = A( I, J )
  970                CONTINUE
  980             CONTINUE
               ELSE
                  DO 1000 J = 1, N
                     DO 990 I = J, MIN( N, J+KD )
                        V( 1+I-J, J ) = A( I, J )
  990                CONTINUE
 1000             CONTINUE
               END IF
*
               CALL SSBEVX( 'N', 'A', UPLO, N, KD, V, LDU, U, LDU, VL,
     $                      VU, IL, IU, ABSTOL, M3, WA3, Z, LDU, WORK,
     $                      IWORK, IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'SSBEVX(N,A,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 1020
                  END IF
               END IF
*
*              Do test 42.
*
               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 1010 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( WA2( J ) ), ABS( WA3( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( WA2( J )-WA3( J ) ) )
 1010          CONTINUE
               RESULT( NTEST ) = TEMP2 / MAX( UNFL,
     $                           ULP*MAX( TEMP1, TEMP2 ) )
*
 1020          CONTINUE
               NTEST = NTEST + 1
               IF( IUPLO.EQ.1 ) THEN
                  DO 1040 J = 1, N
                     DO 1030 I = MAX( 1, J-KD ), J
                        V( KD+1+I-J, J ) = A( I, J )
 1030                CONTINUE
 1040             CONTINUE
               ELSE
                  DO 1060 J = 1, N
                     DO 1050 I = J, MIN( N, J+KD )
                        V( 1+I-J, J ) = A( I, J )
 1050                CONTINUE
 1060             CONTINUE
               END IF
*
               CALL SSBEVX( 'V', 'I', UPLO, N, KD, V, LDU, U, LDU, VL,
     $                      VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, WORK,
     $                      IWORK, IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'SSBEVX(V,I,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 1110
                  END IF
               END IF
*
*              Do tests 43 and 44.
*
               CALL SSYT22( 1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU,
     $                      V, LDU, TAU, WORK, RESULT( NTEST ) )
*
               NTEST = NTEST + 2
*
               IF( IUPLO.EQ.1 ) THEN
                  DO 1080 J = 1, N
                     DO 1070 I = MAX( 1, J-KD ), J
                        V( KD+1+I-J, J ) = A( I, J )
 1070                CONTINUE
 1080             CONTINUE
               ELSE
                  DO 1100 J = 1, N
                     DO 1090 I = J, MIN( N, J+KD )
                        V( 1+I-J, J ) = A( I, J )
 1090                CONTINUE
 1100             CONTINUE
               END IF
*
               CALL SSBEVX( 'N', 'I', UPLO, N, KD, V, LDU, U, LDU, VL,
     $                      VU, IL, IU, ABSTOL, M3, WA3, Z, LDU, WORK,
     $                      IWORK, IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'SSBEVX(N,I,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 1110
                  END IF
               END IF
*
*              Do test 45.
*
               TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
               TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
               IF( N.GT.0 ) THEN
                  TEMP3 = MAX( ABS( WA1( 1 ) ), ABS( WA1( N ) ) )
               ELSE
                  TEMP3 = ZERO
               END IF
               RESULT( NTEST ) = ( TEMP1+TEMP2 ) /
     $                           MAX( UNFL, TEMP3*ULP )
*
 1110          CONTINUE
               NTEST = NTEST + 1
               IF( IUPLO.EQ.1 ) THEN
                  DO 1130 J = 1, N
                     DO 1120 I = MAX( 1, J-KD ), J
                        V( KD+1+I-J, J ) = A( I, J )
 1120                CONTINUE
 1130             CONTINUE
               ELSE
                  DO 1150 J = 1, N
                     DO 1140 I = J, MIN( N, J+KD )
                        V( 1+I-J, J ) = A( I, J )
 1140                CONTINUE
 1150             CONTINUE
               END IF
*
               CALL SSBEVX( 'V', 'V', UPLO, N, KD, V, LDU, U, LDU, VL,
     $                      VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, WORK,
     $                      IWORK, IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'SSBEVX(V,V,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 1200
                  END IF
               END IF
*
*              Do tests 46 and 47.
*
               CALL SSYT22( 1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU,
     $                      V, LDU, TAU, WORK, RESULT( NTEST ) )
*
               NTEST = NTEST + 2
*
               IF( IUPLO.EQ.1 ) THEN
                  DO 1170 J = 1, N
                     DO 1160 I = MAX( 1, J-KD ), J
                        V( KD+1+I-J, J ) = A( I, J )
 1160                CONTINUE
 1170             CONTINUE
               ELSE
                  DO 1190 J = 1, N
                     DO 1180 I = J, MIN( N, J+KD )
                        V( 1+I-J, J ) = A( I, J )
 1180                CONTINUE
 1190             CONTINUE
               END IF
*
               CALL SSBEVX( 'N', 'V', UPLO, N, KD, V, LDU, U, LDU, VL,
     $                      VU, IL, IU, ABSTOL, M3, WA3, Z, LDU, WORK,
     $                      IWORK, IWORK( 5*N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'SSBEVX(N,V,' // UPLO //
     $               ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( NTEST ) = ULPINV
                     GO TO 1200
                  END IF
               END IF
*
               IF( M3.EQ.0 .AND. N.GT.0 ) THEN
                  RESULT( NTEST ) = ULPINV
                  GO TO 1200
               END IF
*
*              Do test 48.
*
               TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
               TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
               IF( N.GT.0 ) THEN
                  TEMP3 = MAX( ABS( WA1( 1 ) ), ABS( WA1( N ) ) )
               ELSE
                  TEMP3 = ZERO
               END IF
               RESULT( NTEST ) = ( TEMP1+TEMP2 ) /
     $                           MAX( UNFL, TEMP3*ULP )
*
 1200          CONTINUE
*
 1210       CONTINUE
*
*           End of Loop -- Check for RESULT(j) > THRESH
*
            NTESTT = NTESTT + NTEST
            CALL SLAFTS( 'SST', N, N, JTYPE, NTEST, RESULT, IOLDSD,
     $                   THRESH, NOUNIT, NERRS )
*
 1220    CONTINUE
 1230 CONTINUE
*
*     Summary
*
      CALL ALASVM( 'SST', NOUNIT, NERRS, NTESTT, 0 )
*
 9999 FORMAT( ' SDRVST: ', A, ' returned INFO=', I6, '.', / 9X, 'N=',
     $      I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' )
*
      RETURN
*
*     End of SDRVST
*
      END
