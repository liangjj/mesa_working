      SUBROUTINE ZTGEVC( JOB, SIDE, SELECT, N, A, LDA, B, LDB, VL, LDVL,
     $                   VR, LDVR, MM, M, WORK, RWORK, INFO )
*
*  -- LAPACK routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      CHARACTER          JOB, SIDE
      INTEGER            INFO, LDA, LDB, LDVL, LDVR, M, MM, N
*     ..
*     .. Array Arguments ..
      LOGICAL            SELECT( * )
      DOUBLE PRECISION   RWORK( N, * )
      COMPLEX*16         A( LDA, * ), B( LDB, * ), VL( LDVL, * ),
     $                   VR( LDVR, * ), WORK( N, * )
*     ..
*
*
*  Purpose
*  =======
*
*  ZTGEVC computes selected left and/or right generalized eigenvectors
*  of a pair of complex upper triangular matrices (A,B).  The
*  j-th generalized left and right eigenvectors are  y  and  x, resp.,
*  such that:
*          H
*         y  (A - wB) = 0       (A - wB)x = 0
*
*                                                                      H
*  Note: the left eigenvector is sometimes defined as the row vector  y
*        but ZTGEVC computes the column vector y.
*
*  Arguments
*  =========
*
*  JOB     (input) CHARACTER*1
*          = 'A': compute All (left/right/left+right) generalized
*                 eigenvectors of (A,B);
*          = 'S': compute Selected (left/right/left+right) generalized
*                 eigenvectors of (A,B) -- see the description of the
*                 argument SELECT;
*          = 'B' or 'T': compute all (left/right/left+right) generalized
*                 eigenvectors of (A,B), and Back Transform them
*                 using the initial contents of VL/VR -- see the
*                 descriptions of the arguments VL and VR.
*
*  SIDE    (input) CHARACTER*1
*          Specifies for which side eigenvectors are to be computed:
*          = 'R': compute right eigenvectors only;
*          = 'L': compute left eigenvectors only;
*          = 'B': compute both right and left eigenvectors.
*
*  SELECT  (input) LOGICAL array, dimension (N)
*          If JOB='S', then SELECT specifies the (generalized)
*          eigenvectors to be computed.  To get the eigenvector
*          corresponding to the j-th eigenvalue, set SELECT(j) to
*          .TRUE..  If JOB='A', 'B', or 'T', SELECT is not referenced,
*          and all eigenvectors are selected.
*
*  N       (input) INTEGER
*          The order of the matrices A and B.  N >= 0.
*
*  A       (input) COMPLEX*16 array, dimension (LDA,N)
*          One of the pair of matrices whose generalized eigenvectors
*          are to be computed.  It must be upper triangular.
*
*  LDA     (input) INTEGER
*          The leading dimension of array A.  LDA >= max(1, N).
*
*  B       (input) COMPLEX*16 array, dimension (LDB,N)
*          The other of the pair of matrices whose generalized
*          eigenvectors are to be computed.  It must be upper
*          triangular with real diagonal elements.
*
*  LDB     (input) INTEGER
*          The leading dimension of array B.  LDB >= max(1, N).
*
*  VL      (input/output) COMPLEX*16 array, dimension (LDVL,MM)
*          The left eigenvectors (column vectors -- see the note in
*          "Purpose".)
*          If JOB='A', then all left eigenvectors of (A,B) will be
*             computed and stored in VL.
*          If JOB='S', then only the eigenvectors selected by SELECT
*             will be computed; the first selected eigenvector will go
*             in column 1, the second in column 2, etc.
*          If JOB='B' or 'T', then all left eigenvectors of (A,B)
*             will be computed and multiplied (on the left) by the
*             matrix found in VL on entry to ZTGEVC.  Usually, this
*             will be the Q matrix computed by ZGGHRD and ZHGEQZ,
*             so that on exit, VL will contain the left eigenvectors
*             of the original matrix pair.
*          In any case, each eigenvector will be scaled so the largest
*          component of each vector has
*          abs(real part) + abs(imag. part)=1, *unless*  A and B both
*          have a zero in the diagonal entry corresponding to the
*          eigenvector, in which case the eigenvector will be zero.
*          If SIDE = 'R', VL is not referenced.
*
*  LDVL    (input) INTEGER
*          The leading dimension of array VL.  LDVL >= 1; if SIDE = 'B'
*          or 'L', LDVL >= N.
*
*  VR      (input/output) COMPLEX*16 array, dimension (LDVR,MM)
*          The right eigenvectors.
*          If JOB='A', then all right eigenvectors of (A,B) will be
*             computed and stored in VR.
*          If JOB='S', then only the eigenvectors selected by SELECT
*             will be computed; the first selected eigenvector will go
*             in column 1, the second in column 2, etc.
*          If JOB='B' or 'T', then all right eigenvectors of (A,B)
*             will be computed and multiplied (on the left) by the
*             matrix found in VR on entry to ZTGEVC.  Usually, this
*             will be the Z matrix computed by ZGGHRD and ZHGEQZ,
*             so that on exit, VR will contain the right eigenvectors
*             of the original matrix pair.
*          In any case, each eigenvector will be scaled so the largest
*          component of each vector has
*          abs(real part) + abs(imag. part)=1, *unless*  A and B both
*          have a zero in the diagonal entry corresponding to the
*          eigenvector, in which case the eigenvector will be zero.
*          If SIDE = 'L', VR is not referenced.
*
*  LDVR    (input) INTEGER
*          The leading dimension of array VR.  LDVR >= 1; if SIDE = 'B'
*          or 'R', LDVR >= N.
*
*  MM      (input) INTEGER
*          The number of columns in VL and/or VR.
*          If JOB='A', 'B', or 'T', then MM >= N.
*          If JOB='S', then MM must be at least the number of .TRUE.
*             values in SELECT(1:N).
*
*  M       (output) INTEGER
*          The number of columns in VL and/or VR actually
*          used to store the eigenvectors.
*
*  WORK    (workspace) COMPLEX*16 array, dimension ( 2, N )
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension ( 2, N )
*
*  INFO    (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = 0.0D+0, CONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            COMPL, COMPR, ILALL, ILBACK, ILBBAD, ILCOMP,
     $                   LSA, LSB
      INTEGER            I, IBEG, IEIG, IEND, IJOB, IM, ISIDE, ISRC, J,
     $                   JE, JR
      DOUBLE PRECISION   ACOEFA, ACOEFF, ANORM, ASCALE, BCOEFA, BIG,
     $                   BIGNUM, BNORM, BSCALE, DMIN, SAFMIN, SBETA,
     $                   SCALE, SMALL, TEMP, ULP, XMAX
      COMPLEX*16         BCOEFF, CA, CB, D, SALPHA, SUM, SUMA, SUMB, X
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLABAD, XERBLA, ZGEMV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCMPLX, DCONJG, DIMAG, MAX, MIN
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   ABS1
*     ..
*     .. Statement Function definitions ..
      ABS1( X ) = ABS( DBLE( X ) ) + ABS( DIMAG( X ) )
*     ..
*     .. Executable Statements ..
*
*     Decode and Test the input parameters
*
      IF( LSAME( JOB, 'A' ) ) THEN
         IJOB = 1
         ILALL = .TRUE.
         ILBACK = .FALSE.
      ELSE IF( LSAME( JOB, 'S' ) ) THEN
         IJOB = 2
         ILALL = .FALSE.
         ILBACK = .FALSE.
      ELSE IF( LSAME( JOB, 'B' ) .OR. LSAME( JOB, 'T' ) ) THEN
         IJOB = 3
         ILALL = .TRUE.
         ILBACK = .TRUE.
      ELSE
         IJOB = -1
      END IF
*
      IF( LSAME( SIDE, 'R' ) ) THEN
         ISIDE = 1
         COMPL = .FALSE.
         COMPR = .TRUE.
      ELSE IF( LSAME( SIDE, 'L' ) ) THEN
         ISIDE = 2
         COMPL = .TRUE.
         COMPR = .FALSE.
      ELSE IF( LSAME( SIDE, 'B' ) ) THEN
         ISIDE = 3
         COMPL = .TRUE.
         COMPR = .TRUE.
      ELSE
         ISIDE = -1
      END IF
*
*     Count the number of eigenvectors
*
      IF( .NOT.ILALL ) THEN
         IM = 0
         DO 10 J = 1, N
            IF( SELECT( J ) )
     $         IM = IM + 1
   10    CONTINUE
      ELSE
         IM = N
      END IF
*
*     Check diagonal of B
*
      ILBBAD = .FALSE.
      DO 20 J = 1, N
         IF( DIMAG( B( J, J ) ).NE.ZERO )
     $      ILBBAD = .TRUE.
   20 CONTINUE
*
      INFO = 0
      IF( IJOB.LT.0 ) THEN
         INFO = -1
      ELSE IF( ISIDE.LT.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( ILBBAD ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( COMPL .AND. LDVL.LT.N .OR. LDVL.LT.1 ) THEN
         INFO = -10
      ELSE IF( COMPR .AND. LDVR.LT.N .OR. LDVR.LT.1 ) THEN
         INFO = -12
      ELSE IF( MM.LT.IM ) THEN
         INFO = -13
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZTGEVC', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      M = IM
      IF( N.EQ.0 )
     $   RETURN
*
*     Machine Constants
*
      SAFMIN = DLAMCH( 'Safe minimum' )
      BIG = ONE / SAFMIN
      CALL DLABAD( SAFMIN, BIG )
      ULP = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' )
      SMALL = SAFMIN*N / ULP
      BIG = ONE / SMALL
      BIGNUM = ONE / ( SAFMIN*N )
*
*     Compute the 1-norm of each column of the strictly upper triangular
*     part of A and B to check for possible overflow in the triangular
*     solver.
*
      ANORM = ABS1( A( 1, 1 ) )
      BNORM = ABS1( B( 1, 1 ) )
      RWORK( 1, 1 ) = ZERO
      RWORK( 1, 2 ) = ZERO
      DO 40 J = 2, N
         RWORK( J, 1 ) = ZERO
         RWORK( J, 2 ) = ZERO
         DO 30 I = 1, J - 1
            RWORK( J, 1 ) = RWORK( J, 1 ) + ABS1( A( I, J ) )
            RWORK( J, 2 ) = RWORK( J, 2 ) + ABS1( B( I, J ) )
   30    CONTINUE
         ANORM = MAX( ANORM, RWORK( J, 1 )+ABS1( A( J, J ) ) )
         BNORM = MAX( BNORM, RWORK( J, 2 )+ABS1( B( J, J ) ) )
   40 CONTINUE
*
      ASCALE = ONE / MAX( ANORM, SAFMIN )
      BSCALE = ONE / MAX( BNORM, SAFMIN )
*
*     Left eigenvectors
*
      IF( COMPL ) THEN
         IEIG = 0
*
*        Main loop over eigenvalues
*
         DO 140 JE = 1, N
            IF( ILALL ) THEN
               ILCOMP = .TRUE.
            ELSE
               ILCOMP = SELECT( JE )
            END IF
            IF( ILCOMP ) THEN
               IEIG = IEIG + 1
*
               IF( ABS1( A( JE, JE ) ).LE.SAFMIN .AND.
     $             ABS( DBLE( B( JE, JE ) ) ).LE.SAFMIN ) THEN
*
*                 Singular matrix pencil -- return zero eigenvector
*
                  DO 50 JR = 1, N
                     VL( JR, IEIG ) = CZERO
   50             CONTINUE
                  GO TO 140
               END IF
*
*              Non-singular eigenvalue:
*              Compute coefficients  a  and  b  in
*                   H
*                 y  ( a A - b B ) = 0
*
               TEMP = ONE / MAX( ABS1( A( JE, JE ) )*ASCALE,
     $                ABS( DBLE( B( JE, JE ) ) )*BSCALE, SAFMIN )
               SALPHA = ( TEMP*A( JE, JE ) )*ASCALE
               SBETA = ( TEMP*DBLE( B( JE, JE ) ) )*BSCALE
               ACOEFF = SBETA*ASCALE
               BCOEFF = SALPHA*BSCALE
*
*              Scale to avoid underflow
*
               LSA = ABS( SBETA ).GE.SAFMIN .AND. ABS( ACOEFF ).LT.SMALL
               LSB = ABS1( SALPHA ).GE.SAFMIN .AND. ABS1( BCOEFF ).LT.
     $               SMALL
*
               SCALE = ONE
               IF( LSA )
     $            SCALE = ( SMALL / ABS( SBETA ) )*MIN( ANORM, BIG )
               IF( LSB )
     $            SCALE = MAX( SCALE, ( SMALL / ABS1( SALPHA ) )*
     $                    MIN( BNORM, BIG ) )
               IF( LSA .OR. LSB ) THEN
                  SCALE = MIN( SCALE, ONE /
     $                    ( SAFMIN*MAX( ONE, ABS( ACOEFF ),
     $                    ABS1( BCOEFF ) ) ) )
                  IF( LSA ) THEN
                     ACOEFF = ASCALE*( SCALE*SBETA )
                  ELSE
                     ACOEFF = SCALE*ACOEFF
                  END IF
                  IF( LSB ) THEN
                     BCOEFF = BSCALE*( SCALE*SALPHA )
                  ELSE
                     BCOEFF = SCALE*BCOEFF
                  END IF
               END IF
*
               ACOEFA = ABS( ACOEFF )
               BCOEFA = ABS1( BCOEFF )
               XMAX = ONE
               DO 60 JR = 1, N
                  WORK( JR, 1 ) = CZERO
   60          CONTINUE
               WORK( JE, 1 ) = CONE
               DMIN = MAX( ULP*ACOEFA*ANORM, ULP*BCOEFA*BNORM, SAFMIN )
*
*                                              H
*              Triangular solve of  (a A - b B)  y = 0
*
*                                      H
*              (rowwise in  (a A - b B) , or columnwise in a A - b B)
*
               DO 100 J = JE + 1, N
*
*                 Compute
*                       j-1
*                 SUM = sum  conjg( a*A(k,j) - b*B(k,j) )*x(k)
*                       k=je
*                 (Scale if necessary)
*
                  TEMP = ONE / XMAX
                  IF( ACOEFA*RWORK( J, 1 )+BCOEFA*RWORK( J, 2 ).GT.
     $                BIGNUM*TEMP ) THEN
                     DO 70 JR = JE, J - 1
                        WORK( JR, 1 ) = TEMP*WORK( JR, 1 )
   70                CONTINUE
                     XMAX = ONE
                  END IF
                  SUMA = CZERO
                  SUMB = CZERO
*
                  DO 80 JR = JE, J - 1
                     SUMA = SUMA + DCONJG( A( JR, J ) )*WORK( JR, 1 )
                     SUMB = SUMB + DCONJG( B( JR, J ) )*WORK( JR, 1 )
   80             CONTINUE
                  SUM = ACOEFF*SUMA - DCONJG( BCOEFF )*SUMB
*
*                 Form x(j) = - SUM / conjg( a*A(j,j) - b*B(j,j) )
*
*                 with scaling and perturbation of the denominator
*
                  D = DCONJG( ACOEFF*A( J, J )-BCOEFF*B( J, J ) )
                  IF( ABS1( D ).LE.DMIN )
     $               D = DCMPLX( DMIN )
*
                  IF( ABS1( D ).LT.ONE ) THEN
                     IF( ABS1( SUM ).GE.BIGNUM*ABS1( D ) ) THEN
                        TEMP = ONE / ABS1( SUM )
                        DO 90 JR = JE, J - 1
                           WORK( JR, 1 ) = TEMP*WORK( JR, 1 )
   90                   CONTINUE
                        XMAX = TEMP*XMAX
                        SUM = TEMP*SUM
                     END IF
                  END IF
                  WORK( J, 1 ) = -SUM / D
                  XMAX = MAX( XMAX, ABS1( WORK( J, 1 ) ) )
  100          CONTINUE
*
*              Back transform eigenvector if JOB='B' or 'T'.
*
               IF( ILBACK ) THEN
                  CALL ZGEMV( 'N', N, N+1-JE, CONE, VL( 1, JE ), LDVL,
     $                        WORK( JE, 1 ), 1, CZERO, WORK( 1, 2 ), 1 )
                  ISRC = 2
                  IBEG = 1
               ELSE
                  ISRC = 1
                  IBEG = JE
               END IF
*
*              Copy and scale eigenvector into column of VL
*
               XMAX = ZERO
               DO 110 JR = IBEG, N
                  XMAX = MAX( XMAX, ABS1( WORK( JR, ISRC ) ) )
  110          CONTINUE
*
               IF( XMAX.GT.SAFMIN ) THEN
                  TEMP = ONE / XMAX
                  DO 120 JR = IBEG, N
                     VL( JR, IEIG ) = TEMP*WORK( JR, ISRC )
  120             CONTINUE
               ELSE
                  IBEG = N + 1
               END IF
*
               DO 130 JR = 1, IBEG - 1
                  VL( JR, IEIG ) = CZERO
  130          CONTINUE
*
            END IF
  140    CONTINUE
      END IF
*
*     Right eigenvectors
*
      IF( COMPR ) THEN
         IEIG = IM + 1
*
*        Main loop over eigenvalues
*
         DO 250 JE = N, 1, -1
            IF( ILALL ) THEN
               ILCOMP = .TRUE.
            ELSE
               ILCOMP = SELECT( JE )
            END IF
            IF( ILCOMP ) THEN
               IEIG = IEIG - 1
*
               IF( ABS1( A( JE, JE ) ).LE.SAFMIN .AND.
     $             ABS( DBLE( B( JE, JE ) ) ).LE.SAFMIN ) THEN
*
*                 Singular matrix pencil -- return zero eigenvector
*
                  DO 150 JR = 1, N
                     VR( JR, IEIG ) = CZERO
  150             CONTINUE
                  GO TO 250
               END IF
*
*              Non-singular eigenvalue:
*              Compute coefficients  a  and  b  in
*
*              ( a A - b B ) x  = 0
*
               TEMP = ONE / MAX( ABS1( A( JE, JE ) )*ASCALE,
     $                ABS( DBLE( B( JE, JE ) ) )*BSCALE, SAFMIN )
               SALPHA = ( TEMP*A( JE, JE ) )*ASCALE
               SBETA = ( TEMP*DBLE( B( JE, JE ) ) )*BSCALE
               ACOEFF = SBETA*ASCALE
               BCOEFF = SALPHA*BSCALE
*
*              Scale to avoid underflow
*
               LSA = ABS( SBETA ).GE.SAFMIN .AND. ABS( ACOEFF ).LT.SMALL
               LSB = ABS1( SALPHA ).GE.SAFMIN .AND. ABS1( BCOEFF ).LT.
     $               SMALL
*
               SCALE = ONE
               IF( LSA )
     $            SCALE = ( SMALL / ABS( SBETA ) )*MIN( ANORM, BIG )
               IF( LSB )
     $            SCALE = MAX( SCALE, ( SMALL / ABS1( SALPHA ) )*
     $                    MIN( BNORM, BIG ) )
               IF( LSA .OR. LSB ) THEN
                  SCALE = MIN( SCALE, ONE /
     $                    ( SAFMIN*MAX( ONE, ABS( ACOEFF ),
     $                    ABS1( BCOEFF ) ) ) )
                  IF( LSA ) THEN
                     ACOEFF = ASCALE*( SCALE*SBETA )
                  ELSE
                     ACOEFF = SCALE*ACOEFF
                  END IF
                  IF( LSB ) THEN
                     BCOEFF = BSCALE*( SCALE*SALPHA )
                  ELSE
                     BCOEFF = SCALE*BCOEFF
                  END IF
               END IF
*
               ACOEFA = ABS( ACOEFF )
               BCOEFA = ABS1( BCOEFF )
               XMAX = ONE
               DO 160 JR = 1, N
                  WORK( JR, 1 ) = CZERO
  160          CONTINUE
               WORK( JE, 1 ) = CONE
               DMIN = MAX( ULP*ACOEFA*ANORM, ULP*BCOEFA*BNORM, SAFMIN )
*
*              Triangular solve of  (a A - b B) x = 0  (columnwise)
*
*              WORK(1:j-1) contains sums w,
*              WORK(j+1:JE) contains x
*
               DO 170 JR = 1, JE - 1
                  WORK( JR, 1 ) = ACOEFF*A( JR, JE ) -
     $                            BCOEFF*B( JR, JE )
  170          CONTINUE
               WORK( JE, 1 ) = CONE
*
               DO 210 J = JE - 1, 1, -1
*
*                 Form x(j) := - w(j) / d
*                 with scaling and perturbation of the denominator
*
                  D = ACOEFF*A( J, J ) - BCOEFF*B( J, J )
                  IF( ABS1( D ).LE.DMIN )
     $               D = DCMPLX( DMIN )
*
                  IF( ABS1( D ).LT.ONE ) THEN
                     IF( ABS1( WORK( J, 1 ) ).GE.BIGNUM*ABS1( D ) ) THEN
                        TEMP = ONE / ABS1( WORK( J, 1 ) )
                        DO 180 JR = 1, JE
                           WORK( JR, 1 ) = TEMP*WORK( JR, 1 )
  180                   CONTINUE
                     END IF
                  END IF
*
                  WORK( J, 1 ) = -WORK( J, 1 ) / D
*
                  IF( J.GT.1 ) THEN
*
*                    w = w + x(j)*(a A(*,j) - b B(*,j) ) with scaling
*
                     IF( ABS1( WORK( J, 1 ) ).GT.ONE ) THEN
                        TEMP = ONE / ABS1( WORK( J, 1 ) )
                        IF( ACOEFA*RWORK( J, 1 )+BCOEFA*
     $                      RWORK( J, 2 ).GE.BIGNUM*TEMP ) THEN
                           DO 190 JR = 1, JE
                              WORK( JR, 1 ) = TEMP*WORK( JR, 1 )
  190                      CONTINUE
                        END IF
                     END IF
*
                     CA = ACOEFF*WORK( J, 1 )
                     CB = BCOEFF*WORK( J, 1 )
                     DO 200 JR = 1, J - 1
                        WORK( JR, 1 ) = WORK( JR, 1 ) + CA*A( JR, J ) -
     $                                  CB*B( JR, J )
  200                CONTINUE
                  END IF
  210          CONTINUE
*
*              Back transform eigenvector if JOB='B' or 'T'.
*
               IF( ILBACK ) THEN
                  CALL ZGEMV( 'N', N, JE, CONE, VR, LDVR, WORK, 1,
     $                        CZERO, WORK( 1, 2 ), 1 )
                  ISRC = 2
                  IEND = N
               ELSE
                  ISRC = 1
                  IEND = JE
               END IF
*
*              Copy and scale eigenvector into column of VR
*
               XMAX = ZERO
               DO 220 JR = 1, IEND
                  XMAX = MAX( XMAX, ABS1( WORK( JR, ISRC ) ) )
  220          CONTINUE
*
               IF( XMAX.GT.SAFMIN ) THEN
                  TEMP = ONE / XMAX
                  DO 230 JR = 1, IEND
                     VR( JR, IEIG ) = TEMP*WORK( JR, ISRC )
  230             CONTINUE
               ELSE
                  IEND = 0
               END IF
*
               DO 240 JR = IEND + 1, N
                  VR( JR, IEIG ) = CZERO
  240          CONTINUE
*
            END IF
  250    CONTINUE
      END IF
*
      RETURN
*
*     End of ZTGEVC
*
      END
