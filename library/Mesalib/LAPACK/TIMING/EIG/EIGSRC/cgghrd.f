      SUBROUTINE CGGHRD( COMPQ, COMPZ, N, ILO, IHI, A, LDA, B, LDB, Q,
     $                   LDQ, Z, LDZ, INFO )
*
*  -- LAPACK routine (instrumented to count operations, version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      CHARACTER          COMPQ, COMPZ
      INTEGER            IHI, ILO, INFO, LDA, LDB, LDQ, LDZ, N
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * ), Q( LDQ, * ),
     $                   Z( LDZ, * )
*     ..
*     ---------------------- Begin Timing Code -------------------------
*     Common block to return operation count and iteration count
*     ITCNT is initialized to 0, OPS is only incremented
*     OPST is used to accumulate small contributions to OPS
*     to avoid roundoff error
*     .. Common blocks ..
      COMMON             / LATIME / OPS, ITCNT
*     ..
*     .. Scalars in Common ..
      REAL               ITCNT, OPS
*     ..
*     ----------------------- End Timing Code --------------------------
*
*
*  Purpose
*  =======
*
*  CGGHRD reduces a pair of complex matrices (A,B) to generalized upper
*  Hessenberg form using unitary similarity transformations, where A is
*  a (generally non-symmetric) square matrix and B is upper
*  triangular.  More precisely, CGGHRD simultaneously decomposes  A
*  into  Q H Z* and  B  into  Q T Z* , where H is upper Hessenberg, T
*  is upper triangular, Q and Z are unitary, and * means conjugate
*  transpose.
*
*  If COMPQ and COMPZ are 'V' or 'I', then the unitary transformations
*  used to reduce (A,B) are accumulated into the arrays Q and Z s.t.:
*
*       Q(in) A(in) Z(in)* = Q(out) A(out) Z(out)*
*       Q(in) B(in) Z(in)* = Q(out) B(out) Z(out)*
*
*  CGGHRD implements an unblocked form of the method, there being no
*  blocked form known at this time.
*
*  Arguments
*  =========
*
*  COMPQ   (input) CHARACTER*1
*          = 'N': do not compute Q.
*          = 'V': accumulate the row transformations into Q.
*          = 'I': overwrite the array Q with the row transformations.
*
*  COMPZ   (input) CHARACTER*1
*          = 'N': do not compute Z.
*          = 'V': accumulate the column transformations into Z.
*          = 'I': overwrite the array Z with the column transformations.
*
*  N       (input) INTEGER
*          The number of rows and columns in the matrices A, B, Q, and
*          Z.  N >= 0.
*
*  ILO     (input) INTEGER
*          Columns 1 through ILO-1 of A are assumed to be in upper
*          triangular form already, and will not be modified.  ILO must
*          be at least 1.
*
*  IHI     (input) INTEGER
*          Rows IHI+1 through N of A are assumed to be in upper
*          triangular form already, and will not be touched.  IHI may
*          not be greater than N.
*
*  A       (input/output) COMPLEX array, dimension (LDA, N)
*          On entry, the N x N general matrix to be reduced.
*          On exit, the upper triangle and the first subdiagonal of A
*          are overwritten with the Hessenberg matrix H, and the rest
*          is set to zero.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max( 1, N ).
*
*  B       (input/output) COMPLEX array, dimension (LDB, N)
*          On entry, the N x N upper triangular matrix B.
*          On exit, the upper triangular matrix T = Q* B Z.  The
*          entries below the diagonal are set to zero.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max( 1, N ).
*
*  Q       (input/output) COMPLEX array, dimension (LDQ, N)
*          If COMPQ='N', then Q will not be referenced.
*          If COMPQ='V', then the conjugate transpose of the Givens
*             transformations which are applied to A and B on the left
*             will be applied to the array Q on the right.
*          If COMPQ='I', then Q will first be overwritten with an
*             identity matrix, and then the Givens transformations will
*             be applied to Q as in the case COMPQ='V'.
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q.  LDQ must be at least
*          1.  If COMPQ='V' or 'I', then LDQ must also be at least N.
*
*  Z       (input/output) COMPLEX array, dimension (LDZ, N)
*          If COMPZ='N', then Z will not be referenced.
*          If COMPZ='V', then the Givens transformations which are
*             applied to A and B on the right will be applied to the
*             array Z on the right.
*          If COMPZ='I', then Z will first be overwritten with an
*             identity matrix, and then the Givens transformations will
*             be applied to Z as in the case COMPZ='V'.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ must be at least
*          1.  If COMPZ='V' or 'I', then LDZ must also be at least N.
*
*  INFO    (output) INTEGER
*          = 0:  normal return.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  This routine reduces A to Hessenberg and B to triangular form by
*  an unblocked reduction, as described in _Matrix_Computations_,
*  by Golub and van Loan.
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX            CONE, CZERO
      PARAMETER          ( CONE = 1.0E+0, CZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ILQ, ILZ
      INTEGER            ICOMPQ, ICOMPZ, JCOL, JROW
      REAL               C, TEMP
      COMPLEX            CTEMP, S
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLARTG, CLAZRO, CROT, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX, REAL
*     ..
*     .. Executable Statements ..
*
*     Decode COMPQ
*
      IF( LSAME( COMPQ, 'N' ) ) THEN
         ILQ = .FALSE.
         ICOMPQ = 1
      ELSE IF( LSAME( COMPQ, 'V' ) ) THEN
         ILQ = .TRUE.
         ICOMPQ = 2
      ELSE IF( LSAME( COMPQ, 'I' ) ) THEN
         ILQ = .TRUE.
         ICOMPQ = 3
      ELSE
         ICOMPQ = 0
      END IF
*
*     Decode COMPZ
*
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ILZ = .FALSE.
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ILZ = .TRUE.
         ICOMPZ = 2
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ILZ = .TRUE.
         ICOMPZ = 3
      ELSE
         ICOMPZ = 0
      END IF
*
*     Test the input parameters.
*
      INFO = 0
      IF( ICOMPQ.LE.0 ) THEN
         INFO = -1
      ELSE IF( ICOMPZ.LE.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( ILO.LT.1 ) THEN
         INFO = -4
      ELSE IF( IHI.GT.N .OR. IHI.LT.ILO-1 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( ( ILQ .AND. LDQ.LT.N ) .OR. LDQ.LT.1 ) THEN
         INFO = -11
      ELSE IF( ( ILZ .AND. LDZ.LT.N ) .OR. LDZ.LT.1 ) THEN
         INFO = -13
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CGGHRD', -INFO )
         RETURN
      END IF
*
*     Initialize Q and Z if desired.
*
      IF( ICOMPQ.EQ.3 )
     $   CALL CLAZRO( N, N, CZERO, CONE, Q, LDQ )
      IF( ICOMPZ.EQ.3 )
     $   CALL CLAZRO( N, N, CZERO, CONE, Z, LDZ )
*
*     Quick return if possible
*
      IF( N.LE.1 )
     $   RETURN
*
*     Zero out lower triangle of B
*
      DO 20 JCOL = 1, N - 1
         DO 10 JROW = JCOL + 1, N
            B( JROW, JCOL ) = CZERO
   10    CONTINUE
   20 CONTINUE
*
*     Reduce A and B
*
      DO 40 JCOL = ILO, IHI - 2
*
         DO 30 JROW = IHI, JCOL + 2, -1
*
*           Step 1: rotate rows JROW-1, JROW to kill A(JROW,JCOL)
*
            CTEMP = A( JROW-1, JCOL )
            CALL CLARTG( CTEMP, A( JROW, JCOL ), C, S,
     $                   A( JROW-1, JCOL ) )
            A( JROW, JCOL ) = CZERO
            CALL CROT( N-JCOL, A( JROW-1, JCOL+1 ), LDA,
     $                 A( JROW, JCOL+1 ), LDA, C, S )
            CALL CROT( N+2-JROW, B( JROW-1, JROW-1 ), LDB,
     $                 B( JROW, JROW-1 ), LDB, C, S )
            IF( ILQ )
     $         CALL CROT( N, Q( 1, JROW-1 ), 1, Q( 1, JROW ), 1, C,
     $                    CONJG( S ) )
*
*           Step 2: rotate columns JROW, JROW-1 to kill B(JROW,JROW-1)
*
            CTEMP = B( JROW, JROW )
            CALL CLARTG( CTEMP, B( JROW, JROW-1 ), C, S,
     $                   B( JROW, JROW ) )
            B( JROW, JROW-1 ) = CZERO
            CALL CROT( IHI, A( 1, JROW ), 1, A( 1, JROW-1 ), 1, C, S )
            CALL CROT( JROW-1, B( 1, JROW ), 1, B( 1, JROW-1 ), 1, C,
     $                 S )
            IF( ILZ )
     $         CALL CROT( N, Z( 1, JROW ), 1, Z( 1, JROW-1 ), 1, C, S )
   30    CONTINUE
   40 CONTINUE
*
*     ---------------------- Begin Timing Code -------------------------
*     Operation count:                                          factor
*     * number of calls to CLARTG   TEMP                          *32
*     * total number of rows/cols
*       rotated in A and B          TEMP*[6n + 2(ihi-ilo) + 5]/6  *20
*     * rows rotated in Q           TEMP*n/2                      *20
*     * rows rotated in Z           TEMP*n/2                      *20
*
      TEMP = REAL( IHI-ILO )*REAL( IHI-ILO-1 )
      JROW = 20*N + 7*( IHI-ILO ) + 27 + 32
      IF( ILQ )
     $   JROW = JROW + 10*N
      IF( ILZ )
     $   JROW = JROW + 10*N
      OPS = OPS + ( REAL( JROW )-REAL( IHI-ILO+1 ) / REAL( 3 ) )*TEMP
      ITCNT = 0
*
*     ----------------------- End Timing Code --------------------------
*
      RETURN
*
*     End of CGGHRD
*
      END
