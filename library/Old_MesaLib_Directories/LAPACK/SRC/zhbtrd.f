      SUBROUTINE ZHBTRD( VECT, UPLO, N, KD, AB, LDAB, D, E, Q, LDQ,
     $                   WORK, INFO )
*
*  -- LAPACK routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO, VECT
      INTEGER            INFO, KD, LDAB, LDQ, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * )
      COMPLEX*16         AB( LDAB, * ), Q( LDQ, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZHBTRD reduces a complex Hermitian band matrix A to real symmetric
*  tridiagonal form T by a unitary similarity transformation:
*  Q**H * A * Q = T.
*
*  Arguments
*  =========
*
*  VECT    (input) CHARACTER*1
*          = 'N': do not form Q;
*          = 'V': form Q.
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  KD      (input) INTEGER
*          The number of superdiagonals of the matrix A if UPLO = 'U',
*          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
*
*  AB      (input/output) COMPLEX*16 array, dimension (LDAB,N)
*          On entry, the upper or lower triangle of the Hermitian band
*          matrix A, stored in the first KD+1 rows of the array.  The
*          j-th column of A is stored in the j-th column of the array AB
*          as follows:
*          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
*          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
*          On exit, the diagonal elements of A are overwritten by the
*          diagonal elements of the tridiagonal matrix T; if KD > 0, the
*          elements on the first superdiagonal (if UPLO = 'U') or the
*          first subdiagonal (if UPLO = 'L') are overwritten by the
*          offdiagonal elements of T; the rest of A is overwritten by
*          values generated during the reduction.
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= KD+1.
*
*  D       (output) DOUBLE PRECISION array, dimension (N)
*          The diagonal elements of the tridiagonal matrix T.
*
*  E       (output) DOUBLE PRECISION array, dimension (N-1)
*          The off-diagonal elements of the tridiagonal matrix T:
*          E(i) = T(i,i+1) if UPLO = 'U'; E(i) = T(i+1,i) if UPLO = 'L'.
*
*  Q       (output) COMPLEX*16 array, dimension (LDQ,N)
*          If VECT = 'V', the N-by-N unitary matrix Q.
*          If VECT = 'N', the array Q is not referenced.
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q.
*          LDQ >= max(1,N) if VECT = 'V'.
*
*  WORK    (workspace) COMPLEX*16 array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER, WANTQ
      INTEGER            I, J, J1, J2, K, KD1, KDN, L, NR, NRT
      DOUBLE PRECISION   ABST
      COMPLEX*16         T, TEMP
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLACGV, ZLAR2V, ZLARGV, ZLARTG, ZLARTV,
     $                   ZLAZRO, ZROT, ZSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCONJG, MAX, MIN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      WANTQ = LSAME( VECT, 'V' )
      UPPER = LSAME( UPLO, 'U' )
      KD1 = KD + 1
      INFO = 0
      IF( .NOT.WANTQ .AND. .NOT.LSAME( VECT, 'N' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( KD.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.KD1 ) THEN
         INFO = -6
      ELSE IF( LDQ.LT.MAX( 1, N ) .AND. WANTQ ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZHBTRD', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Initialize Q to the unit matrix, if needed
*
      IF( WANTQ )
     $   CALL ZLAZRO( N, N, ZERO, ONE, Q, LDQ )
*
*     Wherever possible, plane rotations are generated and applied in
*     vector operations of length NR over the index set J1:J2:KD1.
*
*     The real cosines and complex sines of the plane rotations are
*     stored in the arrays D and WORK.
*
      KDN = MIN( N-1, KD )
      IF( UPPER ) THEN
*
         IF( KD.GT.1 ) THEN
*
*           Reduce to complex Hermitian tridiagonal form, working with
*           the upper triangle
*
            NR = 0
            J1 = KDN + 2
            J2 = 1
*
            AB( KD1, 1 ) = DBLE( AB( KD1, 1 ) )
            DO 60 I = 1, N - 2
*
*              Reduce i-th row of matrix to tridiagonal form
*
               DO 50 K = KDN + 1, 2, -1
                  J1 = J1 + KDN
                  J2 = J2 + KDN
*
                  IF( NR.GT.0 ) THEN
*
*                    generate plane rotations to annihilate nonzero
*                    elements which have been created outside the band
*
                     CALL ZLARGV( NR, AB( 1, J1-1 ), KD1*LDAB,
     $                            WORK( J1 ), KD1, D( J1 ), KD1 )
*
*                    apply rotations from the right
*
                     DO 10 L = 1, KD - 1
                        CALL ZLARTV( NR, AB( L+1, J1-1 ), KD1*LDAB,
     $                               AB( L, J1 ), KD1*LDAB, D( J1 ),
     $                               WORK( J1 ), KD1 )
   10                CONTINUE
                  END IF
*
                  IF( K.GT.2 ) THEN
                     IF( K.LE.N-I+1 ) THEN
*
*                       generate plane rotation to annihilate a(i,i+k-1)
*                       within the band
*
                        CALL ZLARTG( AB( KD-K+3, I+K-2 ),
     $                               AB( KD-K+2, I+K-1 ), D( I+K-1 ),
     $                               WORK( I+K-1 ), TEMP )
                        AB( KD-K+3, I+K-2 ) = TEMP
*
*                       apply rotation from the right
*
                        CALL ZROT( K-3, AB( KD-K+4, I+K-2 ), 1,
     $                             AB( KD-K+3, I+K-1 ), 1, D( I+K-1 ),
     $                             WORK( I+K-1 ) )
                     END IF
                     NR = NR + 1
                     J1 = J1 - KDN - 1
                  END IF
*
*                 apply plane rotations from both sides to diagonal
*                 blocks
*
                  IF( NR.GT.0 )
     $               CALL ZLAR2V( NR, AB( KD1, J1-1 ), AB( KD1, J1 ),
     $                            AB( KD, J1 ), KD1*LDAB, D( J1 ),
     $                            WORK( J1 ), KD1 )
*
*                 apply plane rotations from the left
*
                  CALL ZLACGV( NR, WORK( J1 ), KD1 )
                  DO 20 L = 1, KD - 1
                     IF( J2+L.GT.N ) THEN
                        NRT = NR - 1
                     ELSE
                        NRT = NR
                     END IF
                     IF( NRT.GT.0 )
     $                  CALL ZLARTV( NRT, AB( KD-L, J1+L ), KD1*LDAB,
     $                               AB( KD-L+1, J1+L ), KD1*LDAB,
     $                               D( J1 ), WORK( J1 ), KD1 )
   20             CONTINUE
*
                  IF( WANTQ ) THEN
*
*                    accumulate product of plane rotations in Q
*
                     DO 30 J = J1, J2, KD1
                        CALL ZROT( N, Q( 1, J-1 ), 1, Q( 1, J ), 1,
     $                             D( J ), DCONJG( WORK( J ) ) )
   30                CONTINUE
                  END IF
*
                  IF( J2+KDN.GT.N ) THEN
*
*                    adjust J2 to keep within the bounds of the matrix
*
                     NR = NR - 1
                     J2 = J2 - KDN - 1
                  END IF
*
                  DO 40 J = J1, J2, KD1
*
*                    create nonzero element a(j-1,j+kd) outside the band
*                    and store it in WORK
*
                     WORK( J+KD ) = WORK( J )*AB( 1, J+KD )
                     AB( 1, J+KD ) = D( J )*AB( 1, J+KD )
   40             CONTINUE
   50          CONTINUE
   60       CONTINUE
         END IF
*
         IF( KD.GT.0 ) THEN
*
*           make off-diagonal elements real and copy them to E
*
            DO 70 I = 1, N - 1
               T = AB( KD, I+1 )
               ABST = ABS( T )
               AB( KD, I+1 ) = ABST
               E( I ) = ABST
               IF( ABST.NE.ZERO ) THEN
                  T = T / ABST
               ELSE
                  T = ONE
               END IF
               IF( I.LT.N-1 )
     $            AB( KD, I+2 ) = AB( KD, I+2 )*T
               IF( WANTQ ) THEN
                  CALL ZSCAL( N, DCONJG( T ), Q( 1, I+1 ), 1 )
               END IF
   70       CONTINUE
         ELSE
*
*           set E to zero if original matrix was diagonal
*
            DO 80 I = 1, N - 1
               E( I ) = ZERO
   80       CONTINUE
         END IF
*
*        copy diagonal elements to D
*
         DO 90 I = 1, N
            D( I ) = AB( KD1, I )
   90    CONTINUE
*
      ELSE
*
         IF( KD.GT.1 ) THEN
*
*           Reduce to complex Hermitian tridiagonal form, working with
*           the lower triangle
*
            NR = 0
            J1 = KDN + 2
            J2 = 1
*
            AB( 1, 1 ) = DBLE( AB( 1, 1 ) )
            DO 150 I = 1, N - 2
*
*              Reduce i-th column of matrix to tridiagonal form
*
               DO 140 K = KDN + 1, 2, -1
                  J1 = J1 + KDN
                  J2 = J2 + KDN
*
                  IF( NR.GT.0 ) THEN
*
*                    generate plane rotations to annihilate nonzero
*                    elements which have been created outside the band
*
                     CALL ZLARGV( NR, AB( KD1, J1-KD1 ), KD1*LDAB,
     $                            WORK( J1 ), KD1, D( J1 ), KD1 )
*
*                    apply plane rotations from one side
*
                     DO 100 L = 1, KD - 1
                        CALL ZLARTV( NR, AB( KD1-L, J1-KD1+L ),
     $                               KD1*LDAB, AB( KD1-L+1, J1-KD1+L ),
     $                               KD1*LDAB, D( J1 ), WORK( J1 ),
     $                               KD1 )
  100                CONTINUE
                  END IF
*
                  IF( K.GT.2 ) THEN
                     IF( K.LE.N-I+1 ) THEN
*
*                       generate plane rotation to annihilate a(i+k-1,i)
*                       within the band
*
                        CALL ZLARTG( AB( K-1, I ), AB( K, I ),
     $                               D( I+K-1 ), WORK( I+K-1 ), TEMP )
                        AB( K-1, I ) = TEMP
*
*                       apply rotation from the left
*
                        CALL ZROT( K-3, AB( K-2, I+1 ), LDAB-1,
     $                             AB( K-1, I+1 ), LDAB-1, D( I+K-1 ),
     $                             WORK( I+K-1 ) )
                     END IF
                     NR = NR + 1
                     J1 = J1 - KDN - 1
                  END IF
*
*                 apply plane rotations from both sides to diagonal
*                 blocks
*
                  IF( NR.GT.0 )
     $               CALL ZLAR2V( NR, AB( 1, J1-1 ), AB( 1, J1 ),
     $                            AB( 2, J1-1 ), KD1*LDAB, D( J1 ),
     $                            WORK( J1 ), KD1 )
*
*                 apply plane rotations from the right
*
                  CALL ZLACGV( NR, WORK( J1 ), KD1 )
                  DO 110 L = 1, KD - 1
                     IF( J2+L.GT.N ) THEN
                        NRT = NR - 1
                     ELSE
                        NRT = NR
                     END IF
                     IF( NRT.GT.0 )
     $                  CALL ZLARTV( NRT, AB( L+2, J1-1 ), KD1*LDAB,
     $                               AB( L+1, J1 ), KD1*LDAB, D( J1 ),
     $                               WORK( J1 ), KD1 )
  110             CONTINUE
*
                  IF( WANTQ ) THEN
*
*                    accumulate product of plane rotations in Q
*
                     DO 120 J = J1, J2, KD1
                        CALL ZROT( N, Q( 1, J-1 ), 1, Q( 1, J ), 1,
     $                             D( J ), WORK( J ) )
  120                CONTINUE
                  END IF
*
                  IF( J2+KDN.GT.N ) THEN
*
*                    adjust J2 to keep within the bounds of the matrix
*
                     NR = NR - 1
                     J2 = J2 - KDN - 1
                  END IF
*
                  DO 130 J = J1, J2, KD1
*
*                    create nonzero element a(j+kd,j-1) outside the
*                    band and store it in WORK
*
                     WORK( J+KD ) = WORK( J )*AB( KD1, J )
                     AB( KD1, J ) = D( J )*AB( KD1, J )
  130             CONTINUE
  140          CONTINUE
  150       CONTINUE
         END IF
*
         IF( KD.GT.0 ) THEN
*
*           make off-diagonal elements real and copy them to E
*
            DO 160 I = 1, N - 1
               T = AB( 2, I )
               ABST = ABS( T )
               AB( 2, I ) = ABST
               E( I ) = ABST
               IF( ABST.NE.ZERO ) THEN
                  T = T / ABST
               ELSE
                  T = ONE
               END IF
               IF( I.LT.N-1 )
     $            AB( 2, I+1 ) = AB( 2, I+1 )*T
               IF( WANTQ ) THEN
                  CALL ZSCAL( N, T, Q( 1, I+1 ), 1 )
               END IF
  160       CONTINUE
         ELSE
*
*           set E to zero if original matrix was diagonal
*
            DO 170 I = 1, N - 1
               E( I ) = ZERO
  170       CONTINUE
         END IF
*
*        copy diagonal elements to D
*
         DO 180 I = 1, N
            D( I ) = AB( 1, I )
  180    CONTINUE
      END IF
*
      RETURN
*
*     End of ZHBTRD
*
      END
