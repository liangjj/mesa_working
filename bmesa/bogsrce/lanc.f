C*****************************************************************
C*****************************************************************
C           AN ITERATIVE BLOCK LANCZOS METHOD FOR THE
C              SOLUTION OF LARGE SPARSE SYMMETRIC
C                       EIGENPROBLEMS
C                            by
C                   Richard Ray Underwood
C
C
C     if further details concerning the content of the code are
C     required, you are advised to refer to R. R. Underwood's
C     PhD thesis - "An Iterative Block Lanczos Method for the
C     Solution of Large Sparse Symmetric Eigenproblems", 1975.
C*****************************************************************
C*****************************************************************
      SUBROUTINE MINVAL(SIZE,IA,JA,SA,N,Q,PINIT,R,MMAX,EPS,M,
     1D,X,IECODE)
C=================================================================
C     this subroutine is the main subroutine implementing the
C     iterative block lanczos method for computing the
C     eigenvalues and eigenvectors of symmetric matrices.
C=================================================================
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER SIZE,N,Q,PINIT,R,MMAX
      DOUBLE PRECISION EPS
      DIMENSION D(Q),X(N,Q)
      INTEGER IA(SIZE),JA(SIZE) 
      DOUBLE PRECISION SA(SIZE)
      INTEGER IECODE
      DIMENSION E(25),C(25,25)
      parameter (nmax = 800000)
      DIMENSION U(nmax),V(nmax)
      INTEGER P,S,PS
C
C------
C     description of parameters:
C     N:     integer variable. The order of the symmetric
C            matrix A whose eigenvalues and eigenvectors are
C            being computed. The value of N should be less than
C            or equal to nmax and greater than or equal to  2.
C
C     Q:     integer variable. The number of vectors of length
C            N contained in the array X. The value of Q should
C            be less than or equal to 25, at least one greater
C            than the value of R and less than or equal to N.
C
C     PINIT: integer variable. The initial block size to be used
C            in the block lanczos method. If PINIT is negative,
C            then -PINIT is used for the block size and columns
C            M+L,...,M+(-PINIT) of the array X are assumed to be
C            initialized to the initial matrix used to start the
C            block lanczos method. If the subroutine terminates
C            with a value of M less than R, then PINIT is assigned
C            a value -P where P is the final block size chosen.
C            In this circumstance, columns M+1,...M+P will contain
C            the most recent set of eigenvector approximations
C            which can be used to restart the subroutine if desired.
C
C     R:     integer variable. The number of eigenvalues and
C            eigenvectors being computed. That is, MINVAL attempts
C            to compute accurate approximations to the R least
C            eigenvalues and eigenvectors of the matrix A. The value
C            of R should be greater than zero and less than Q.
C
C     MMAX:  integer variable. The maximum number of matrix-vector
C            products A*X (where X is a vector) that are allowed
C            during one call of this subroutine to complete its task
C            of computing R eigenvalues and eigenvectors. Unless the
C            problem indicates otherwise, MMAX should be given a very
C            large value.
C
C     EPS:   DOUBLE PRECISION variable. Initially, EPS should contain a value
C            indicating the relative precision to which MINVAL will
C            attempt to compute the eigenvalues and eigenvectors of A.
C            For eigenvalues less in modulus than 1, EPS will be an
C            absolute tolerance. Because of the way this method works,
C            it may happen that the later eigenvalues cannot be 
C            computed to the same relative precision as those less in
C            value.
C
C     OP:    subroutine name. The actual argument corresponding to OP
C            should be the name of a subroutine used to define the
C            matrix A. This subroutine should have three arguments
C            N, U, and V, say, where N is an integer variable giving
C            the order of A, and U and V are two one-dimensional
C            arrays of length N. If W denotes the vector of order N
C            stored in U, then the statement
C                    CALL OP(N,U,V)
C            should result in the vector A*W being computed and stored
C            in V. The contents of U can be modified by this call.
C
C     M:     integer variable. M gives the number of eigenvalues and
C            eigenvectors already computed. Thus, initially, M should
C            be zero. If M is greater than zero, then columns one
C            through M of the array X are assumed to contain the 
C            computed approximations to the M least eigenvalues and
C            eigenvectors of A. On exit, M contains a value equal to
C            the total number of eigenvalues and eigenvectors
C            computed including any already computed when MINVAL was
C            entered. Thus, on exit, the first M elements of D and the
C            first M columns of X will contain approximations to the
C            M least eigenvalues of A.
C
C     D:     DOUBLE PRECISION array. D contains the computed eigenvalues. D
C            should be a one-dimensional array with at least Q
C            elements.
C
C     X:     DOUBLE PRECISION array. X contains the computed eigenvectors. X
C            should be an array containing at least N*Q elements. X
C            is used not only to store the eigenvectors computed by
C            MINVAL, but also as working storage for the block lanczos
C            method. On exit, the first N*M elements of X contain the
C            eigenvector approximations - the first vector in the 
C            first N elements, the second in the second N elements,
C            etc...
C
C
C     IECODE:integer variable. The value of IECODE indicates whether
C            MINVAL terminated successfully, and if not, the reason
C            why.
C
C               IECODE=0 : successful termination.
C               IECODE=1 : the value of N is less than 2.
C               IECODE=2 : the value of N exceeds nmax.
C               IECODE=3 : the value of R is less than 1.
C               IECODE=4 : the value of Q is less than or equal to R.
C               IECODE=5 : the value of Q is greater than 25.
C               IECODE=6 : the value of Q exceeds N.
C               IECODE=7 : the value of MMAX was exceeded before R
C                          eigenvalues and eigenvectors were
C                          computed.
C
C      Note that the subroutine has been designed to allow initial
C      approximations to the eigenvectors corresponding to the least
C      eigenvalues to be utilised if they are known (by storing them
C      in X and assigning PINIT minus the value of their number).
C      Furthermore, it has also been designed to allow restarting if
C      it stops with IECODE=7. Thus, the user of this program can
C      restart it after examining any partial results without loss of
C      previous work.
C------
C
C------
C     check that the initial values of the subroutine parameters
C     are in range.
C------ 
      IF (N.LT.2) GO TO 901
      IF (N.GT.nmax) GO TO 902
      IF (R.LT.1) GO TO 903
      IF (Q.LE.R) GO TO 904
      IF (Q.GT.25) GO TO 905
      IF (Q.GT.N) GO TO 906

C
C------
C     choose initial values for the block size P, the number of
C     steps that the block lanczos method is carried out, and
C     choose an initial N-by-P orthonormal matrix X1 used to start
C     the block lanczos method.
C------      
      P = PINIT
      IF (P.LT.0) P = -P
      S = (Q-M)/P
      IF (S.GT.2) GO TO 100
      S = 2
      P = Q/2
  100 IF (PINIT.LT.0) GO TO 150
      DO 120 K = 1,P
             CALL RANDOM(N,Q,K,X)
  120 CONTINUE
  150 IF (M.GT.0) GO TO 200
      CALL ORTHG(N,Q,M,P,C,X)
C
C------
C     rotate the initial N-by-P matrix X1 so that
C        X1'*A*X1=diag(D1,D2,...,DP)
C     where DI is stored in D(I), I=1,...,P.
C------ 
      CALL SECTN(SIZE,IA,JA,SA,N,Q,M,P,X,C,D,U,V)
      ERRC = 0.D0
  200 ITER = 0
      IMM=0
C
C------
C     the main body of the subroutine starts here. IMM
C     counts the number of matrix-vector products computed,
C     which is the number of times the subroutine named by
C     OP is called. ERRC measures the accumulated error in
C     the eigenvalues and eigenvectors.
C------
C
  300 IF (M.GE.R) GO TO 900
      IF (IMM.GT.MMAX) GO TO 907
      ITER = ITER + 1
      PS = P*S
C
C------
C     BKLANC carries out the block lanczos method and stores
C     the resulting block tridiagonal matrix MS in C and the
C     N-by-PS orthonormal matrix XS in X. The initial N-by-P
C     orthonormal matrix is assumed to be stored in columns
C     M+1 through M+PS of X. The residuals for these vectors
C     and the eigenvalue approximations in D are computed and
C     stored in E.
C------
      CALL BKLANC(SIZE,IA,JA,SA,N,Q,M,P,S,D,C,X,E,U,V)
C
C------
C     EIGEN solves the eigenproblem for MS, storing the eigenvalues
C     in elements M+1 through M+PS of D and the eigenvectors in the
C     first P*S rows and columns of C (overwriting MS, possibly).
C------
      CALL EIGEN(Q,M,P,PS,C,D)
C
C------
C     CNVTST determines how many of the eigenvalues and eigenvectors
C     have converged using the error estimates stored in E. The number
C     that have converged is stored in NCONV. If NCONV=0, then none
C     have converged.
C------
      CALL CNVTST(N,Q,M,P,ERRC,EPS,D,E,NCONV)
C     write(6,*) nconv
C
C------
C     PCH chooses new values for P and S, the block size and the
C     number of steps for the block lanczos subprogram, respectively.
C------
      CALL PCH(N,Q,M,R,NCONV,P,S)
C
C------
C     ROTATE computes the eigenvectors of the restricted matrix
C     using XS stored in X and the eigenvectors of MS stored in C.
C     These vectors serve both as eigenvector approximations and
C     to form the matrix used to start the block lanczos method in
C     the next iteration.
C------
      CALL ROTATE(N,Q,M,PS,NCONV+P,C,X)
C
      M = M + NCONV
      IMM=IMM+P*S
C      PRINT 1001,ITER,IMM,P,PS
C 1001 FORMAT(' =>ITER,IMM,P,PS =',4I5)
      GO TO 300
C
C------
C     this is the end of the main body of the subroutine. Now set
C     the value of IECODE and EXIT.
C------
  900 IECODE = 0
      RETURN
  901 IECODE = 1
      RETURN
  902 IECODE = 2
      RETURN
  903 IECODE = 3
      RETURN
  904 IECODE = 4
      RETURN
  905 IECODE = 5
      RETURN
  906 IECODE = 6
      RETURN
  907 IECODE = 7
      PINIT = -P
C
      RETURN
      END
C
C
C
C
      SUBROUTINE BKLANC(SIZE,IA,JA,SA,N,Q,M,P,S,D,C,X,E,U,V)
C====================================================================
C     this subroutine implements the block lanczos method
C     with reorthogonalization. BKLANC computes a block
C     tridiagonal matrix MS which it stores in rows and
C     columns M+1 through M+P*S of the array C, and an
C     orthonormal matrix XS which it stores in columns M+1
C     through M+P*S of the N-by-Q array X. MS is a symmetric
C     matrix with P-by-P symmetric matrices M(1),...,M(S) on
C     its diagonal and P-by-P upper triangular matrices
C     R(2),...,R(S) along its lower diagonal. Since MS is 
C     symmetric and banded, only its lower triangle (P+1
C     diagonals) is stored in C. XS is composed of S N-by-P
C     orthonormal matrices X(1),...,X(S) where X(1) is given
C     and should be stored in columns M+1 through M+P of X.
C     Furthermore, X(1) is assumed to satisfy X(1)*A*X(1) =
C     diag(D(M+1),D(M+2),...,D(M+P)), and if M>0, then X(1) is
C     assumed to be orthogonal to the vectors stored in columns
C     1 through M of X. OP is the name of an external subroutine
C     used to define the matrix A. During the first step, the
C     subroutine ERR is called and the quantities EJ are computed
C     where EJ=||A*X1J-D(M+J)*X1J||, X1J is the J-th column of x(1),
C     and ||*|| denotes the Euclidean norm. EJ is stored in E(M+J),
C     J=1,2,...,P. U and V are auxilliary vectors used by OP.
C====================================================================
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER SIZE,N,Q,M,P,S
      DIMENSION D(Q),C(Q,Q),X(N,Q)
      INTEGER IA(SIZE),JA(SIZE)
      DOUBLE PRECISION SA(SIZE)
      DIMENSION E(Q),U(N),V(N)
C
      MP1 = M+1
      MPPS = M+P*S
      DO 90 L = 1,S
            LL = M+(L-1)*P+1
            LU = M + L*P
            DO 70 K = LL,LU
                  DO 10 I = 1,N
   10                   U(I) = X(I,K)
                  CALL OP(SIZE,IA,JA,SA,N,U,V)
                  IF (L.GT.1) GO TO 19
                  DO 12 I=K,LU
   12                   C(I,K) = 0.D0
                  C(K,K) = D(K)
                  DO 14 I = 1,N
   14                   V(I) = V(I) - D(K)*X(I,K)
                  GO TO 61
   19             DO 30 I = K,LU
                        T = 0
                        DO 20 J = 1,N
   20                         T = T+V(J)*X(J,I)
   30             C(I,K) = T
                  IT = K-P
                  DO 60 I = 1,N
                        T = 0
                        DO 40 J = IT,K
   40                         T = T + X(I,J)*C(K,J)
                  IF (K.EQ.LU) GO TO 60
                  KP1 = K+1
                  DO 50 J = KP1,LU
   50             T = T+X(I,J)*C(J,K)
   60             V(I) = V(I)-T
   61             IF (L.EQ.S) GO TO 70
                  DO 63 I = 1,N
   63                   X(I,K+P) = V(I)
   70    CONTINUE
         IF (L.EQ.1) CALL ERR(N,Q,M,P,X,E)
         IF (L.EQ.S) GO TO 90
         CALL ORTHG(N,Q,LU,P,C,X)
         IL = LU+1
         IT = LU
         DO 80 J = 1,P
               IT = IT+1
               DO 80 I = IL,IT
   80                C(I,IT-P) = C(I,IT)
   90 CONTINUE
C
      RETURN
      END
C
C
C
C
      SUBROUTINE PCH(N,Q,M,R,NCONV,P,S)
C====================================================================
C     based on the values of N, Q, M, R and NCONV, PCH chooses new
C     values for P and S, the block size and number of steps for the
C     block lanczos method. The strategy used here is to choose P to
C     be the smaller of the two following values: 
C            1) the previous block size
C     and,   2) the number of values left to be computed. S is chosen
C     as large as possible subject to the constraints imposed by the
C     limits of storage. In any event, S is greater than or equal to
C     2. N is the order of the problem and Q is the number of vectors
C     available for storing eigenvectors and applying the block
C     lanczos method. M is the number of eigenvalues and eigenvectors
C     that have already been computed and R is the required number.
C     Finally, NCONV is the number of eigenvalues and eigenvectors
C     that have converged in the current iteration.
C====================================================================
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER N,Q,M,R,NCONV,P,S
      INTEGER PT,ST
C
      MT = M+NCONV
      PT = R-MT
      IF (PT.GT.P) PT = P
      IF (PT.GT.0) GO TO 101
      P = 0
      RETURN
  101 CONTINUE
      ST = (Q-MT)/PT
      IF (ST.GT.2) GO TO 110
      ST = 2
      PT = (Q-MT)/2
  110 P = PT
      S = ST
C
      RETURN
      END
C
C
C
C
      SUBROUTINE ERR(N,Q,M,P,X,E)
C================================================
C     ERR computes the Euclidean lengths of the 
C     vectors stored in the columns M+P+1 through
C     M+P+P of the N-by-Q array X and stores them
C     in elements M+1 through M+P of E.
C================================================
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER N,Q,M,P
      DIMENSION X(N,Q),E(Q)
C
      MP1 = M+P+1
      MPP = M+P+P
      DO 200 K = MP1,MPP
             T = 0.D0
             DO 100 I = 1,N
  100               T = T+X(I,K)**2
  200        E(K-P) = DSQRT(T)
C
      RETURN
      END
C
C
C
C
      SUBROUTINE CNVTST(N,Q,M,P,ERRC,EPS,D,E,NCONV)
C===========================================================
C     CNVTST determines which of the P eigenvalues stored
C     in elements M+1 through M+P of D have converged. ERRC
C     is a measure of the accumulated error in the M
C     previously computed eigenvalues and eigenvectors. ERRC
C     is updated if more approximations have converged. The 
C     norms of the residual vectors are stored in elements
C     M+1 through M+P of E. EPS is the precision to which we
C     are computing the approximations. Finally, NCONV is the
C     number that have converged. If NCONV=0, then none have
C     converged.
C============================================================
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER Q,M,P
      DOUBLE PRECISION EPS
      DIMENSION D(Q),E(Q)
      DATA CHEPS/2.22D-16/
C
      K = 0
      DO 100 I = 1,P
         T = DABS(D(M+I))
         IF (T.LT.1.D0) T = 1.D0
      IF (E(M+I).GT.T*(EPS+10.D0*N*CHEPS)+ERRC) GO TO 200
  100 K = I
  200 NCONV = K
      IF (K.EQ.0) RETURN
      T = 0.D0
      DO 300 I = 1,K
             T = T+E(M+I)**2
  300 CONTINUE
      ERRC = DSQRT(ERRC**2+T)
C
      RETURN
      END
C
C
C
C
      SUBROUTINE EIGEN(Q,M,P,PS,C,D)
C=======================================================
C     EIGEN solves the eigenproblem for the symmetric
C     matrix MS of order PS stored in rows and columns
C     M+1 through M+PS of C. The eigenvalues of MS are
C     stored in elements M+1 through M+PS of D and the
C     eigenvactors are stored in rows and columns 1 
C     through PS of C possibly overwriting MS. EIGEN
C     simply re-stores MS in a manner acceptable to the
C     subroutines TRED2 and TQL2. These two routines are
C     available through eispack.
C=======================================================
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER M,P,Q,PS,PP1
      DIMENSION C(Q,Q),D(Q)
      DIMENSION DD(25),V(25)
C
      DO 150 I=1,PS
             LIM = I-P
             IF (I.LE.P) LIM = 1
             IF (LIM.LE.1) GO TO 130
             LM1 = LIM-1
         DO 120 J = 1,LM1
  120           C(I,J) = 0.D0
  130    DO 140 J = LIM,I
  140           C(I,J) = C(I+M,J+M)
  150 CONTINUE
C
C------
C     call eispack routines here
C------
      CALL TRED2(Q,PS,C,DD,V,C)
      CALL TQL2(Q,PS,DD,V,C,IERR)
C
C      PRINT 601,PS,(DD(I),I=1,PS)
C  601 FORMAT(' => ORDER =',I4,/,('    EIGENVALUES =',10d10.4))
      PP1 = P+1
ccc      DO 155 J = 1,PP1
ccc  155 PRINT 602,J,(C(I,J),I=1,PS)
  602 FORMAT(' => J =',I4,/,('    EIGENVECTORS =',10d10.4))
      DO 160 I = 1,PS
  160        D(M+I) = DD(I)
C
      RETURN
      END
C
C
C
C
      SUBROUTINE SECTN(SIZE,IA,JA,SA,N,Q,M,P,X,C,D,U,V)
C==============================================================
C     SECTN transforms the N-by-P orthonormal matrix X1,
C     say, stored in columns M+1 through M+P of the N-by-Q
C     array X so that X1'*A*X1 = diag(D1,D2,...,DP), where
C     ' denotes transpose and A is a symmetric matrix of 
C     order N defined by the subroutine OP. The values D1,...
C     ,DP are stored in elements M+1 through M+P of D. SECTN
C     forms the matrix X1'*A*X1 = CP, storing CP in the array
C     C. The values D1,D2,...,DP and the eigenvectors QP of CP
C     are computed by EIGEN and stored in D and C respectively.
C     ROTATE then carries out the transformation X1<=X1*QP.
C==============================================================
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER SIZE,N,Q,M,P
      DIMENSION X(N,Q),C(Q,Q),D(Q),U(N),V(N)
      INTEGER IA(SIZE),JA(SIZE)
      DOUBLE PRECISION SA(SIZE)
C
      ICOL1 = M
      DO 300 J = 1,P
             ICOL1 = ICOL1 + 1
         DO 100 I = 1,N
  100           U(I) = X(I,ICOL1)
                CALL OP(SIZE,IA,JA,SA,N,U,V)
                ICOL2 = M
             DO 300 I = 1,J
                    ICOL2 = ICOL2+1
                    T = 0.D0
                DO 200 K = 1,N
  200                  T = T + V(K)*X(K,ICOL2)
  300 C(ICOL1,ICOL2) = T
      CALL EIGEN(Q,M,P,P,C,D)
      CALL ROTATE(N,Q,M,P,P,C,X)
C
      RETURN
      END
C
C
C
C
      SUBROUTINE ROTATE(N,Q,M,PS,L,C,X)
C======================================================
C     ROTATE computes the first L columns of the matrix
C     XS*QS where XS is an N-by-PS orthonormal matrix 
C     stored in columns M+1 through M+PS of the N-by-Q
C     array X and QS is a PS-by-PS orthonormal matrix
C     stored in rows and columns 1 through Ps of the
C     array C. The result is stored in columns M+1
C     through M+L of X overwriting part of XS.
C======================================================
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER N,Q,M,PS,L
      DIMENSION C(Q,Q),X(N,Q)
      DIMENSION V(25)
C
      DO 300 I = 1,N
         DO 200 K = 1,L
                T = 0.D0
            DO 100 J = 1,PS
  100              T = T+X(I,M+J)*C(J,K)
  200       V(K) = T
         DO 300 K = 1,L
  300 X(I,M+K) = V(K)
C
      RETURN
      END
C
C
C
C
      SUBROUTINE ORTHG(N,Q,F,P,B,X)
C=========================================================
C     ORTHG reorthogonalizes the N-by-P matrix Z stored
C     in columns F+1 through F+P of the N-by-Q array X
C     with respect to the vectors stored in the first F
C     columns of X and then decomposes the resulting
C     matrix into the product of an N-by-P orthonormal
C     matrix XORTH, say, and a P-by-P upper triangular
C     matrix R. XORTH is stored over Z and the upper
C     triangle of R is stored in rows and columns F+1
C     through F+P of the Q-by-Q array B. A stable variant
C     of the Gram-Schmidt orthogonalization method is
C     utilised.
C=========================================================
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER N,Q,F,P
      DIMENSION B(Q,Q),X(N,Q)
      INTEGER FP1,FPP
      LOGICAL*1 ORIG
C
      IF (P.EQ.0) RETURN
      FP1 = F+1
      FPP = F+P
      DO 50 K = FP1,FPP
            ORIG = .TRUE.
            KM1 = K-1
   10       T = 0.D0
            IF (KM1.LT.1) GO TO 25
         DO 20 I = 1,KM1
               S = 0.D0
            DO 15 J = 1,N
   15          S = S+X(J,I)*X(J,K)
            IF (ORIG.AND.I.GT.F) B(I,K) = S
            T = T + S*S
         DO 20 J = 1,N
   20          X(J,K) = X(J,K) - S*X(J,I)
   25    S = 0.D0
         DO 30 J = 1,N
   30          S = S + X(J,K)*X(J,K)
         T = T+S
         IF (S.GT.T/100.D0) GO TO 40
         ORIG = .FALSE.
         GO TO 10
   40    S = DSQRT(S)
         B(K,K) = S
         IF (S.NE.0) S = 1.D0/S
         DO 50 J = 1,N
   50 X(J,K) = S*X(J,K)
C
      RETURN
      END
C
C
C
C
      SUBROUTINE RANDOM(N,Q,L,X)
C======================================================
C     RANDOM computes and stores a sequence of N
C     pseudo-random numbers in the L-th column of the
C     N-by-Q array X. RANDOM generates two sequences of
C     pseudo-random numbers, filling an array with one
C     sequence and using the second to access the array 
C     in a random fashion.
C======================================================
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER N,Q,L,FT,X1,F1,F2
      INTEGER A,C,X0
      DIMENSION X(N,Q)
      DIMENSION T(100)
      DATA F1/71416/,F2/27183/
      DATA A/6821/,C/5327/,X0/5328/
C
      DO 100 I = 1,100
            X1 = A *X0+C
            IF (X1.GE.10000) X1 = X1 - 10000
            T(I) = X1/9999.D0 - 0.5D0
  100 X0 = X1
      DO 200 I = 1,N
             FT = F1+F2
             IF (FT.GE.1000000) FT = FT-1000000
             F1 = F2
             F2 = FT
             K = FT/1.D6*100 + 1
             X(I,L) = T(K)
             X1 = A*X0 + C
             IF (X1.GE.10000) X1 = X1 - 10000
             T(K) = X1/9999D0 - 0.5D0
  200    X0 = X1
C
      RETURN
      END
C
C
C
C
      SUBROUTINE OP(SIZE,IA,JA,SA,N,U,V)
C======================================================
C     user dependent subroutine. It is defined to be
C     the actual argument corresponding to OP and is
C     used to define the matrix A. This subroutine
C     should have three arguments N, U and V, say,
C     where N is an integer variable giving the order
C     of A, and U and V are two one-dimensional arrays
C     of length N.
C======================================================
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER SIZE,N
      DIMENSION U(N),V(N)
      INTEGER IA(SIZE),JA(SIZE)
      DOUBLE PRECISION SA(SIZE)

C  
C     The non-zero elements of A are: 
C            row of A is in IA[i] 
C            column of A is in JA[i]
C
C     The actual matrix element is found in SA[i].
C
C     Note that the FORTRAN-C interface is inefficient at passing
C     two-dimenionsal arrays; therefore, all information is found in
C     three one-dimensional arrays.
C
      DO 100 I=1,N
100	 V(I)=0.0
      DO 90 I=1,SIZE
         V(IA(I)) = V(IA(I)) + SA(I)*U(JA(I))
90       IF (IA(I).NE.JA(I)) V(JA(I)) = V(JA(I)) + SA(I)*U(IA(I))
C
      RETURN
      END
      double precision function epslon (x)
      double precision x
c
c     estimate unit roundoff in quantities of size x.
c
      double precision a,b,c,eps
c
c     this program should function properly on all systems
c     satisfying the following two assumptions,
c        1.  the base used in representing floating point
c            numbers is not a power of three.
c        2.  the quantity  a  in statement 10 is represented to
c            the accuracy used in floating point variables
c            that are stored in memory.
c     the statement number 10 and the go to 10 are intended to
c     force optimizing compilers to generate code satisfying
c     assumption 2.
c     under these assumptions, it should be true that,
c            a  is not exactly equal to four-thirds,
c            b  has a zero for its last bit or digit,
c            c  is not exactly equal to one,
c            eps  measures the separation of 1.0 from
c                 the next larger floating point number.
c     the developers of eispack would appreciate being informed
c     about any systems where these assumptions do not hold.
c
c     this version dated 4/6/83.
c
      a = 4.0d0/3.0d0
   10 b = a - 1.0d0
      c = b + b + b
      eps = dabs(c-1.0d0)
      if (eps .eq. 0.0d0) go to 10
      epslon = eps*dabs(x)
      return
      end

