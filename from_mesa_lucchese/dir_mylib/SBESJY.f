C---------------------------------------------------------------------
C      END OF SUBROUTINE RICBES
C---------------------------------------------------------------------
C---------------------------------------------------------------------
      SUBROUTINE SBESJY  (X,LMAX, J,Y,JP,YP, IFAIL )
C---------------------------------------------------------------------
C   REAL SPHERICAL BESSEL FUNCTIONS AND X DERIVATIVES
C      J , Y , J', Y'      FROM L=0 TO L=LMAX
C      FOR REAL X > SQRT(ACCUR) (E.G. 1D-7)    AND INTEGER LMAX
C
C  J (L)  =      J/L/(X) STORES   REGULAR SPHERICAL BESSEL FUNCTION:
C  JP(L)  = D/DX J/L/(X)      J(0) =  SIN(X)/X
C  Y (L)  =      Y/L/(X) STORES IRREGULAR SPHERICAL BESSEL FUNCTION:
C  YP(L)  = D/DX Y/L/(X)      Y(0) = -COS(X)/X
C
C    IFAIL = -1 FOR ARGUMENTS OUT OF RANGE
C          =  0 FOR ALL RESULTS SATISFACTORY
C
C   USING LENTZ-THOMPSON EVALUATION OF CONTINUED FRACTION CF1,
C   AND TRIGONOMETRIC FORMS FOR L = 0 SOLUTIONS.
C   LMAX IS LARGEST L NEEDED AND MUST BE <= MAXL, THE ARRAY INDEX.
C   MAXL CAN BE DELETED AND ALL THE ARRAYS DIMENSIONED (0:*)
C   SMALL IS MACHINE DEPENDENT, ABOUT SQRT(MINIMUM REAL NUMBER),
C         SO 1D-150 FOR DOUBLE PRECISION ON VAX, PCS ETC.
C   PRECISION:  RESULTS TO WITHIN 2-3 DECIMALS OF "MACHINE ACCURACY"
C   IN OSCILLATING REGION X .GE.  [ SQRT{LMAX*(LMAX+1)} ]
C   I.E. THE TARGET ACCURACY ACCUR SHOULD BE 100 * ACC8 WHERE ACC8
C   IS THE SMALLEST NUMBER WITH 1+ACC8.NE.1 FOR OUR WORKING PRECISION
C   THIS RANGES BETWEEN 4E-15 AND 2D-17 ON CRAY, VAX, SUN, PC FORTRANS
C   SO CHOOSE A SENSIBLE  ACCUR = 1.0D-14
C   IF X IS SMALLER THAN [ ] ABOVE, THE ACCURACY BECOMES STEADLY WORSE:
C   THE VARIABLE ERR IN COMMON /STEED/ HAS AN ESTIMATE OF ACCURACY.
C
C   NOTE: FOR X=1 AND L=100  J = 7.4 E-190     Y = -6.7+E186
C---------------------------------------------------------------------
C   AUTHOR :   A.R.BARNETT      MANCHESTER    12 MARCH 1990.
C                               AUCKLAND      12 MARCH 1991.
C---------------------------------------------------------------------
      INTEGER     LIMIT,      MAXL,      LMAX, IFAIL, NFP, L
      PARAMETER ( LIMIT = 20000, MAXL = 1001 )
      DOUBLE PRECISION  J(0:MAXL), Y(0:MAXL), JP(0:MAXL), YP(0:MAXL)
      DOUBLE PRECISION  ZERO,ONE,TWO,THREE,SMALL, ACCUR, TK,SL, ERR
      DOUBLE PRECISION  X,XINV, CF1,DCF1, DEN, C,D, OMEGA, TWOXI
      PARAMETER ( ZERO  = 0.0D0  , ONE   = 1.0D0 , TWO = 2.0D0 )
      PARAMETER ( SMALL = 1.D-150, THREE = 3.0D0 )
      COMMON /STEED/    ERR,NFP    
C----
      ACCUR = 1.D-14      
      IFAIL = -1      
      IF (X .LT. SQRT(ACCUR) ) GOTO 50
C---- TAKES CARE OF NEGATIVE X ... USE REFLECTION FORMULA
C---- BEGIN CALCULATION OF CF1 UNLESS LMAX = 0, WHEN SOLUTIONS BELOW
      XINV = ONE / X
      IF (LMAX .GT. 0) THEN
       TWOXI = XINV + XINV
       SL = REAL(LMAX)* XINV 
       TK = TWO * SL + XINV * THREE
       CF1 = SL 
       DEN = ONE 
       IF ( ABS(CF1) .LT. SMALL ) CF1 = SMALL
C---- INVERSE RATIO OF A CONVERGENTS
       C = CF1 
C---- DIRECT  RATIO OF B CONVERGENTS
       D = ZERO 
       DO 10 L = 1,LIMIT
        C = TK - ONE / C
        D = TK - D
        IF ( ABS(C) .LT. SMALL ) C = SMALL
        IF ( ABS(D) .LT. SMALL ) D = SMALL
        D = ONE / D
        DCF1= D * C
        CF1 = CF1 * DCF1
        IF ( D .LT. ZERO ) DEN = - DEN
        IF ( ABS(DCF1 - ONE) .LE. ACCUR ) GOTO 20
        TK = TK + TWOXI
10     CONTINUE
C---- ERROR EXIT, NO CONVERGENCE
       GOTO 50
20     NFP = L 
C---- ERROR ESTIMATE
       ERR = ACCUR*SQRT(DBLE(NFP))
       J (LMAX) = DEN 
       JP(LMAX) = CF1 * DEN
C---- DOWNWARD RECURSION TO L=0  AS SPHERICAL BESSEL FUNCTIONS
       DO 30 L = LMAX , 1, -1
        J (L-1) = (SL + XINV) * J(L) + JP(L)
        SL = SL - XINV
        JP(L-1) = SL * J(L-1) - J(L)
30     CONTINUE
       DEN = J(0)
C---- END LOOP FOR LMAX > 0
       ENDIF
C---- CALCULATE THE L=0 SPHERICAL BESSEL FUNCTIONS DIRECTLY
       J (0) = XINV * SIN(X)
       Y (0) = -XINV * COS(X)
       JP(0) = -Y(0) - XINV * J(0)
       YP(0) = J(0) - XINV * Y(0)
       IF (LMAX .GT. 0) THEN
        OMEGA = J(0) / DEN
        SL = ZERO
        DO 40 L = 1 , LMAX
         J (L) = OMEGA * J (L)
         JP(L) = OMEGA * JP(L)
         Y (L) = SL * Y(L-1) - YP(L-1)
         SL = SL + XINV
         YP(L) = Y(L-1) - (SL + XINV) * Y(L)
40      CONTINUE
       ENDIF
C---- CALCULATIONS SUCCESSFUL
       IFAIL = 0 
       RETURN
C---------------------------------------------------------------------
C---- ERROR TRAPS
C---------------------------------------------------------------------
50     IF (X .LT. ZERO) THEN
        WRITE(6,1000) X
       ELSE
     -  IF (X .EQ. ZERO) THEN
         IFAIL = 0
         J(0) = ONE
         DO 60 L = 1, LMAX
          J(L) = ZERO 
60       CONTINUE
         ELSE 
         WRITE(6,1001) X
       ENDIF
1000  FORMAT(' X NEGATIVE !',1P,E15.5,' ... USE REFLECTION FORMULA'/)          
1001  FORMAT(' WITH X = ',1P,E15.5,'    TRY SMALL-X SOLUTIONS',/,              
     X  '    J/L/(X)  ->   X**L / (2L+1)!!      AND',/,                         
     X  '    Y/L/(X)  ->  -(2L-1)!! / X**(L+1)'/)                               
      RETURN
      END
