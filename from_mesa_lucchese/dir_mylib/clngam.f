*DECK CLNGAM
      COMPLEX FUNCTION CLNGAM (ZIN)
C***BEGIN PROLOGUE  CLNGAM
C***PURPOSE  Compute the logarithm of the absolute value of the Gamma
C            function.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C7A
C***TYPE      COMPLEX (ALNGAM-S, DLNGAM-D, CLNGAM-C)
C***KEYWORDS  ABSOLUTE VALUE, COMPLETE GAMMA FUNCTION, FNLIB, LOGARITHM,
C             SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C CLNGAM computes the natural log of the complex valued gamma function
C at ZIN, where ZIN is a complex number.  This is a preliminary version,
C which is not accurate.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  C9LGMC, CARG, CLNREL, R1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   780401  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  CLNGAM
      COMPLEX ZIN, Z, CORR, CLNREL, C9LGMC
      LOGICAL FIRST
      external carg
      SAVE PI, SQ2PIL, BOUND, DXREL, FIRST
      DATA PI / 3.1415926535 8979324E0 /
      DATA SQ2PIL / 0.9189385332 0467274E0 /
      DATA FIRST /.TRUE./

C***FIRST EXECUTABLE STATEMENT  CLNGAM
      IF (FIRST) THEN
         N = -0.30*LOG(R1MACH(3))
C BOUND = N*(0.1*EPS)**(-1/(2*N-1))/(PI*EXP(1))
         BOUND = 0.1171*N*(0.1*R1MACH(3))**(-1./(2*N-1))
         DXREL = SQRT (R1MACH(4))
      ENDIF
      FIRST = .FALSE.
C
      Z = ZIN
      X = REAL(ZIN)
      Y = AIMAG(ZIN)
C
      CORR = (0.0, 0.0)
      CABSZ = ABS(Z)
      IF (X.GE.0.0 .AND. CABSZ.GT.BOUND) GO TO 50
      IF (X.LT.0.0 .AND. ABS(Y).GT.BOUND) GO TO 50
C
      IF (CABSZ.LT.BOUND) GO TO 20
C
C USE THE REFLECTION FORMULA FOR REAL(Z) NEGATIVE, ABS(Z) LARGE, AND
C ABS(AIMAG(Y)) SMALL.
C
      IF (Y.GT.0.0) Z = CONJG (Z)
      CORR = EXP (-CMPLX(0.0,2.0*PI)*Z)
      IF (REAL(CORR) .EQ. 1.0 .AND. AIMAG(CORR) .EQ. 0.0) CALL XERMSG
     +   ('SLATEC', 'CLNGAM', 'Z IS A NEGATIVE INTEGER', 3, 2)
C
      CLNGAM = SQ2PIL + 1.0 - CMPLX(0.0,PI)*(Z-0.5) - CLNREL(-CORR)
     1  + (Z-0.5)*LOG(1.0-Z) - Z - C9LGMC(1.0-Z)
      IF (Y.GT.0.0) CLNGAM = CONJG (CLNGAM)
      RETURN
C
C USE THE RECURSION RELATION FOR ABS(Z) SMALL.
C
 20   IF (X.GE.(-0.5) .OR. ABS(Y).GT.DXREL) GO TO 30
      IF (ABS((Z-AINT(X-0.5))/X) .LT. DXREL) CALL XERMSG ('SLATEC',
     +   'CLNGAM',
     +   'ANSWER LT HALF PRECISION BECAUSE Z TOO NEAR NEGATIVE INTEGER',
     +   1, 1)
C
 30   N = SQRT (BOUND**2 - Y**2) - X + 1.0
      ARGSUM = 0.0
      CORR = (1.0, 0.0)
      DO 40 I=1,N
        ARGSUM = ARGSUM + CARG(Z)
        CORR = Z*CORR
        Z = 1.0 + Z
 40   CONTINUE
C
      IF (REAL(CORR) .EQ. 0.0 .AND. AIMAG(CORR) .EQ. 0.0) CALL XERMSG
     +   ('SLATEC', 'CLNGAM', 'Z IS A NEGATIVE INTEGER', 3, 2)
      CORR = -CMPLX (LOG(ABS(CORR)), ARGSUM)
C
C USE STIRLING-S APPROXIMATION FOR LARGE Z.
C
 50   CLNGAM = SQ2PIL + (Z-0.5)*LOG(Z) - Z + CORR + C9LGMC(Z)
      RETURN
C
      END