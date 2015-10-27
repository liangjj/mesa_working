MODULE Extract_mod
   IMPLICIT NONE
   INTEGER, PARAMETER :: XR = 8
   INTEGER, PARAMETER :: InpLineL = 256
   REAL (KIND = XR) :: PClose = 1.0e-6_XR
   REAL (KIND = XR), PARAMETER :: PI = 3.141592653589793E0_XR
   REAL (KIND = XR), PARAMETER :: XFEVAU = 27.2113834E0_XR
   REAL (KIND = XR), PARAMETER :: XFC = 137.03600E0_XR
   REAL (KIND = XR), PARAMETER :: xfauang = 0.5291772083E0_XR
   REAL (KIND = XR), PARAMETER :: XFAUMB = 100.0_XR*xfauang*xfauang
   TYPE IntVec
      INTEGER, POINTER, DIMENSION(:) :: v
   END TYPE IntVec

CONTAINS
   CHARACTER (LEN = 10) FUNCTION Int4ToCharL(IntIn)
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::IntIn
      INTEGER :: fac, i, j, js
      INTEGER :: resid, newdig
      CHARACTER (LEN = 1), DIMENSION(0:9), PARAMETER :: dig = (/ "0", "1", "2", "3", "4", "5", "6", "7", "8", "9" /)

      Int4ToCharL = " "
      IF (IntIn == 0) THEN
         Int4ToCharL(1:1) = "0"
      ELSE
         fac = 1000000000
         resid = IntIn
         j = 0
         IF (resid < 0) THEN
            Int4ToCharL(1:1) = "-"
            j = 1
            resid = -resid
         END IF
         js = j
         DO i = 1, 10
            newdig = resid/fac
            resid = resid-newdig*fac
            fac = fac/10
            IF (newdig /= 0 .OR. j /= js) THEN
               j = j+1
               Int4ToCharL(j:j) = dig(newdig)
            END IF
         END DO
      END IF
   END FUNCTION Int4ToCharL
   INTEGER FUNCTION GetSymTyp(Label, List)
      IMPLICIT NONE
      CHARACTER (LEN = *) :: Label
      CHARACTER (LEN = *), DIMENSION(:) :: List
      INTEGER :: i

      DO i = 1, SIZE(List)
         IF (Label == List(i)) THEN
            GetSymTyp = i
            EXIT
         END IF
      END DO
      IF (i > SIZE(List)) THEN
         WRITE (UNIT = 6, FMT = "('GetSymTyp, could not find ', a5, '  in:', /, (10a5))") Label, List
         STOP
      END IF
   END FUNCTION GetSymTyp
   FUNCTION MakeList(Vals)
      IMPLICIT NONE

      CHARACTER (LEN = InpLineL) :: MakeList
      INTEGER, INTENT(IN), DIMENSION(:) :: Vals
      
      INTEGER :: i

      MakeList = "("
      DO i = 1, SIZE(Vals)
         IF (i > 1) THEN
            MakeList = TRIM(MakeList)//","
         END IF
         MakeList = TRIM(MakeList)//TRIM(Int4ToCharL(Vals(i)))
      END DO

      MakeList = TRIM(MakeList)//")"
   END FUNCTION MakeList
   SUBROUTINE ToCaps(StringIn, StringOut)
      IMPLICIT NONE
      INTEGER :: i, j
      CHARACTER (LEN = *), INTENT(IN) :: StringIn
      CHARACTER (LEN = *), INTENT(OUT) :: StringOut
      CHARACTER (LEN = 26) :: Lower = 'abcdefghijklmnopqrstuvwxyz', Upper = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
      
      StringOut = StringIn
      
      DO i = 1, MIN(LEN(StringIn), LEN(StringOut))
         DO j = 1, LEN(Lower)
            IF (StringIn(i:i) == Lower(j:j)) THEN
               StringOut(i:i) = Upper(j:j)
               EXIT
            END IF
         END DO
      END DO

   END SUBROUTINE ToCaps
   SUBROUTINE AppendFile(StartFile, AddFile, UnitTmp2, UnitTmp3)
      IMPLICIT NONE
      CHARACTER (LEN = *) :: StartFile, AddFile
      INTEGER (KIND = 4) :: UnitTmp2, UnitTmp3
      INTEGER (KIND = 4) :: ios
      CHARACTER (LEN = InpLineL) :: Line

      OPEN (UNIT = UnitTmp2, FILE = StartFile, POSITION = "APPEND")
      OPEN (UNIT = UnitTmp3, FILE = AddFile, POSITION = "REWIND")

      DO
         READ (UNIT = UnitTmp3, FMT = "(a)", IOSTAT = ios) Line
         IF (ios /= 0) THEN
            EXIT
         END IF
         WRITE (UNIT = UnitTmp2, FMT = "(a)") TRIM(Line)
      END DO

      CLOSE (UNIT = UnitTmp2)
      CLOSE (UNIT = UnitTmp3)

   END SUBROUTINE AppendFile
   SUBROUTINE AddToDRTLine(Line, NOrb, Group, NS, LineCnt, LineCntX, UnitOut)
      IMPLICIT NONE
      CHARACTER (LEN = *) :: Line
      INTEGER, INTENT(IN) :: NOrb, Group, NS, LineCntX
      INTEGER (KIND = 4), INTENT(IN) :: UnitOut
      INTEGER, INTENT(INOUT) :: LineCnt
      CHARACTER (LEN = 10) :: CNOrb, CGroup, CNS

      WRITE (UNIT = CNOrb, FMT = "(i10)") NOrb
      WRITE (UNIT = CGroup, FMT = "(i10)") Group
      WRITE (UNIT = CNS, FMT = "(i10)") NS

      Line = TRIM(LINE)//" "//TRIM(ADJUSTL(CNOrb))//"typ"//TRIM(ADJUSTL(CGroup))&
           &//";"//TRIM(ADJUSTL(CNS))
      Line = ADJUSTL(Line)
      LineCnt = LineCnt+1
      IF (LineCnt >= LineCntX) THEN
         WRITE (UNIT = UnitOut, FMT = "(a)") TRIM(Line)
         LineCnt = 0
         LIne = " "
      END IF
   END SUBROUTINE AddToDRTLine
      FUNCTION CLGAMI(ZR, ZI)
! **      USE SetUp_par
      IMPLICIT NONE

      REAL (KIND = XR) :: ZR, ZI
      REAL (KIND = XR) :: CLGAMI

      REAL (KIND = XR) :: V(2), H(2), R(2), S(2), S2(2), CL(2), S3(2)
      EQUIVALENCE (VR,V(1)), (VI,V(2)), (HR,H(1)), (HI,H(2))
      EQUIVALENCE (CLR,CL(1)), (CLI,CL(2)), (SI,S(2)), (SR,S(1))
      REAL (KIND = XR), DIMENSION(2) :: CONE = (/ 1.0_XR,0.0_XR/), P5 = (/0.5_XR,0.0_XR/)
      REAL (KIND = XR), DIMENSION(10) :: B = (/ &
              +8.33333333333333E-2_XR, -2.77777777777778E-3_XR, &
              +7.93650793650794E-4_XR, -5.95238095238095E-4_XR, &
              +8.41750841750842E-4_XR, -1.91752691752692E-3_XR, &
              +6.41025641025641E-3_XR, -2.95506535947712E-2_XR, &
              +1.79644372368831E-1_XR, -1.39243221690590_XR /)
      REAL (KIND = XR) :: X, T, F, VR, VI, HR, HI, C, D, A, SR, SI
      REAL (KIND = XR) :: E, CLR, CLI
      INTEGER :: N, II, I

      X = ZR
      T = ZI
      IF(ABS(AINT(X)+ABS(X)) .LT. 1.0E-10_XR .AND.ABS(T) .LT. 1.0E-10_XR) &
            GOTO 5
      F = ABS(T)
      VR = X
      VI = F
      IF (X .LT. 0.0_XR) CALL  CMIN(V,CONE,V)
      HR = 0.0_XR
      HI = 0.0_XR
      C = VR
      IF (C .GE. 7.0_XR) GO TO 3
      N = 6 - INT(C)
      HR = VR
      HI = VI
      D = VI
      A = ATAN2(D, C)
      IF (N .EQ. 0) GO TO 2

      DO 1 I = 1, N
         C = C + 1.0_XR
         VR = C
         VI = D
         CALL  CMULX(H, V)
         A = A + ATAN2(D, C)
 1    CONTINUE

    2 HR = 0.5_XR * LOG(HR**2 + HI**2)
      HI = A
      VR = VR + 1.0_XR
    3 CALL CDIV(S, CONE, V)
      CALL CMUL(R, S, S)
      SR = B(10)
      SI = 0.0_XR

      DO 12 II = 1, 9
         I = 10 - II
         CALL CMULX(S, R)
         SR = SR + B(I)
 12   CONTINUE

      CALL CDIV(S2, S, V)
      CALL CMIN(S, S2, H)
      CALL CMIN(S, S, V)
      CALL CMIN(S2, V, P5)
      CALL CLOX(S3, V)
      CALL CMULX(S2, S3)
      CALL CADD2(S, S, S2)
      CLR = 0.918938533204673_XR + SR
      CLI = SI
      IF (X .GE. 0.0_XR) GO TO 4
      A = AINT(X) - 1.0_XR
      C = PI * (X-A)
      D = PI * F
      E = EXP(-2.0_XR*D)
      F = SIN(C)
      E = D+0.5_XR*LOG(E*F**2+0.25_XR*(1.0_XR-E)**2)
      F = ATAN2(COS(C)*TANH(D),F)-A*PI
      CLR = 1.144729885849400_XR-E-CLR
      CLI = - F - CLI
    4 IF (T .LT. 0.0_XR) CLI = - CLI
      CLGAMI = CLI
      RETURN
    5 WRITE (6, 100) X
      CLGAMI = 0.0_XR
      STOP 'Bad cal to CLOGAM'
      RETURN
  100 FORMAT(1X,'CLOGAM ... ARGUMENT IS NON-POSITIVE INTEGER = ',F20.2)
      END FUNCTION CLGAMI
! ==============================================================
      SUBROUTINE CDIV(A,B,C)
! **       USE SetUp
      IMPLICIT NONE

      REAL (KIND = XR) :: A(2),B(2),C(2)
      REAL (KIND = XR) :: R, AR, AI

      R=C(1)**2+C(2)**2
      AR=(B(1)*C(1)+B(2)*C(2))/R
      AI=(C(1)*B(2)-B(1)*C(2))/R
      A(1)=AR
      A(2)=AI

      RETURN
      END SUBROUTINE CDIV
      SUBROUTINE CLOX(A,B)
! **      USE SetUp
      IMPLICIT REAL (KIND = XR) (A-H, O-Z)
      DIMENSION A(2),B(2)
      AR=LOG(SQRT(B(1)**2+B(2)**2))
      AI=ATAN2(B(2),B(1))
      A(1)=AR
      A(2)=AI
      RETURN


   END SUBROUTINE CLOX
      SUBROUTINE CMIN(A,B,C)
! **      USE SetUp
      IMPLICIT REAL (KIND = XR) (A-H, O-Z)
      DIMENSION A(2),B(2),C(2)
      AR=B(1)-C(1)
      AI=B(2)-C(2)
      A(1)=AR
      A(2)=AI
      RETURN


   END SUBROUTINE CMIN
      SUBROUTINE CMUL(A,B,C)
! **      USE SetUp
      IMPLICIT REAL (KIND = XR) (A-H, O-Z)
      DIMENSION A(2),B(2),C(2)
      AR=B(1)*C(1)-B(2)*C(2)
      AI=B(1)*C(2)+C(1)*B(2)
      A(1)=AR
      A(2)=AI
      RETURN


   END SUBROUTINE CMUL
      SUBROUTINE CMULX(B,C)
! **      USE SetUp
      IMPLICIT REAL (KIND = XR) (A-H, O-Z)
      DIMENSION B(2),C(2)
      AR=B(1)*C(1)-B(2)*C(2)
      AI=B(1)*C(2)+C(1)*B(2)
      B(1)=AR
      B(2)=AI
      RETURN


   END SUBROUTINE CMULX
      SUBROUTINE CADD2(A,B,C)
! **       USE SetUp
      IMPLICIT REAL (KIND = XR) (A-H, O-Z)
      DIMENSION A(2),B(2),C(2)
      AR=B(1)+C(1)
      AI=B(2)+C(2)
      A(1)=AR
      A(2)=AI
      RETURN


   END SUBROUTINE CADD2

   SUBROUTINE RawPhse(tkin, LMaxIdy, Phase)

! ** ipx, LMax, var, nrdimTarg)

! **       USE SetUp_par
      IMPLICIT NONE

      COMPLEX (KIND = XR), INTENT(OUT), DIMENSION(0:LMaxIdy) :: Phase
      REAL (KIND = XR), INTENT(IN) :: tkin
      INTEGER :: LMaxIdy

! ** local variables

      COMPLEX (KIND = XR) :: PH, SIGL
      REAL (KIND = XR) :: C

      INTEGER :: LNX
      REAL (KIND = XR) :: ETA

!     COMPUTE PHASE FACTOR

      LNX = 0
      ETA = -1.0_XR/TKIN
      SIGL = CLGAMI(REAL(LNX+1, KIND = XR), ETA)
      PH = LNX * 0.5_XR * PI - SIGL
      C = SQRT(2.0_XR/PI)
      PHASE(LNX) = C * EXP((0.0_XR, 1.0_XR)*PH)

      DO LNX = 1, LMaxIdy
         SIGL = SIGL + ATAN2(ETA, REAL(LNX, KIND = XR))
         PH = LNX * 0.5_XR * PI - SIGL
         PHASE(LNX) = C * EXP((0.0_XR, 1.0_XR)*PH)
      END DO

   END SUBROUTINE RawPhse
END MODULE Extract_mod
