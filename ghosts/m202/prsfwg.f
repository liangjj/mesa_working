*deck %W%  %G%
      SUBROUTINE PRSFWG(FWG,IPOS,ICHAR,JCHAR,NATSS)
C***Begin prologue     PrsFWG
C***Date written       850601  yymmdd
C***Revision date      yymmdd  yymmdd
C***Keywords           Symmetry
C***Author             Gauss82
C***Source             %W%   %G%
C***Purpose            Parses the framework group string.
C***Description
C     Call PrsFWG(FWG,IPos,IChar,JChar,NAtSS)
C
C      PARSE FRAMEWORK GROUP STRING ONE SUBSPACE AT A TIME FOR FUNCTION
C      NUMDOF.
C      INPUT   FWG,   UNPACKED FRAMEWORK GROUP STRING.  THE FIRST
C                     FOUR CHARACTERS SHOULD BE THE POINT GROUP WITH
C                     CHARACTER FOUR USED EVEN FOR N<10.
C              IPOS,  THE CURRENT POSITION IN FWG.  THIS SHOULD EITHER
C                     BE ZERO (FOR THE FIRST CALL) OR THE FIRST
C                     CHARACTER OF A SUBSPACE UPON ENTRY TO PRSFWG.
C                     UPON EXIT IT IS SET TO THE FIRST CHARACTER OF THE
C                     NEXT SUBSPACE.
C      OUTPUT  ICHAR  FIRST CHARACTER IN SYMMETRIC SUBSPACE IDENTIFIER
C                     ('O','C','S','X'), SET TO 'E' UPON AN ERROR.
C              JCHAR  UNDEFINED FOR ICHAR = 'O' OR 'E',
C                     ORDER OF SUBSPACE FOR ICHAR = 'C',
C                     ORDER OF POINT GROUP FOR ICHAR = 'X',
C                     'H' FOR ICHAR = S IN SIGMAH, ETC.
C              NATSS  NUMBER OF ATOMS IN THE SUBSPACE.
C***References
C***Routines called    DigTst(M202), Numer(Symm)
C***End prologue       PrsFWG
      Implicit Real*8(A-H,O-Z)
C
      INTEGER FWG(1), RBrack, RParen
      LOGICAL DIGIT
      Common/IO/Inp,IOut
      Data IBlank/1h /, LBrack/1h[/, RBrack/1h]/, LParen/1h(/,
     $     RParen/1h)/, IComma/1h,/, IStar/1h*/,
     $     IHC/1hC/, IHD/1hD/, IHE/1hE/, IHH/1hH/, IHO/1hO/,
     $     IHS/1hS/, IHT/1hT/, IHX/1hX/
 1010 FORMAT(45H PRSFWG-- UNABLE TO FIND LEFT SQUARE BRACKET.)
 1020 FORMAT(51H PRSFWG-- NON-DIGIT FOLLOWS "C" IN LINEAR SUBSPACE.)
 1030 FORMAT(41H PRSFWG-- UNRECOGNIZED SUBSPACE, ICHAR  ",A4,1H")
 1040 FORMAT(45H PRSFWG-- UNRECOGNIZED POINT GROUP, ICHAR="X")
 1050 FORMAT(37H PRSFWG-- LEFT PARENTHESIS NOT FOUND.)
 1060 FORMAT(38H PRSFWG-- RIGHT PARENTHESIS NOT FOUND.)
C
C
C                                       IPOS=0 UPON THE INITIAL ENTRY.
      NATSS = 0
      IF (IPOS .NE. 0) GOTO 60
         DO 20 I=1,100
            IF (FWG(I) .NE. LBrack) GOTO 20
               IPOS = I + 1
               GOTO 60
   20       CONTINUE
         WRITE (IOUT,1010)
         ICHAR = IHE
         RETURN
C                                       IPOS IS POSITIONED AT START OF
C                                       A SUBSPACE.  CHECK FOR A DIGIT.
   60 CONTINUE
      MULT = 1
      CALL DIGTST(FWG,IPOS,DIGIT,IVAL)
      IF (DIGIT) MULT = IVAL
      ICHAR = FWG(IPOS)
      IPOS = IPOS + 1
C                                       DETERMINE JCHAR.
      IF (ICHAR .EQ. RBrack) GOTO 1000
      IF (ICHAR .EQ. IHO) GOTO 200
      IF (ICHAR .NE. IHC) GOTO 80
         CALL DIGTST(FWG,IPOS,DIGIT,IVAL)
         IF (.NOT. DIGIT .AND. FWG(IPOS).NE.IStar) GOTO 70
            JCHAR = IVAL
            GOTO 200
   70    ICHAR = IHE
         WRITE (IOUT,1020)
         RETURN
   80 CONTINUE
      IF (ICHAR .NE. IHS) GOTO 100
         JCHAR = FWG(IPOS+1)
         IPOS = IPOS + 1
         IF (JCHAR .NE. LParen) GOTO 200
            JCHAR = IBlank
            IPOS = IPOS - 1
            GOTO 200
  100 CONTINUE
      IF (ICHAR .NE. IHX) GOTO 160
         NPRIN = NUMER(FWG)
         IF (FWG(1) .NE. IHC) GOTO 120
            JCHAR = 2 * NPRIN
            IF (FWG(4) .EQ. IBlank) JCHAR = NPRIN
            IF (NPRIN .EQ. 0) JCHAR = 2
            GOTO 200
  120    IF (FWG(1) .NE. IHD) GOTO 125
            JCHAR = 4 * NPRIN
            IF (FWG(4) .EQ. IBlank) JCHAR = 2 * NPRIN
            GOTO 200
  125    IF (FWG(1) .NE. IHS) GOTO 130
            JCHAR = NPRIN
            GOTO 200
  130    IF (FWG(1) .NE. IHT) GOTO 140
            JCHAR = 12
            IF (FWG(2) .EQ. IHD) JCHAR = 24
            GOTO 200
  140    IF (FWG(1) .NE. IHO) GOTO 150
            JCHAR = 24
            IF (FWG(2) .EQ. IHH) JCHAR = 48
            GOTO 200
  150    WRITE (IOUT,1040)
         ICHAR = IHE
         RETURN
  160 CONTINUE
      WRITE (IOUT,1030) ICHAR
      ICHAR = IHE
      RETURN
C                                         FIND LEFT PARENTHESIS.
  200 CONTINUE
         IF (IPOS .LE. 100) GOTO 210
            ICHAR = IHE
            WRITE (IOUT,1050)
            RETURN
  210    IF (FWG(IPOS) .EQ. LParen) GOTO 220
         IPOS = IPOS + 1
         GOTO 200
C                                         SUM ATOMS UNTIL RIGHT
C                                         PARENTHESIS IS ENCOUNTERED.
  220 CONTINUE
         IPOS = IPOS + 1
         IF (IPOS .LE. 100) GOTO 225
            ICHAR = IHE
            WRITE (IOUT,1060)
            RETURN
  225    IF (FWG(IPOS) .EQ. RParen) GOTO 240
  230    CALL DIGTST(FWG,IPOS,DIGIT,IVAL)
         IF (.NOT. DIGIT) GOTO 220
            IPOS = IPOS - 1
            NATSS = NATSS + IVAL*MULT
            GOTO 220
C                                         MOVE IPOS TO THE NEXT SUBSPACE
  240 CONTINUE
      IPOS = IPOS + 1
      IF (FWG(IPOS) .EQ. IComma) IPOS = IPOS + 1
      GOTO 1000
C                                         NORMAL EXIT.
 1000 CONTINUE
      RETURN
      END
