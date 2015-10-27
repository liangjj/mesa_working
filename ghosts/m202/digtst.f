*deck %W%  %G%
      SUBROUTINE DIGTST(FWG,IPOS,DIGIT,IVAL)
C***Begin prologue     DigTst
C***Date written       850601  yymmdd
C***Revision date      yymmdd  yymmdd
C***Keywords           Framework group
C***Author             Gauss82
C***Source             %W%   %G%
C***Purpose            Tests the current position in the FWG for a number.
C***Description
C     Call DigTst(FWG,IPos,Digit,IVal)
C***References
C***Routines called    (None)
C***End prologue       DigTst
C
C      TEST THE CURRENT POSITION IN FWG FOR A NUMBER.
C      INPUT   FWG,   THE FRAMEWORK GROUP STRING.
C              IPOS,  THE CHARACTER IN FWG TO TEST.  IF A DIGIT IS
C                     FOUND THEN IPOS IS LEFT POINTING AT THE FIRST
C                     CHARACTER PASSED THE NUMBER (INDEPENDENT OF HOW
C                     LONG IT IS.)
C      OUTPUT  DIGIT, A LOGICAL VARIABLE SET TRUE IF FWG(IPOS) IS A
C                     DIGIT.
C              IVAL,  THE VALUE OF THE DIGIT IN THE RANGE 0- .
C
      Implicit Integer(A-Z)
      Dimension FWG(*), Num(10)
      Logical Digit
C
      DATA NUM/1H0, 1H1, 1H2, 1H3, 1H4, 1H5, 1H6, 1H7, 1H8, 1H9/
C
C
      IVAL = 0
      DIGIT = .FALSE.
C
   10 CONTINUE
      DO 20 I=1,10
         IF (FWG(IPOS) .NE. NUM(I)) GOTO 20
            DIGIT = .TRUE.
            IPOS = IPOS + 1
            IVAL = 10*IVAL + I - 1
            GOTO 10
   20    CONTINUE
      RETURN
C
      END
