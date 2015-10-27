*deck %W%  %G%
      SUBROUTINE SUBGRP (NOP,ITR,ISBGRP,ERROR)
C***Begin prologue     SubGrp
C***Date written       850601  yymmdd
C***Revision date      yymmdd  yymmdd
C***Keywords           Symmetry
C***Author             Gauss82
C***Source             %W%   %G%
C***Purpose            Determines the highest Abelian subgroup within a
C                      point group.
C***Description
C     Call SubGrp(NOp,ITr,ISbGrp,Error)
C
C     INPUT
C           NOP ...... NUMBER OF TWO-FOLD SYMMETRY OPERATIONS.
C           ITR ...... INTEGER FORM OF THE 3X3 TRANSFORMATION
C                      MATRICES, ONE FOR EACH OF THE NOP OPERATIONS.
C     OUTPUT
C           ISBGRP ... THE SCHONFLIES SYMBOL OF THE ABELIAN SUBGROUP
C                      DEFINED BY THE OPERATIONS INPUT.  THIS IS
C                      RETURNED AS AN UNPACKED STRING OF LENGTH 3.
C           ERROR .... SET TRUE IF THE SUBGROUP CAN NOT BE IDENTIFIED.
C                      NOTE THAT ONLY THE TYPES OF OPERATIONS ARE
C                      EXAMINED; THEIR SPATIAL RELATIONSHIPS ARE NOT.
C                      THUS ANY GROUP CONSISTING OF THE IDENTITY, TWO
C                      PLANES AND AN AXIS IS CLASSIFIED AS C2V, WHETHER
C                      OR NOT THE PLANES BOTH CONTAIN THE C2 AXIS.
C
C***References
C***Routines called    (None)
C***End prologue       SubGrp
C
      IMPLICIT INTEGER (A-Z)
      LOGICAL ERROR
C
      DIMENSION ISBGRP(1), ITR(3,8)
      DIMENSION GROUPS (3,9)
C
      DATA GROUPS /1HC, 1H1, 1H ,
     $             1HC, 1HI, 1H ,
     $             1HC, 1H2, 1H ,
     $             1HC, 1HS, 1H ,
     $             1HC, 1H2, 1HV,
     $             1HD, 1H2, 1H ,
     $             1HC, 1H2, 1HH,
     $             1HD, 1H2, 1HH,
     $             1H?, 1H?, 1H?/
C
C
C      THE EIGHT POSSIBILITIES ARE  C1, CI, C2, CS, C2V, D2, C2H, D2H
C      THE POSSIBLE VALUES OF NOP ARE 1, 2, 4, AND 8
C      BRANCH ON NOP
C
      IF (NOP.LT.1 .OR. NOP.GT.8) GOTO 1000
      NGROUP = 9
      GOTO (100,200,1000,400,1000,1000,1000,800), NOP
C
C      ONE OPERATION   THE SUBGROUP IS C1
C
  100 CONTINUE
      NGROUP = 1
      GOTO 1000
C
C      TWO OPERATIONS   THE POSSIBILITIES ARE CI, CS, C1
C
  200 CONTINUE
      ITST = ITR(1,2) + ITR(2,2) + ITR(3,2)
      IF (ITST .EQ. -3) NGROUP = 2
      IF (ITST .EQ. -1) NGROUP = 3
      IF (ITST .EQ. +1) NGROUP = 4
      GOTO 1000
C
C      FOUR OPERATIONS   THE POSSIBILITIES ARE C2H, C2V, D2
C
  400 CONTINUE
      ITST = ITR(1,2) + ITR(2,2) + ITR(3,2)
      JTST = ITR(1,3) + ITR(2,3) + ITR(3,3)
      KTST = ITR(1,4) + ITR(2,4) + ITR(3,4)
      LTST = ITST + JTST + KTST
      IF (LTST .EQ. 1) NGROUP = 5
      IF (LTST.EQ.-3 .AND. ITST.EQ.JTST .AND. JTST.EQ.KTST) NGROUP = 6
      IF (LTST.EQ.-3 .AND. NGROUP.EQ.9) NGROUP = 7
      GOTO 1000
C
C      EIGHT OPERATIONS   THE SUBGROUP IS D2H
C
  800 CONTINUE
      NGROUP = 8
      GOTO 1000
C
C
 1000 CONTINUE
      ISBGRP(1) = GROUPS (1,NGROUP)
      ISBGRP(2) = GROUPS (2,NGROUP)
      ISBGRP(3) = GROUPS (3,NGROUP)
      ERROR = NGROUP .EQ. 9
C
      RETURN
      END
