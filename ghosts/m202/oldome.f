*deck %W%  %G%
      SUBROUTINE OLDOME(MAXAP3,LENFOR,MOLFOR,LENFWG,FWG,ISYMM,
     $                 NOP,TRANS,NPERM,MAXOP,NATOMS,IPRINT,IDUMP,
     $                 TRVEC,RotMat,NGRP,IAN,ITRANS,TOANG,C,
     $                 JPRINT,NEQATM,NEQATL,PrtSym)
C***Begin prologue     Omega
C***Date written       850601  yymmdd
C***Revision date      yymmdd  yymmdd
C***Keywords           Symmetry
C***Author             Martin, Richard (LANL)
C***Source             %W%   %G%
C***Purpose            The final routine in the symmetry package.
C                      Communicates symmetry information to the RWF.
C***Description
C     Call Omega(MaxAP3,LenFor,MolFor,LenFWG,FWG,ISymm,
C                NOp,Trans,NPerm,MaxOp,NAtoms,IPrint,IDump,
C                TrVec,RotMat,NGrp,IAN,ITrans,ToAng,C,
C                JPrint,NEqAtm,NEqAtL,PrtSym)
C***References
C***Routines called    IOSys(IO), NoOnes(M202), NumDOF(M202), Numer(Symm),
C                      DeOrnt(M202), SubGrp(M202),
C                      CorPr1(Util)
C***End prologue       Omega
C***Begin prologue     OutRep
      Implicit Real*8(A-H,O-Z)
C
C      OMEGA IS THE FINAL ROUTINE IN THE SYMMETRY PACKAGE.
C     IT PLACES THE CHARACTER STRINGS MOLFOR AND FWG ON R/W FILE
C     SO THAT THEY MAY BE READ BY THE ARCHIVER.
C     IT MONITORS THE FRAMEWORK GROUP DURING OPTIMIZATIONS AND TURNS
C     SYMMETRY OFF IF THERE IS ANY CHANGE FROM POINT TO POINT.
C     USING THE TRANSLATION VECTOR AND ROTATION MATRIX PROVIDED BY
C     PTGRP, IT REORINTS THE COORDINATES IN BLANK COMMON (ARRAY C
C     IN THE CALLING ARGUMENTS).
C     IT EDITS THE LIST OF OPERATIONS SO THAT ONLY CARTESIAN TWO-FOLD
C     OPERATIONS REMAIN.  IT ORDERS THE OPERATIONS WRITTEN TO THE R/W
C     FILE SUCH THAT THE FIRST NOP1 DEFINE THE LARGEST CONCISE ABELIAN
C     SUBGROUP WHILE THE FIRST NOP2 (NOP2 >= NOP1) DEFINE THE LONGEST
C     ABELIAN SUBGROUP.
C
      INTEGER FWG(1), CIN(132)
      Logical Error,PrtSym
      character*3 answer
      Common/IO/Inp,IOut
      COMMON/TOL/ TOLER,TOL2
      COMMON/SYMINF/NOP1,NOP2,JTRANS(3,8),T(3,3)
      Dimension MOLFOR(1), TRANS(3,3,MAXOP), NPERM(MAXAP3,MAXOP),
     $          RotMat(3,3), TRVEC(1), NGRP(1), IAN(1),
     $       ITRANS(3,MAXOP),JPRINT(NAtoms),NEQATM(NAtoms,8),NEQATL(1),
     $          ISBGRP(3), IFLAG(8), C(3,1)
C     Note that NEQATM and NEQATL are equivalenced in the calling routine.
      Data Zero/0.D0/, Fuzzy/1.D-06/, IBlnk/1H /,
     $     IOBrak/1h[/, IHC/1hC/, IHI/1hI/
 1000 FORMAT(30H OMEGA--  SYMMETRY TURNED OFF.)
 1010 FORMAT(34H  NON-TWO-FOLD OPERATIONS DELETED./)
 1020 FORMAT(1X,44I3)
 1030 FORMAT(1X,(/1X,I2,1H:,3X,3I4))
 1040 FORMAT(29H OMEGA-- NO USEABLE SYMMETRY.)
 1050 FORMAT(24H FINAL SET OF OPERATIONS/)
 1060 Format(1X,'Standard orientation:')
 1100 FORMAT(31H PREVIOUS ROTATION MATRIX USED.)
 1110 FORMAT(57h0OMEGA--  *** WARNING *** CHANGE IN POINT GROUP DETECTED
     $./)
 1120 FORMAT(17H ROTATION MATRIX:)
 1130 FORMAT(1X,22X,3F12.6)
 1140 Format(1X,'Symmetry:')
 1150 Format(5X,'Stoichiometry       ',100A1)
 1160 Format(5X 'Framework group     ',100A1)
 1180 Format(5X,'Degrees of freedom',I3)
 1190 Format(5X,'Full point group                  ',3A1,5X,'NOp',I3)
 1200 Format(5X,'Largest Abelian subgroup          ',3A1,5X,'NOp',I3)
 1210 Format(5X,'Largest concise Abelian subgroup  ',3A1,5X,'NOp',I3)
 1220 FORMAT(87h *** WARNING *** TROUBLE REMOVING REDUNDANT OPERATIONS I
     $N OMEGA.   SYMMETRY TURNED OFF.)
 1230 FORMAT(44H *** WARNING *** UNABLE TO IDENTIFY SUBGROUP)
C
C     TREAD THE STOICHIOMETRIC FORMULA AND FRAMEWORK GROUP IF THE EXIST
C     FROM A PREVIOUS CALC'N.
C
C     CALL RTRACE(6HOMEGA ,1)
      LMOUT = intowp(9) + 26
      Call IOSys('does fwgrp exist on rwf',0,0,0,Answer)
      If(Answer.eq.'yes') then
         Call iosys('read real fwgrp from rwf',-1,CIn,0,' ')
         LRWFG=132
      Else
         LRWFWG=0
      EndIf
C
C      Store and print the current information.  Ones are stripped off
C     OF STOICHIOMETRIC FORMULAS BEFORE PRINTING.  NOTE THAT LENFOR,
C     MOLFOR, LENFWG, AND FWG ARE EQUIVALENCED TO A 66 DOUBLE-WORD
C     ARRAY IN SYMM.
C
      If(PrtSym) Write(IOut,1140)
      LEN = INTOWP(132)
      Call iosys('write real fwgrp on rwf',Len,MolFor,0,' ')
      CALL NOONES(LENFOR,MOLFOR,LENPR,JPRINT)
      If(PrtSym) WRITE(IOUT,1150) (JPRINT(I),I=1,LENPR)
      CALL NOONES(LENFWG,FWG,LENPR,JPRINT)
      If(PrtSym) WRITE(IOUT,1160) (JPRINT(J),J=1,LENPR)
      ISBGRP(1) = JPRINT(1)
      ISBGRP(2) = JPRINT(2)
      ISBGRP(3) = JPRINT(3)
      If(JPrint(3).eq.IOBrak) ISBGRP(3) = IBLNK
      If(JPrint(2).eq.IOBrak) ISBGRP(2) = IBLNK
      If(JPrint(2).eq.IOBrak) ISBGRP(3) = IBLNK
C
C     DETERMINE THE NUMBER OF DEGREES OF FREEDOM AVAILABLE TO THE
C     MOLECULE WITHIN THIS FRAMEWORK GROUP AND PRINT IT.
C
      NDOF = NUMDOF(FWG,NATOMS)
      If(PrtSym) then
         WRITE (IOUT,1180) NDOF
         WRITE (IOUT,1190) (ISBGRP(I),I=1,3),NOP
      EndIf
C
C     IF THE NOSYM FLAG IS SET, DO NO MORE.  IF IT IS SET TO 2 THEN
C     READ THE OLD TRANSFORMATION MATRIX AND ROTATE THE MOLECULE
C     USING IT.
C
      Call iosys('read integer nosym from rwf',1,NoSym,0,' ')
      IF (NOSYM .NE. 1) GO TO 30
      CALL DEORNT(RotMat,0)
      Goto 998
   30 IF (NOSYM .NE. 2) GOTO 45
      Call iosys('read integer "symmetry information" from rwf',
     $     -1,NOp1,0,' ')
      WRITE(IOUT,1100)
      NSYMOP = 0
      GOTO 480
C
C     COMPARE THE CURRENT FRAMEWORK GROUP WITH THE PREVIOUS ONE
C     (IF ANY).  IF THEY'RE NOT THE SAME SET NOSYM TO 2 AND USE THE
C     ROTATION MATRIX ON THE R/W FILES FOR THE REST OF THE JOB.
C
   45 IF (LRWFWG.eq.0) GOTO 57
      DO 50 I=1,100
      IF (CIN(I+31) .NE. FWG(I)) GOTO 55
   50 CONTINUE
      GOTO 57
   55 Call iosys('write integer nosym to rwf',1,2,0,' ')
      WRITE(IOUT,1110)
      WRITE(IOUT,1000)
      Call iosys('read integer "symmetry information" from rwf',
     $     -1,NOp1,0,' ')
      NSYMOP = 0
      GOTO 480
C
C     REMOVE NON-TWO-FOLD OPERATIONS.  IN OTHER WORDS, REMOVE THOSE
C     OPERATIONS WHERE THE TRANSFORMATION MATRIX DOES NOT CONSIST OF
C     PLUS OR MINUS ONE ON THE DIAGONAL AND ZERO OFF THE DIAGONAL.
C     THE 8 OPERATIONS WHICH CAN SURVIVE ARE E, C2(X), C2(Y), C2(Z),
C     SIGMA(XY), SIGMA(XZ), SIGMA(YZ), AND I.  THE NOP2 OPERATIONS
C     WHICH CAN SURVIVE DEFINE THE LARGES ABELIAN SUBGROUP OF THE
C     POINT GROUP.
C
C     FIRST STORE THE DIAGONAL ELEMENTS OF THE 3X3 MATRICES IN AN
C     INTEGER ARRAY, ITRANS.
C
   57 DO 80 IOP=1,NOP
      DO 80 IXYZ=1,3
      TEMP=TRANS(IXYZ,IXYZ,IOP)+Sign(FUZZY,TRANS(IXYZ,IXYZ,IOP))
      ITrans(IXYZ,IOP)=Int(Temp)
   80 CONTINUE
C
C     NOW MOVE ONLY THE TWO FOLD OPERATIONS FROM ITRANS TO JTRANS
C     AND FROM NPERM TO NEQATM.
C
      NSYMOP = 0
      DO 140 IOP=1,NOP
      JTST = IABS(ITRANS(1,IOP)) +
     $       IABS(ITRANS(2,IOP)) +
     $       IABS(ITRANS(3,IOP))
      IF (JTST .NE. 3) GOTO 140
      NSYMOP = NSYMOP + 1
      DO 100 IXYZ=1,3
      JTRANS(IXYZ,NSYMOP) = ITRANS(IXYZ,IOP)
  100 CONTINUE
      DO 120 IAT=1,NATOMS
      NEQATM(IAT,NSYMOP) = NPERM(IAT,IOP)
  120 CONTINUE
  140 CONTINUE
C
C     DEBUG PRINTING OF REMAINING OPERATIONS.
C
      IF (IPRINT.eq.0) GOTO 200
      WRITE(IOUT,1010)
      WRITE(IOUT,1020) (I,I=1,NSYMOP)
      WRITE(IOUT,1020)
      DO 160 IAT=1,NATOMS
      WRITE(IOUT,1020) (NEQATM(IAT,IOP),IOP=1,NSYMOP)
  160 CONTINUE
      DO 180 IOP=1,NSYMOP
      WRITE(IOUT,1030) IOP,(JTRANS(IXYZ,IOP),IXYZ=1,3)
  180 CONTINUE
C
  200 CONTINUE
      CALL SUBGRP (NSYMOP,JTRANS,ISBGRP,ERROR)
      IF (.NOT.ERROR) GOTO 205
         Call iosys('write integer nosym to rwf',1,1,0,' ')
         WRITE (IOUT,1230)
         WRITE (IOUT,1000)
         CALL DEORNT(RotMat,0)
         Goto 998
  205 CONTINUE
      If(PrtSym) WRITE (IOUT,1200) (ISBGRP(I),I=1,3),NSYMOP
      NOP2 = NSYMOP
C
C     EDIT THE LIST OF OPERATIONS SO THAT EACH PERMUTATION OF ATOMS
C     OCCURS ONLY ONCE IN THE FIRST NOP1 ENTRIES ON THE R/W FILE.  THIS
C     SET DEFINES THE LARGEST CONCISE ABELIAN SUBGROUP.  THE REDUNDANT
C     OPERATIONS WHICH ARE A PART OF THE LARGEST ABELIAN SUBGROUP ARE
C     MOVED TO THE END OF THE LIST.
C
C     FIRST MARK THE REDUNDANT OPERATIONS.
C
      NOP1 = NOP2
      IF (NSYMOP .LE. 1) GOTO 380
      DO 210 I=1,8
  210    IFLAG(I) = 0
      I2 = NSYMOP - 1
      DO 260 IOP=1,I2
      IF (IFLAG(IOP).eq.-99) GOTO 260
      J1 = IOP + 1
      DO 240 JOP=J1,NSYMOP
      DO 220 IAT=1,NATOMS
      IF (NEQATM(IAT,JOP) .NE. NEQATM(IAT,IOP)) GOTO 240
  220 CONTINUE
      IFLAG(JOP) = -99
      NOP1 = NOP1 - 1
  240 CONTINUE
  260 CONTINUE
C
C     NOW REMOVE THE REDUNDANT OPERATIONS.
C
      IF (NOP2 .LE. 2) GOTO 380
      JOPSAV = NOP1 + 1
      DO 360 IOP=1,NOP1
         IF (IFLAG(IOP) .NE. -99) GOTO 360
         JOPLO = JOPSAV
         DO 340 JOP=JOPLO,NOP2
            IF (IFLAG(JOP).eq.-99) GOTO 340
            IFLAG(JOP) = -99
            IFLAG(IOP) = 0
            JOPSAV = JOP + 1
            DO 315 IAT=1,NATOMS
               IJKL = NEQATM(IAT,IOP)
               NEQATM(IAT,IOP) = NEQATM(IAT,JOP)
               NEQATM(IAT,JOP) = IJKL
315            CONTINUE
            DO 325 IXYZ=1,3
               IJKL = JTRANS(IXYZ,IOP)
               JTRANS(IXYZ,IOP) = JTRANS(IXYZ,JOP)
               JTRANS(IXYZ,JOP) = IJKL
  325          CONTINUE
            GOTO 360
  340    CONTINUE
         WRITE (IOUT,1220)
         Call iosys('write integer nosym to rwf',1,1,0,' ')
  360 CONTINUE
C
  380 CONTINUE
      CALL SUBGRP (NOP1,JTRANS,ISBGRP,ERROR)
      IF (.NOT. ERROR) GOTO 385
         Call iosys('write integer nosym to rwf',1,1,0,' ')
         WRITE (IOUT,1230)
         WRITE (IOUT,1000)
         CALL DEORNT(RotMat,0)
         Goto 998
  385 CONTINUE
      If(PrtSym) WRITE (IOUT,1210) (ISBGRP(I),I=1,3),NOP1
C
C     Debug printing of the final operations.
C
      If(IPrint.eq.0) goto 440
      WRITE(IOUT,1050)
      WRITE(IOUT,1020) (I,I=1,NOP2)
      WRITE(IOUT,1020)
      DO 410 IAT=1,NATOMS
      WRITE(IOUT,1020) (NEQATM(IAT,IOP),IOP=1,NOP2)
  410 CONTINUE
      DO 420 IOP=1,NOP2
      WRITE(IOUT,1030) IOP,(JTRANS(IXYZ,IOP),IXYZ=1,3)
  420 CONTINUE
C
C     COMPRESS NEQATM INTO LINEAR FORM AND WRIT IT OUT.
  440 K=0
      DO 470 J=1,NOP2
         DO 470 I=1,NATOMS
            K=K+1
  470       NEQATL(K)=NEQATM(I,J)
      LEN=NOP2*NATOMS
      Call iosys('write integer neqatm on rwf',Len,NeqAtL,0,' ')
C
C     Transform the coordinates using the specified translation
C     and rotation.  Skip the rotation if the point group is CI.
C
  480 Do 481 I=1,3
         Do 481 J=1,3
  481 T(J,I)=RotMat(J,I)
      Call iosys('write integer "symmetry information" on rwf',
     $     LMOut,NOp1,0,' ')
      Do 500 IAt = 1, NAtoms
          Do 500 IXYZ = 1, 3
  500         C(IXYZ,IAT) = C(IXYZ,IAT) + TRVEC(IXYZ)
      If(NGrp(1).eq.IHC.and.NGrp(2).eq.IHI) GOTO 540
      Do 520 IAT=1,NATOMS
          CX = T(1,1)*C(1,IAT) + T(1,2)*C(2,IAT) + T(1,3)*C(3,IAT)
          CY = T(2,1)*C(1,IAT) + T(2,2)*C(2,IAT) + T(2,3)*C(3,IAT)
          CZ = T(3,1)*C(1,IAT) + T(3,2)*C(2,IAT) + T(3,3)*C(3,IAT)
          If(Abs(CX).LT.Toler) CX = ZERO
          If(Abs(CY).LT.Toler) CY = ZERO
          If(Abs(CZ).LT.Toler) CZ = ZERO
          C(1,IAT) = CX
          C(2,IAT) = CY
  520     C(3,IAT) = CZ
C
C     PRINT THE COORDINATES.
C
  540 Continue
      If(PrtSym) then
         Write(IOut,1060)
         Call CorPr1(IOut,NAtoms,IAN,C,ToAng)
      EndIf
C
C     PERHAPS PRINT THE ROTATION MATRIX, THEN RETURN.
C
      If(IPrint.eq.0) goto 998
      WRITE(IOUT,1120)
      DO 560 I=1,3
  560     WRITE(IOUT,1130) (T(I,J),J=1,3)
  998 Return
      End
