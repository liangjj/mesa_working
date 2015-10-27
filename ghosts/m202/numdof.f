*deck %W%  %G%
      Integer Function NumDOF(Frame,NAtoms)
C***Begin prologue     NumDOF
C***Date written       850601  yymmdd
C***Revision date      yymmdd  yymmdd
C***Keywords           Symmetry, Degrees of freedom
C***Author             Martin, Richard (LANL)
C***Source             %W%   %G%
C***Purpose            Determines the number of degrees of freedom within
C                      a given symmetry constraint.
C***Description
C     NumDOF is an integer function used as:
C       Number=NumDOF(Frame,NAtoms)
C         Frame   The unpacked framework group string.
C         NAtoms  The number of atoms.
C         NumDOF  The number of degrees of freedom within the given
C                 symmetry constraint.  A zero is returned and a
C                 message printed if an error is detected.
C
C                 To determine the number of equivalent sets of atoms in
C                 each symmetric subspace you divide the total number
C                 of atoms in the subspace by the number of atoms in the
C                 subspace that are equivalent under the operations of the
C                 group.  The divisor is dependent on the subspace and
C                 the point group.
C***References         Pople,Satay,and Halevi, Israel J. Chem.(1979)
C***Routines called    PrsFWG(M202)
C***End prologue       NumDOF
      Implicit Real*8(A-H,O-Z)
      INTEGER FRAME(1), FWG(100)
      Common/IO/Inp,IOut
      Data IHK/1hK/, IHI/1hI/, IHT/1hT/, IBlank/1h /, IHD/1hD/,
     $     ICBrak/1h]/, IHO/1hO/, IHC/1hC/, IHH/1hH/, IHS/1hS/,
     $     IHX/1hX/, IHV/1hV/
 1010 FORMAT(49H NUMDOF-- NOT CODED TO HANDLE GROUPS TH, I, OR IH)
 1020 FORMAT(50H NUMDOF-- UNRECOGNIZED SYMMETRIC SUBSPACE, ICHAR= ,
     $  A4,1H")
 1030 FORMAT(15H NUMDOF-- NAT= ,I5,9H NATOMS= ,I5)
C
      DO 10 I=1,100
         FWG(I) = FRAME(I)
   10    CONTINUE
C
C                                   NPRIN  ORDER OF THE PRINCIPAL AXIS.
C                                   IPOS   CURRENT POSITION IN FWG.
C                                   NAT    NUMBER OF ATOMS ACCORDING TO
C                                          FWG, TO BE COMPARED WITH
C                                          NATOMS.
C
      NUMDOF = 0
      IPOS   = 0
      NAT    = 0
      NPRIN  = NUMER(FWG)
C                                    TEST FOR ATOMS
      IF (FWG(1) .EQ. IHK) RETURN
C                                    TEST FOR CUBIC GROUPS
      IF (FWG(1).NE.IHI .AND. FWG(2).NE.IHT) GOTO 20
      IF ((FWG(1).EQ.IHT .AND. FWG(2).EQ.IBlank) .OR.
     $    (FWG(1).EQ.IHT .AND. FWG(2).EQ.IHD)) GOTO 20
         WRITE (IOUT,1010)
         NUMDOF = -1
         GOTO 200
C                                    FWG IS PARSED BY SUBROUTINE PRSFWG.
C                                    IT FINDS THE NEXT SYMMETRIC
C                                    SUBSPACE AND RETURNS ITS TYPE AND
C                                    THE NUMBER OF ATOMS IS CONTAINS.
C                                    LOOP OVER SYMMETRIC SUBSPACES
   20 CONTINUE
         CALL PRSFWG(FWG,IPOS,ICHAR,JCHAR,NATSS)
         IF (ICHAR .EQ. ICBrak) GOTO 100
         NAT = NAT + NATSS
C                                    POINT SUBSPACE.
         IF (ICHAR .EQ. IHO) GOTO 20
C                                    LINEAR SUBSPACE.
C                                    JCHAR IS THE ORDER OF THE SUBSPACE.
C                                        CASE 1- SN, D*H, CNH, DN
C                                        CASE 2- C2 IN D GROUPS
C                                        CASE 3- CN, CNV, C*V
C                                        CASE 4- C3 IN T, TD
C                                        CASE 5- C4 IN O, OH
C                                        CASE 6- C2 IN T, TD, O, OH
C                                                C3 IN O, OH
         IF (ICHAR .NE. IHC) GOTO 40
            NDOF = NATSS / 2
            IF (JCHAR.EQ.2 .AND. FWG(1).EQ.IHD) NDOF = NATSS / NPRIN
            IF (FWG(1).EQ.IHC .AND. FWG(4).NE.IHH) NDOF = NATSS
            IF (FWG(1).EQ.IHT .AND. JCHAR.EQ.3) NDOF = NATSS / 4
            IF (FWG(1).EQ.IHO .AND. JCHAR.EQ.4) NDOF = NATSS / 6
            IF ((FWG(1).EQ.IHO .AND. (JCHAR.EQ.2.OR.JCHAR.EQ.3)) .OR.
     $          (FWG(1).EQ.IHT.AND.JCHAR.EQ.2)) NDOF = NATSS / 8
            NUMDOF = NUMDOF + NDOF
            GOTO 20
C                                  PLANAR SUBSPACES.
C                                  JCHAR CONTAINS 'H' FOR SIGMAH, ETC.
C                                    CASE 1- SIGMAV/SIGMAD IN DNH, DND.
C                                            SIGMA/SIGMA'/SIGMA'' IN D2H
C                                            SIGMAH IN DNH
C                                    CASE 2- SIGMA IN CS
C                                            SIGMAH IN  CNH
C                                            SIGMAV/SIGMAD IN CNV.
C                                    CASE 3- SIGMAD IN TD, SIGMAH IN OH
C                                            SIGMAD IN OH
   40    CONTINUE
         IF (ICHAR .NE. IHS) GOTO 60
            IF (NPRIN .EQ. 0) NPRIN = 1
            NDOF = NATSS / NPRIN
            IF (FWG(1) .EQ. IHC) NDOF = 2 * NATSS / NPRIN
            IF (FWG(1).EQ.IHT .OR. FWG(1).EQ.IHO) NDOF = NATSS / 12
            NUMDOF = NUMDOF + NDOF
            GOTO 20
C                                    ASSYMETRIC SUBSPACE.
C                                    JCHAR IS THE ORDER OF THE PT GRP.
   60    CONTINUE
         IF (ICHAR .NE. IHX) GOTO 80
            NUMDOF = NUMDOF  +  3 * (NATSS/JCHAR)
            GOTO 20
C                                    UNRECOGNIZED SUBSPACE.
   80    CONTINUE
            NUMDOF = -1
            WRITE (IOUT,1020) ICHAR
            GOTO 200
C                                    END OF LOOP OVER SUBSPACES.
C                                    CHECK NAT, IF OK THEN SUBTRACT
C                                    FOR TRANSLATION/ROTATION AND
C                                    RETURN.
  100 CONTINUE
      IF (NAT .EQ. NATOMS) GOTO 120
         NUMDOF = 0
         WRITE (IOUT,1030) NAT,NATOMS
         GOTO 200
  120 CONTINUE
      IF (FWG(1).NE.IHC .AND. FWG(1).NE.IHS) GOTO 200
         NS = 0
         IF (FWG(1) .EQ. IHS) NS = 1
         IF (FWG(4) .EQ. IHH) NS = 1
         IF (FWG(4).EQ.IHV .OR. FWG(3).EQ.IHV) NS = 1
         IF (FWG(2) .EQ. IHS) NS = 3
         IF (FWG(2) .EQ. IHI) NS = 3
         IF (NS.EQ.0 .AND. NPRIN.EQ.1) NS = 6
         IF (NS .EQ. 0) NS = 2
         NUMDOF = NUMDOF - NS
         GOTO 200
C
  200 CONTINUE
      RETURN
C
      END
