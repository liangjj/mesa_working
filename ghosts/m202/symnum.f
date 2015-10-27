*deck %W%  %G%
      Function Symnum(Linear,PG)
C***Begin prologue     SymNum
C***Date written       850602  yymmdd
C***Revision date      yymmdd  yymmdd
C***Keywords           Symmetry
C***Author             Gauss82
C***Source             %W%   %G%
C***Purpose            Returns the rotational symmetry number of the molecule.
C***Description
C     SymNum is a real function used as:
C       RotNum=SymNum(Linear,PTGrp)
C
C     PROVIDE THE CALLING ROUTINE WITH THE ROTATIONAL SYMMETRY NUMBER
C     OF THE MOLECULE.  A TABLE LOOK UP IS DONE BASED UPON THE POINT
C     GROUP OF THE MOLECULE.  SEE S. W. BENSON, "THERMOCHEMICAL
C     KINETICS, 2ND ED.", WILEY, NEW YORK, 1976, P49.
C     IN ADDITION, THE LOGICAL VARIABLE LINEAR IS SET TRUE IF THE
C     MOLECULAR POINT GROUP IS D*H OR C*V.  PtGrp is the point group
C     as an unpacked character string.
C***References         S.W.Benson,Thermochemical Kinetics, 2nd. Ed.,49(1976).
C***Routines called    (None)
C***End prologue       SymNum
c
      Implicit Real*8(A-H,O-Z)
C
C
      real*8 symnum
      integer PG(*)
      integer num(10)
      Logical Linear, Test
      Data One,Two,Twelve,F24/1.D0,2.D0,12.D0,24.D0/, IHSt/1h*/
      Data Num/1H0, 1H1, 1H2, 1H3, 1H4, 1H5, 1H6, 1H7, 1H8, 1H9/
      Data IHC/1hC/, IHI/1hI/, IHS/1hS/, IHD/1hD/, IHT/1hT/, IHO/1hO/
C
      Symnum = One
      Linear = PG(2) .EQ. IHSt
      N      = Numer(PG)
      Test = .False.
      Do 10 I=1,10
   10     Test = Test .or. PG(3).eq.Num(I)
      If(.not.Test) N = N / 10
C
C     CI, CS, CN, CNH, CNV.
C
      If(PG(1).ne.IHC) GOTO 20
      If(Linear.or.PG(2).eq.IHI.or.PG(2).eq.IHS) goto 100
      SymNum = Float(N)
      Goto 100
C
C     DN, DNH, DND.
C
   20 If(PG(1).ne.IHD) goto 40
      SymNum = Two * Float(N)
      If(Linear) SymNum = Two
      Goto 100
C
C     SN.
C
   40 If(PG(1).ne.IHS) goto 60
      SYMNUM = Float(N) / TWO
      Goto 100
C
C     T, TD, O, OH.
C
   60 If(PG(1).eq.IHT) SymNum = Twelve
      IF (PG(1).eq.IHO) SymNum = F24
      Goto 100
C
  100 Return
      End
