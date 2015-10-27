*deck %W%  %G%
      Subroutine FillOp(Ops,NewOp)
C***Begin prologue     FillOp
C***Date written       850601  yymmdd
C***Revision date      yymmdd  yymmdd
C***Keywords           Options
C***Author             Martin, Richard (LANL)
C***Source             %W%   %G%
C***Purpose            Adds an option to the option string.
C***Description
C     Call FillOp(Ops,NewOp)
C       Ops     The options string.
C       NewOp   The option to add.
C
C     This routine will add NewOp to the Ops string if it doesn't
C     already exist there.
C***References
C***Routines called    CSkipB(Chr), PakStr(Chr)
C***End prologue       FillOp
      Implicit Integer(A-Z)
      Character*(*) Ops,NewOp
      Character*1 Comma,Blank
C
      Data Comma/','/, Blank/' '/
C
C     First find the end of the options string.
      LenOps=CSkipB(Ops,Blank)
C
C     Pack down the new option.
      Call PakStr(NewOp,LenNew)
      If(LenNew.eq.0) Return
C
C     See if the new option already exists in the string.
      If(Index(Ops,NewOp(:LenNew)).ne.0) Return
C
C     Place it in the options string followed by a comma.
      Ops(LenOps+1:)=NewOp(:LenNew)//Comma
C
C
      Return
      End
