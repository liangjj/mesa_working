*deck %W%  %G%
      Subroutine NoOnes(OldLen,OldStr,NewLen,NewStr)
C***Begin prologue     NoOnes
C***Date written       850601  yymmdd
C***Revision date      yymmdd  yymmdd
C***Keywords           Framework group, Stoichiometry
C***Author             Gauss82
C***Source             %W%   %G%
C***Purpose            Removes delimiters from a framework group or
C                      stoichiometry string and prepares it for printing.
C***Description
C     Call NoOnes(OldLen,OldStr,NewLen,NewStr)
C
C     Input:
C     OldStr ... An unpacked character string containing a framework
C                group or stoichiometry.
C     OldLen ... The number of characters in OldStr.
C
C     Output:
C     NewStr ... An unpacked character string containing the fwg or
C                stoich in a form suitable for printing.
C     NewLen ... The number of characters in NewStr.
C
C     Operation:
C     The first four characters of a fwg string consist of the point
C     group.  Blanks and zeroes are removed from these.  Stoich strings
C     and fwg strings contain ones following the atomic symbols.  These
C     are removed by applying the following rules:  A one is removed
C     unless the next character is a digit, a plus sign, or a minus sign
C     and the previous character is not a digit.
C
C***References
C***Routines called    (None)
C***End prologue       NoOnes
      Implicit Integer(A-Z)
      Logical FWG
      Dimension OldStr(1), NewStr(1)
      Data Zero/1h0/, Blank/1h /, Bra/1h[/, One/1h1/, Plus/1h+/,
     $     Minus/1h-/, Nine/1h9/
C
C
      FWG = OldStr(6).eq.Bra
      NewLen = 0
      Prev = Blank
      I1 = 1
      If(.not.FWG) goto 40
C
C     Deal with the point group part of an fwg string.
C
      Do 20 I = 1, 4
          If(OldStr(I).eq.Blank.or.(OldStr(I).eq.Zero.and.I.ne.3))
     $    goto 20
              NewLen = NewLen + 1
              NewStr(NewLen) = OldStr(I)
   20    Continue
      I1 = 6
C
C     Remove ones.
C
   40 Do 60 I = I1, OldLen
          Next = OldStr(I+1)
          If(I.eq.OldLen) Next = Blank
          If(I.ne.1) Prev = OldStr(I-1)
          If(OldStr(I).eq.One.and.
     $       (Next.lt.Zero.or.Next.gt.Nine).and.
     $       Next.ne.Plus.and.Next.ne.Minus.and.
     $       (Prev.lt.Zero.or.Prev.gt.Nine)) goto 60
                 NewLen = NewLen + 1
                 NewStr(NewLen) = OldStr(I)
   60    Continue
      Return
      End
