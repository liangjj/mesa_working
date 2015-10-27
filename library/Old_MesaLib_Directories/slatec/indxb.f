*deck indxb
      subroutine indxb (i, ir, idx, idp)
c***begin prologue  indxb
c***subsidiary
c***purpose  subsidiary to blktri
c***library   slatec
c***type      integer (indxb-i)
c***author  (unknown)
c***see also  blktri
c***routines called  (none)
c***common blocks    cblkt
c***revision history  (yymmdd)
c   801001  date written
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c   920422  added statement so idx would always be defined.  (wrb)
c***end prologue  indxb
c
      common /cblkt/  npp        ,k          ,eps        ,cnv        ,
     1                nm         ,ncmplx     ,ik
c***first executable statement  indxb
      idx = i
      idp = 0
      if (ir) 107,101,103
  101 if (i-nm) 102,102,107
  102 idx = i
      idp = 1
      return
  103 izh = 2**ir
      id = i-izh-izh
      idx = id+id+(ir-1)*ik+ir+(ik-i)/izh+4
      ipl = izh-1
      idp = izh+izh-1
      if (i-ipl-nm) 105,105,104
  104 idp = 0
      return
  105 if (i+ipl-nm) 107,107,106
  106 idp = nm+ipl-i+1
  107 return
      end
