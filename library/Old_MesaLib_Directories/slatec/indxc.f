*deck indxc
      subroutine indxc (i, ir, idxc, nc)
c***begin prologue  indxc
c***subsidiary
c***purpose  subsidiary to blktri
c***library   slatec
c***type      integer (indxc-i)
c***author  (unknown)
c***see also  blktri
c***routines called  (none)
c***common blocks    cblkt
c***revision history  (yymmdd)
c   801001  date written
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c***end prologue  indxc
      common /cblkt/  npp        ,k          ,eps        ,cnv        ,
     1                nm         ,ncmplx     ,ik
c***first executable statement  indxc
      nc = 2**ir
      idxc = i
      if (idxc+nc-1-nm) 102,102,101
  101 nc = 0
  102 return
      end
