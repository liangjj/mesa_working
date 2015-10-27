*deck inxcc
      subroutine inxcc (i, ir, idxc, nc)
c***begin prologue  inxcc
c***subsidiary
c***purpose  subsidiary to cblktr
c***library   slatec
c***type      integer (inxcc-i)
c***author  (unknown)
c***see also  cblktr
c***routines called  (none)
c***common blocks    ccblk
c***revision history  (yymmdd)
c   801001  date written
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c***end prologue  inxcc
      common /ccblk/  npp        ,k          ,eps        ,cnv        ,
     1                nm         ,ncmplx     ,ik
c***first executable statement  inxcc
      nc = 2**ir
      idxc = i
      if (idxc+nc-1-nm) 102,102,101
  101 nc = 0
  102 return
      end
