*deck inxca
      subroutine inxca (i, ir, idxa, na)
c***begin prologue  inxca
c***subsidiary
c***purpose  subsidiary to cblktr
c***library   slatec
c***type      integer (inxca-i)
c***author  (unknown)
c***see also  cblktr
c***routines called  (none)
c***common blocks    ccblk
c***revision history  (yymmdd)
c   801001  date written
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c***end prologue  inxca
      common /ccblk/  npp        ,k          ,eps        ,cnv        ,
     1                nm         ,ncmplx     ,ik
c***first executable statement  inxca
      na = 2**ir
      idxa = i-na+1
      if (i-nm) 102,102,101
  101 na = 0
  102 return
      end
