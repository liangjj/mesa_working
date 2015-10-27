*deck pimach
      function pimach (dum)
c***begin prologue  pimach
c***subsidiary
c***purpose  subsidiary to hstcsp, hstssp and hwscsp
c***library   slatec
c***type      single precision (pimach-s)
c***author  (unknown)
c***description
c
c     this subprogram supplies the value of the constant pi correct to
c     machine precision where
c
c     pi=3.1415926535897932384626433832795028841971693993751058209749446
c
c***see also  hstcsp, hstssp, hwscsp
c***routines called  (none)
c***revision history  (yymmdd)
c   801001  date written
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c***end prologue  pimach
c
c***first executable statement  pimach
      pimach = 3.14159265358979
      return
      end
