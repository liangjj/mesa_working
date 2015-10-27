*deck dreadp
      subroutine dreadp (ipage, list, rlist, lpage, irec)
c***begin prologue  dreadp
c***subsidiary
c***purpose  subsidiary to dsplp
c***library   slatec
c***type      double precision (sreadp-s, dreadp-d)
c***author  (unknown)
c***description
c
c     read record number irecn, of length lpg, from unit
c     number ipagef into the storage array, list(*).
c     read record  irecn+1, of length lpg, from unit number
c     ipagef into the storage array rlist(*).
c
c     to convert this program unit to double precision change
c     /real (12 blanks)/ to /double precision/.
c
c***see also  dsplp
c***routines called  xermsg
c***revision history  (yymmdd)
c   811215  date written
c   890605  corrected references to xerrwv.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   900510  convert xerrwv calls to xermsg calls.  (rwc)
c***end prologue  dreadp
      integer list(*)
      double precision rlist(*)
      character*8 xern1, xern2
c***first executable statement  dreadp
      ipagef=ipage
      lpg   =lpage
      irecn=irec
      read(ipagef,rec=irecn,err=100)(list(i),i=1,lpg)
      read(ipagef,rec=irecn+1,err=100)(rlist(i),i=1,lpg)
      return
c
  100 write (xern1, '(i8)') lpg
      write (xern2, '(i8)') irecn
      call xermsg ('slatec', 'dreadp', 'in dsplp, lpg = ' // xern1 //
     *   ' irecn = ' // xern2, 100, 1)
      return
      end
