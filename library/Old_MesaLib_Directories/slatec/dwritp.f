*deck dwritp
      subroutine dwritp (ipage, list, rlist, lpage, irec)
c***begin prologue  dwritp
c***subsidiary
c***purpose  subsidiary to dsplp
c***library   slatec
c***type      double precision (swritp-s, dwritp-d)
c***author  (unknown)
c***description
c
c     write record number irecn, of length lpg, from storage
c     array list(*) onto unit number ipagef.
c     write record number irecn+1, of length lpg, onto unit
c     number ipagef from the storage array rlist(*).
c
c     to change this program unit to double precision change
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
c***end prologue  dwritp
      integer list(*)
      double precision rlist(*)
      character*8 xern1, xern2
c***first executable statement  dwritp
      ipagef=ipage
      lpg   =lpage
      irecn =irec
      write(ipagef,rec=irecn,err=100)(list(i),i=1,lpg)
      write(ipagef,rec=irecn+1,err=100)(rlist(i),i=1,lpg)
      return
c
  100 write (xern1, '(i8)') lpg
      write (xern2, '(i8)') irecn
      call xermsg ('slatec', 'dwritp', 'in dsplp, lgp = ' // xern1 //
     *   ' irecn = ' // xern2, 100, 1)
      return
      end
