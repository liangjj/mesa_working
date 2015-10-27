*deck sclosm
      subroutine sclosm (ipage)
c***begin prologue  sclosm
c***subsidiary
c***purpose  subsidiary to splp
c***library   slatec
c***type      all (sclosm-a)
c***author  (unknown)
c***description
c
c     1. unload, release, or close unit number ipagef.
c
c***see also  splp
c***routines called  xermsg
c***revision history  (yymmdd)
c   811215  date written
c   890605  corrected references to xerrwv.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c   900510  convert xerrwv calls to xermsg calls.  (rwc)
c***end prologue  sclosm
      character*8 xern1
c
c***first executable statement  sclosm
      ipagef=ipage
      close(unit=ipagef,iostat=ios,err=100,status='keep')
      return
c
  100 write (xern1, '(i8)') ios
      call xermsg ('slatec', 'sclosm',
     *   'in splp, close has error flag = ' // xern1, 100, 1)
      return
      end
