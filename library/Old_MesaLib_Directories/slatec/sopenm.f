*deck sopenm
      subroutine sopenm (ipage, lpage)
c***begin prologue  sopenm
c***subsidiary
c***purpose  subsidiary to splp
c***library   slatec
c***type      all (sopenm-a)
c***author  (unknown)
c***description
c
c     1. open unit number ipagef as a random access file.
c
c     2. the record length is constant=lpg.
c
c***see also  splp
c***routines called  xermsg
c***revision history  (yymmdd)
c   811215  date written
c   890605  corrected references to xerrwv.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c   900510  convert xerrwv calls to xermsg calls.  (rwc)
c***end prologue  sopenm
      character*8 xern1
c
c***first executable statement  sopenm
      ipagef=ipage
      lpg   =lpage
      open(unit=ipagef,iostat=ios,err=100,status='unknown',
     *access='direct',form='unformatted',recl=lpg)
      return
c
 100  write (xern1, '(i8)') ios
      call xermsg ('slatec', 'sopenm',
     *   'in splp, open has error flag = ' // xern1, 100, 1)
      return
      end
