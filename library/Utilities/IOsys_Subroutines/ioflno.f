*deck @(#)ioflno.f	5.1  11/6/94
      function ioflno(key,nfiles,pt,keylis,error)
c
c***begin prologue     ioflno
c***date written       850125   (yymmdd)
c***revision date      860112   (yymmdd)
c***keywords           iosys dependent routines
c
c***author             saxe, paul,    (lanl)
c***source             @(#)ioflno.f	5.1   11/6/94
c***purpose            to find the position of a file in the file list
c
c***description        #
c
c
c***references
c
c***routines called    (none)
c
c   common blocks:     (none)
c
c***end prologue       ioflno
c
      implicit integer (a-z)
      integer ioflno
c
      character*(*) key,keylis(*)
c
      error=0
      do 1 ioflno=pt+1,pt+nfiles
         if (key.eq.keylis(ioflno)) return
    1 continue
c
      error=1
      return
      end
