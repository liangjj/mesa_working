*deck @(#)ioopn.f	5.1   11/6/94
 
      subroutine ioopn(unit,file,status,iostat,extra)
c
c***begin prologue     ioopn
c***date written       870706   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           iosys dependent routine
c***author             saxe, paul (lanl)
c***source             @(#)ioopn.f	5.1   11/6/94
c
c***purpose            to open a unit to a disc file 
c
c***description        
c
c***references         
c
c***routines called    (none)
c
c***end prologue       ioopn
c
      implicit integer (a-z)
c
      parameter (buflen=2048)
c
      character*(*) file
      character*(*) status
      character*(*) extra
      integer unit
      integer iostat
c
c     ----- open the unit as a fortran unit -----
c
      if (status.eq.'scratch') then
c
c        note that the record length is in 4-byte  words in ultrix.
         open (unit=unit,
     $        recl=buflen,
     $        access='direct',
     $        status=status,
     $        iostat=iostat)
      else
         open (unit=unit,
     $        file=file,
     $        recl=buflen,
     $        access='direct',
     $        status=status,
     $        iostat=iostat)
      end if
c
c
      return
      end
