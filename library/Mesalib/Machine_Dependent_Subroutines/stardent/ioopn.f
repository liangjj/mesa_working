*deck @(#)ioopn.f	5.1  11/6/94
 
      subroutine ioopn(unit,file,status,iostat,extra)
c
c***begin prologue     ioopn
c***date written       870706   (yymmdd)
c***revision date      900710   (yymmdd)
c
c   10 july 1990       rlm at lanl
c     titan operating system puts scratch files on the /tmp directory.
c     this can be rediverted with the tmpdir variable in the shell --
c     be sure to export it to the children.
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
c        note that the record length is in 4-byte words on the Titan.
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
