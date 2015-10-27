*deck @(#)prtnam.f	1.1  11/30/90
      subroutine prtnam(messge,filnam,iout)
c***begin prologue     prtnam
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           print, file name, input/output
c***author             martin, richard (lanl)
c***source
c***purpose            prints a message and the local file name associated
c                      with the one known to mesa.
c***description
c                      call prtnam(messge,filnam,iout)
c                        messge   a message to be printed (character*120).
c                        filnam   the name of the file known to mesa
c                                 'inp','out','chk',etc. (character*8).
c                        iout     unit number of the output file.
c
c***references
c***routines called    streqc(chr)
c***end prologue       prtnam
      implicit integer(a-z)
      character*(*) messge,filnam
      character*128 name
      character strout*120
      logical streqc
c
 1000 format(1x,131a1)
c
c     module to associate an external file name with the internal one,
c     and to print out a message.
c
      if(streqc(filnam,'rwf')) then
         call iosys('read character"read-write filename" from rwf',
     $        0,0,0,name)
      else if(streqc(filnam,'chk')) then
         call iosys('read character "checkpoint filename" from rwf',
     $        0,0,0,name)
      endif
c
      strout=messge//name
      write(iout,1000) (strout(i:i),i=1,len(strout))
c
c
      return
      end
