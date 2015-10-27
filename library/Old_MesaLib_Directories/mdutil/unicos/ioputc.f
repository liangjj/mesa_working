*deck @(#)ioputc.f	5.1  11/6/94
      subroutine ioputc(file,data,nbytes,fbyte,lbyte)
c
c***begin prologue     ioputc
c***date written       870706   (yymmdd)
c***revision date      910617   (yymmdd)
c
c   17 june   1991     rlm at lanl
c      unicos version
c***keywords           i/o write
c***author             saxe, paul (lanl)
c***source             @(#)ioputc.f	5.1   11/6/94
c
c***purpose            to transfer character data to a disc file
c
c***description        ioputc transfers 'nbytes' bytes of data from 
c     'data' to the unit 'file', starting at location 'fbyte' on that 
c     unit. the first position in the file is considered to be byte 0.
c
c     unicos... note that for the purposes of putwa the first word on
c     disk is considered to be 1. also note that putwa expects the name
c     of the file, not its unit number as passed in variable file.
c     n.b. this routine is machine dependent, and may restrict the byte
c     addresses and numbers of bytes to word or double word boundaries
c     and amounts.
c
c on input:
c     file        integer
c                 the unit number of the file to write to.
c
c     data        character 'nbytes' long
c                 the data to write.
c
c     nbytes      integer
c                 the number of bytes to write.
c
c     fbyte       integer
c                 the byte location (starting at 0) to write the data
c                 to the file.
c
c***references         
c
c***routines called    (none)
c
c***end prologue       ioputc
c
c
      implicit integer(a-z)
c
      parameter (itobyt=8,maxun=20)
c
      character*4 itoc
      character*256 unitnm
      character*8 filenm
      character*1 data(*)
      character*1 buffer
      logical print
c
      common /ioqprt/ print
      common /io/     inp,iout
      common /fnames/ unitnm(maxun)
c
c     ----- print if requested data on tranfer -----
c
      if (print) then
         write (iout,1) file,nbytes,fbyte
 1       format (' ioputc: file=',i3,' nbytes=',i10,' fbyte=',i10)
      end if
c
c     ----- determine the first character word to write -----
c
      fword=1+fbyte/itobyt
c
c     ----- check for non-integral starting position -----
c
      if ((fword-1)*itobyt.ne.fbyte) then
         call lnkerr('ioput: non-integral starting position for write')
      end if
c
c     ----- determine the number of integer words to write -----
c
      nwords=nbytes/itobyt
c
c     ----- check for non-integral number of words -----
c
      if (nwords*itobyt.ne.nbytes) then
         call lnkerr('ioput: non-integral number of words to write')
      end if
c
c     ----- write -----
c
      filenm=unitnm(file)
      call putwa(filenm,data,fword,nwords,ierr)
      if(ierr.ne.0) then
         call lnkerr(' putwa i/o error:'//itoc(ierr))
      endif
c
c     ----- and increment a sequential pointer -----
c
      lbyte=fbyte+nbytes
c
c
      return
      end
