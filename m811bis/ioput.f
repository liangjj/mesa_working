*deck @(#)ioput.f	5.1  11/6/94
      subroutine ioput(file,data,nbytes,fbyte,lbyte)
c
c***begin prologue     ioput
c***date written       870706   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           i/o write
c***author             saxe, paul (lanl)
c***source             @(#)ioput.f	5.1   11/6/94
c
c***purpose            to transfer data to a disc file
c
c***description        ioput transfers 'nbytes' bytes of data from 
c     'data' to the unit 'file', starting at location 'fbyte' on that 
c     unit. the first position in the file is considered to be byte 0.
c
c     n.b. this routine is machine dependent, and may restrict the byte
c     addresses and numbers of bytes to word or double word boundaries
c     and amounts.
c
c on input:
c     file        integer
c                 the unit number of the file to write to.
c
c     data        integer, real or logical. 'nbytes' long
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
c***end prologue       ioput
c
c
      implicit integer(a-z)
c
      parameter (buflen=2048)
c      parameter (buflen=4096)
      parameter (itobyt=4)
c
      character*4 itoc
      integer data(*)
      integer buffer
      logical print
c
      common /ioqprt/ print
      common /ioqbuf/ buffer(buflen)
      common /io/     inp,iout
c
c     ----- print if requested data on tranfer -----
c
      if (print) then
         write (iout,1) file,nbytes,fbyte
 1       format (' ioput:  file=',i3,' nbytes=',i10,' fbyte=',i10)
      end if
c
c     ----- determine the first integer word to write -----
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
c     ---- determine which records need to be written ----
c
      first=(fword-1)/buflen + 1
      lword=fword+nwords-1
      last=(lword-1)/buflen + 1
c
c     ----- loop over the records to be written -----
c
      do 1000 record=first,last
c
c        ----- update the buffer -----
c
         if (record.eq.first) then
c
c           ----- read in the first record...noting that it may
c                 not exist
c
            read (unit=file,rec=record,err=10) buffer
 10         continue
c
            offset=fword-(record-1)*buflen - 1
            count=buflen-offset
            if(count.gt.nwords) count=nwords
            do 20 i=1,count
               buffer(i+offset)=data(i)
 20         continue
         else if (record.eq.last) then
c
c           ----- read in the last record...noting that it may
c                 not exist
c
            read (unit=file,rec=record,err=30) buffer
 30         continue
c
            count=lword-(record-1)*buflen
            offset=nwords-count
            do 40 i=1,count
               buffer(i)=data(i+offset)
 40         continue
         else
            offset=buflen*(record-1)-fword+1
            do 50 i=1,buflen
               buffer(i)=data(i+offset)
 50         continue
         end if
c
c        ----- write the buffer onto record 'record' -----
c
         write (unit=file,rec=record,iostat=ierr) buffer
c         if(ierr.ne.0) then
c            call lnkerr('fortran i/o error '//itoc(ierr)//
c     $           ' writing in ioput')
c         endif
c
 1000 continue
c
c     ----- and increment a sequential pointer -----
c
      lbyte=fbyte+nbytes
c
c
      return
      end




