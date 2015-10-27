*deck @(#)iogetc.f	5.1  11/6/94
      subroutine iogetc(file,data,nbytes,fbyte,lbyte)
c
c***begin prologue     iogetc
c***date written       870706   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           i/o read
c***author             saxe, paul (lanl)
c***source             @(#)iogetc.f	5.1   11/6/94
c
c***purpose            to transfer data from a disc file
c
c***description        iogetc transfers 'nbytes' bytes of data from
c     the unit 'file', to 'data' starting at location 'fbyte' on that
c     unit. the first position in the file is considered to be byte 0.
c
c     n.b. this routine is machine dependent, and may restrict the byte
c     addresses and numbers of bytes to word or double word boundaries
c     and amounts.
c
c on input:
c     file        integer
c                 the unit number of the file to read to.
c
c     data        character, 'nbytes' long
c                 the data area to read to.
c
c     nbytes      integer
c                 the number of bytes to read.
c
c     fbyte       integer
c                 the byte location (starting at 0) to read the data to
c                 the file.
c
c***references         
c
c***routines called    (none)
c
c***end prologue       iogetc
c
c
      implicit integer(a-z)
c
      parameter (buflen=8192)
c
      character*4 itoc
      character*1 data(nbytes)
      character*1 buffer(buflen)
      logical print
c
      common /ioqprt/ print
      common /io/     inp,iout
c
c     ----- print if requested data on tranfer -----
c
      if (print) then
         write (iout,1) file,nbytes,fbyte
 1       format (' iogetc: file=',i3,' nbytes=',i10,' fbyte=',i10)
      end if
c
c     ----- determine the first integer word to read -----
c
      fword=1+fbyte
c
c     ----- determine the number of integer words to read -----
c
      nwords=nbytes
c
c     ---- determine which records need to be read ----
c
      first=(fword-1)/buflen + 1
      lword=fword+nwords-1
      last=(lword-1)/buflen + 1
c
c     ----- loop over the records to be read -----
c
      do 1000 record=first,last
c
c        ----- read the buffer onto record 'record' -----
c
         read (unit=file,rec=record,iostat=ierr) buffer
c
         if(ierr.ne.0) then
            call lnkerr('fortran i/o error '//itoc(ierr)//
     $           ' reading in iogetc')
         endif
c
         if (record.eq.first) then
            offset=fword-(record-1)*buflen - 1
            count=buflen-offset
            if(count.gt.nwords) count=nwords
            do 20 i=1,count
               data(i)=buffer(i+offset)
 20         continue
         else if (record.eq.last) then
            count=lword-(record-1)*buflen
            offset=nwords-count
            do 40 i=1,count
               data(i+offset)=buffer(i)
 40         continue
         else
            offset=buflen*(record-1)-fword+1
            do 50 i=1,buflen
               data(i+offset)=buffer(i)
 50         continue
         end if
 1000 continue
c
c     ----- and increment a sequential pointer -----
c
      lbyte=fbyte+nbytes
c
c
      return
      end
