*deck @(#)ioget.f	5.1  11/6/94
      subroutine ioget(file,data,nbytes,fbyte,lbyte)
c
c***begin prologue     ioget
c***date written       870706   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           i/o read
c***author             saxe, paul (lanl)
c***source             @(#)ioget.f	5.1   11/6/94
c
c***purpose            to transfer data from a disc file
c
c***description        ioget transfers 'nbytes' bytes of data from
c     the unit 'file', to 'data' starting at location 'fbyte' on that
c     unit. the first position in the file is considered to be byte 0.
c
c     n.b. this routine is machine dependent, and may restrict the byte
c     addresses and numbers of bytes to word or double word boundaries and
c     amounts.
c
c on input:
c     file        integer
c                 the unit number of the file to read to.
c
c     data        integer, real or logical. 'nbytes' long
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
c***end prologue       ioget
c
c
      implicit integer(a-z)
c
      parameter (buflen=2048)
c      parameter (buflen=8192)
      parameter (itobyt=4)
c
      character*4 itoc
      integer data(*)
      integer buffer
      logical print
c      logical writit
c
      common /ioqprt/ print
      common /ioqbuf/ buffer(buflen)
      common /io/     inp,iout
c      writit=.true.
c
c     ----- print if requested data on tranfer -----
c
      if (print) then
         write (iout,1) file,nbytes,fbyte
 1       format (' ioget:  file=',i3,' nbytes=',i10,' fbyte=',i10)
      end if
c
c     ----- determine the first integer word to read -----
c
      fword=1+fbyte/itobyt
c
c     ----- check for non-integral starting position -----
c
      if ((fword-1)*itobyt.ne.fbyte) then
         call lnkerr('ioget: non-integral starting position for read')
      end if
c
c     ----- determine the number of integer words to read -----
c
      nwords=nbytes/itobyt
c
c     ----- check for non-integral number of words -----
c
      if (nwords*itobyt.ne.nbytes) then
         call lnkerr('ioget: non-integral number of words to read')
      end if
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
c         if(writit) then
c            write(iout,*) file, record, ierr
c         endif
c
c         if(ierr.ne.0) then
c            call lnkerr('fortran i/o error '//itoc(ierr)//
c     $           ' reading in ioget')
c         endif
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
