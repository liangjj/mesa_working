*deck %W%  %G%
      subroutine ioputc(file,data,nbytes,fbyte,lbyte)
c
c***begin prologue     ioputc
c***date written       870706   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           i/o write
c***author             saxe, paul (lanl)
c***source             %W%   %G%
c
c***purpose            to transfer character data to a disc file
c
c***description        ioputc transfers 'nbytes' bytes of data from
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
c..bhl
      common /bhl/ maxwrd(100)
c..bhl
c
      parameter (itobyt=8)
c
      character*4 itoc
      character*1 data(*)
      character*1 buffer
      logical print
c
      common /ioqprt/ print
      common /io/     inp,iout
c
c     ----- print if requested data on tranfer -----
c
      if (print) then
         write (iout,1) file,nbytes,fbyte
 1       format (' ioputc: file=',i3,' nbytes=',i10,' fbyte=',i10)
      end if
c
c     ----- determine the first integer word to write -----
c
      fword=fbyte/itobyt
c
c     ----- check for non-integral starting position -----
c
      if (fword*itobyt.ne.fbyte) then
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
      if(maxwrd(file).gt.fword+nwords+1) then
      call wrabs(file,data,nwords,fword)
      else
c
c     call lnkerr(' cant extend a file size ')
c
       need=fword+nwords+1-maxwrd(file)
       need=( (need-1)/262144 + 1 ) * 262144
       call xtendabs(file,need,newsiz,ierr)
       if(ierr.ne.0) call lnkerr(' cant extend a file size ')
       call wrabs(file,data,nwords,fword)
       maxwrd(file)=newsiz
c
 
      end if
c
c     ----- and wait for i/o to complete -----
c
c     ----- and wait for i/o to complete -----
c
      if (unit(file)) 12,11,10
 10   continue
         call lnkerr('i/o error')
 11   continue
         call lnkerr('eoi or eof error in i/o')
 12   continue
c
c     ----- and increment a sequential pointer -----
c
      lbyte=fbyte+nbytes
c
c
      return
      end
