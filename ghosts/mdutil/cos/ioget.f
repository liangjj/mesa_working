*deck %W%  %G%
      subroutine ioget(file,data,nbytes,fbyte,lbyte)
c
c***begin prologue     ioget
c***date written       870706   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           i/o read
c***author             saxe, paul (lanl)
c***source             %W%   %G%
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
      parameter (itobyt=8)
c
      character*4 itoc
      integer data(*)
      integer buffer
      logical print
c
      common /ioqprt/ print
      common /io/     inp,iout
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
      fword=fbyte/itobyt
c
c     ----- check for non-integral starting position -----
c
      if (fword*itobyt.ne.fbyte) then
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
c     ----- read -----
c
      call rdabsf(file,data,nwords,fword)
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
