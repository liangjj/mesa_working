*deck @(#)pm992.f	1.3  7/30/91
      subroutine pm992
c***begin prologue     m992
c***date written       020731   (yymmdd)
c***revision date               (yymmdd)
c
c***keywords           m992, link 992, io test
c***author             Schneider, Barry (NSF)
c***source             @(#)m992.f
c***purpose            test i/o under new compaq compiler/linux
c***description
c
c***references
c
c***routines called
c***end prologue       pm992
c
      implicit integer (a-z)
c
      parameter (maxnbf=2000)
c
      character*4096 ops
      character*8 rtine
      integer need, ngot
      integer a
      real*8 z
      pointer(p,a(1)), (p,z(1))
c
      common /io/ inp, iout
c
c
 1000 format(1x,'m992:io test',//)
 1001 format(5x,'memory use        ',11x,i9)
      write(iout,1000)      
      read(inp,*) data, reclen
      p=malloc(8*reclen)
      write(iout,*) '    no. 8 byte data elements to write  = ',data
      write(iout,*) '    record length in 4 byte words      = ',reclen
      call iostd(z,a,data,reclen)
c
c     ----- and exit with grace -----
c
c
c
      stop
      end
