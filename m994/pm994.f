*deck @(#)pm994.f	1.3  7/30/91
      subroutine pm994
c***begin prologue     m994
c***date written       020731   (yymmdd)
c***revision date               (yymmdd)
c
c***keywords           m994, link 994, io test
c***author             Schneider, Barry (NSF)
c***source             @(#)m994.f
c***purpose            test i/o under new compaq compiler/linux
c***description
c
c***references
c
c***routines called
c***end prologue       pm994
c
      implicit integer (a-z)
c
      integer need, ngot
      integer a
      real*8 z
      pointer(p,a(1)), (p,z(1))
c
      common /io/     inp,iout
c
c
 1000 format(1x,'m994:io test',//)
 1001 format(5x,'memory use        ',11x,i9)
c
c     ----- recover the options string -----
c
      write(iout,1000)      
      read(inp,*) data,lenbin,reclen
      p=malloc(4*lenbin)
      write(iout,*) '    length of bin               = ',lenbin
      write(iout,*) '    no. data elements to write  = ',data
      call iostd(z(1),a(1),data,lenbin,reclen)
c      call readran(z(1),a(1),data,lenbin,reclen)
      stop
      end
