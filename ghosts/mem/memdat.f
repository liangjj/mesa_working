      subroutine memdat
c***begin prologue     memdat
c***date written       910801    yymmdd  
c***revision date      yymmdd
c
c***keywords           
c***author             martin, richard (lanl) 
c***source             %W%   %G%
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       memdat 
      implicit integer(a-z)
c
      parameter (maxarr=1000)
c
      character*8 handle
      integer narray
      integer addr
      integer length
      logical debug
c
      data debug/.true./
c
      common /mem1/ narray,addr(0:maxarr),length(0:maxarr)
      common /mem2/ handle(0:maxarr)
      common /io/ inp,iout
c
 1000 format(1x,'number',5x,'handle',4x,'address',5x,'length(bytes)')
 1010 format(1x,i6,3x,a8,2x,i8,5x,i13)
c
      write(iout,1000)
      do 100 i=0,narray
         write(iout,1010) i,handle(i),addr(i),length(i)
  100 continue
c
c
      return
      end
