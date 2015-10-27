      function inimem(a,maxsiz)
c***begin prologue     inimem
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
c***end prologue       inimem 
      implicit integer(a-z)
      integer inimem
      integer maxsiz
      integer a(maxsiz)
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
c     initialize array table
      narray=0
      do 100 i=1,maxarr
         addr(i)=0
         length(i)=0
         handle(i)=' '
  100 continue
c
c     initialize expansion array
      addr(0)=loc(a(1))
      length(0)=itobyt(maxsiz)
      handle(0)='free'
c
c
      inimem=0
c
c
      return
      end
