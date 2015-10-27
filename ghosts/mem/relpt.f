      function relpt(name)
c***begin prologue     relpt
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
c***end prologue       relpt 
      implicit integer(a-z)
      integer relpt
      character*(*) name
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
c
      if (debug) then
         write(iout,*) 'relpt:releasing space for ',name
         write(iout,*) 'current array table'
         call memdat
      end if
c
c     find the array
      do 10 i=1,narray
         if(handle(i).eq.name) then
            array=i
            go to 20
         endif
   10 continue
      call lnkerr('relpt:cant find array '//name)
   20 continue
c
c
c     release the space.
      handle(array)='free'
c
c     check on either side to see if two free spaces are now continguous.
      total=narray
      do 100 i=0,total
         if (addr(i)+length(i).eq.addr(array)) then
c           found the one below
            if(handle(i).eq.'free') then
               length(i)=length(i)+length(array)
               if(debug) then
                  write(iout,*) 'merging free space with array ',i
               endif
c              remove the array and pop the stack
               do 30 j=array,narray-1
                  handle(j)=handle(j+1)
                  length(j)=length(j+1)
                  addr(j)=addr(j+1)
   30          continue
               handle(narray)=' '
               addr(narray)=0
               length(narray)=0
c              reset the current array index
               array=i
               narray=narray-1
            endif
         else if (addr(array)+length(array).eq.addr(i)) then
c           found the one above
            if(handle(i).eq.'free') then
               if(debug) then
                  write(iout,*) 'merging free space with array ',i
               endif
               addr(i)=addr(array)
               length(i)=length(i)+length(array)
c              remove the array and pop the stack
               do 40 j=array,narray-1
                  handle(j)=handle(j+1)
                  length(j)=length(j+1)
                  addr(j)=addr(j+1)
   40          continue
               handle(narray)=' '
               addr(narray)=0
               length(narray)=0
c              reset the current array index
               array=i
               narray=narray-1
            endif
         endif
  100 continue
      if(debug) then
         write(iout,*) 'revised array table'
         call memdat
      endif
c
c     return pointer of 0 if successful
      relpt=0
c
c
      return
      end
