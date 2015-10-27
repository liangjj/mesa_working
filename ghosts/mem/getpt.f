      function getpt(type,nwords,name,heap)
c***begin prologue     getpt
c***date written       910801    yymmdd  
c***revision date      yymmdd
c
c***keywords           memory,memory manager,core 
c***author             martin, richard (lanl) 
c***source             %W%   %G%
c***purpose            
c***description
c     
c    
c
c***references         (none)
c
c***routines called
c       start here
c
c***end prologue       getpt
      implicit integer(a-z)
      integer getpt
      character*(*) type
      integer nwords
      character*(*) name
      logical heap
      logical debug
c
      parameter (maxarr=1000) 
c
      character*8 handle
      integer narray,addr,length
c
      data debug/.true./
c
      common /mem1/ narray,addr(0:maxarr),length(0:maxarr)
      common /mem2/ handle(0:maxarr)
      common /io/ inp,iout
c
c
c     how many bytes are requested?
      if(debug) then
         write(iout,*) ' getpt:type nwords name heap? '
         write(iout,*) '     ',type,' ',nwords,' ',name,' ',heap
         write(iout,*) ' current array table'
         call memdat
      endif
      if(type.eq.'real') then
         request=itobyt(wptoin(nwords))
      else if(type.eq.'integer') then
ctemp
c        always request an even number of words for integers
c        this ensures that all array addresses will begin on
c        8-byte boundaries.  this is helpful on some machines
c        which prefer real*8 addresses to satisfy this criterion. 
         if (mod(nwords,2).eq.0) then
            request=itobyt(nwords)
         else
            request=itobyt(nwords+1)
        endif
      else
         call lnkerr(' unrecognized data type in getpt.')
      endif
c
c
      if (heap) then
c        the user has requested the space be allocated from the expansion
c        heap. allocate the space from array '0', the free space area
c        which extends from the top of the heap to the end of
c        the original allocation.
         if(debug) then
            write(iout,*) ' request to allocate from heap' 
            write(iout,*) ' arrayno:',narray+1
         endif
         array=narray+1
         handle(array)=name
         addr(array)=addr(0)
         length(array)=request
         narray=narray+1
c
c        reduce the amount of space left in the last free block,
c        and decrement the amount of available space in array '0'
         handle(0)='free'
         addr(0)=addr(0)+request
         length(0)=length(0)-request
c
c        make sure it fits
         if(length(0).lt.0) then
            write(iout,*) ' memory: request exceeds available space'
            call lnkerr('request exceeds available space')
         endif
      else
c
c        try to find an open space.
         do 10 i=1,narray
            if(handle(i).eq.'free') then
               if(length(i).eq.request) then
c                 it fits exactly into this spot.
                  if(debug) then
                     write(iout,*) ' array fits into old spot ',
     $                             i,' exactly '
                  endif
                  array=i
                  handle(i)=name
                  length(i)=request 
                  addr(i)=addr(i)
                  go to 20
               else if(length(i).gt.request) then
c                 if fits in a region previously released.
                  if (debug) then
                     write(iout,*) ' array fits into old spot ',i
                  endif
                  array=i
                  handle(i)=name
                  left=length(i)-request
                  length(i)=request
                  addr(i)=addr(i)
c                 assign name to remaining space.  
                  narray=narray+1
                  handle(narray)='free'
                  length(narray)=left
                  addr(narray)=addr(i)+request
                  go to 20
               endif
            endif
   10    continue
c
c        it won't fit into an existing hole, so try to put it in
c        array '0'.
         if(request.le.length(0)) then
            if(debug) then
               write(iout,*) ' allocating from free expansion space'
            endif
            narray=narray+1
            array=narray
            handle(array)=name
            addr(array)=addr(0)
            length(array)=request
c           adjust array 0
            addr(0)=addr(0)+request
            length(0)=length(0)-request
         else
c           it won't fit anywhere.
            call lnkerr('wont fit')
         endif
      endif
c
c
   20 continue
      if(debug) then
         write(iout,*) ' new array table'
         call memdat
      endif
c
c     return the pointer
      getpt=addr(array)
c

      return
      end
