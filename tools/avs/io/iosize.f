*deck @(#)iosize.f	4.2  7/16/93
      function iosize(type)
c***begin prologue     iosize.f
c***date written       850125  
c***revision date      7/16/93      
c
c***keywords           io, i/o, files, input, output
c***author             saxe, paul(lanl) 
c***source             @(#)iosize.f	4.2   7/16/93
c***purpose            returns the byte size of word types
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       iosize.f
c
      implicit integer (a-z)
      integer iosize
c
      character*(*) type
c
      if (type(1:1).eq.'r') then
         iosize=wptbyt(1)
      else if (type(1:1).eq.'i') then
         iosize=itobyt(1)
      else if (type(1:1).eq.'c') then
         iosize=1
      else
         call lnkerr('iosys has problems with the size/type of an '//
     #               'entry in the directory...this is serious')
      end if
c
c
      return
      end
