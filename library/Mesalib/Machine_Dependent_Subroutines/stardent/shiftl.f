*deck @(#)shiftl.f	5.1  11/6/94
      function shiftl(i,j)
c***begin prologue     %m%
c***date written       900130  
c***revision date      11/6/94      
c
c***keywords           shift, left shift, bits, packing 
c***author             martin, richard (lanl) 
c***source             @(#)shiftl.f	5.1   11/6/94
c***purpose            left bit shift 
c***description
c     
c                      shiftl(i,j) shifts the bits in the first argument to the
c                         left by j places.  in this implementation for the
c                         titan, bits shifted out from the left are lost, and
c                         zeroes are shifted in from the right side. 
c                         note the intrinsic titan function assumes j >= 0 . 
c
c***references
c
c***routines called
c     %m%
c       start here
c
c***end prologue       %m%
c
      implicit integer(a-z)
      integer shiftl
c
c
      shiftl=ishft(i,j)
c
c
      return
      end
