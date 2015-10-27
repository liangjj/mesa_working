*deck @(#)shiftr.f	5.1  11/6/94
      function shiftr(i,j)
c***begin prologue     %m%
c***date written       900130  
c***revision date      11/6/94      
c
c***keywords           shift, right shift, bits, packing 
c***author             martin, richard (lanl) 
c***source             @(#)shiftr.f	5.1   11/6/94
c***purpose            right bit shift 
c***description
c     
c                      shiftr(i,j) shifts the bits in the first argument to the
c                         right by j places.  in this implementation for the
c                         titan, bits shifted out from the right are lost, and
c                         zeroes are shifted in from the left side. 
c                         note the intrinsic titan function assumes j <= 0 . 
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
      integer shiftr
c
c
      shiftr=ishft(i,-j)
c
c
      return
      end
