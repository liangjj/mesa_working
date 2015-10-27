*deck @(#)ifill.f	1.1  11/30/90
      subroutine ifill(iv,is,n)
c***begin prologue     ifill
c***date written       
c***revision date      
c***keywords           fill, load
c***author             
c***source             
c***purpose            load:  iv=is.
c***description
c                      call ifill(iv,is,n)
c                        iv       output vector of length n.
c                        is       scalar to load.
c                        n        length of vector.
c
c***references
c***routines called    (none)
c***end prologue       ifill
      dimension iv(n)
      do 1 i=1,n
         iv(i)=is
    1 continue
      return
      end
