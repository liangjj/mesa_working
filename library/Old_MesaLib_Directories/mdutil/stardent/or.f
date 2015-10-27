*deck @(#)or.f	5.1  11/6/94
      function or(i,j)
c***begin prologue     %m%
c***date written       900130  
c***revision date      11/6/94      
c
c***keywords           bits, bitwise or, packing 
c***author             martin, richard (lanl) 
c***source             @(#)or.f	5.1   11/6/94
c***purpose            perform a bitwise or on the two arguments. 
c***description
c                      or(i,j) = i.or.j
c                        this is supposedly an ansi generic name, but the titan
c                        requires a function ior(i,j).
c***references
c
c***routines called
c
c***end prologue       %m%
c
      implicit integer(a-z)
      integer or
c
      or = ior(i,j) 
c
      return
      end
