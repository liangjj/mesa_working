*deck @(#)and.f	5.1  11/6/94
      function and(i,j)
c***begin prologue     %m%
c***date written       900130  
c***revision date      11/6/94      
c
c***keywords           bits, bitwise and, packing 
c***author             martin, richard (lanl) 
c***source             @(#)and.f	5.1   11/6/94
c***purpose            perform a bitwise and on the two arguments. 
c***description
c                      and(i,j) = i.and.j
c                        this is supposedly an ansi generic name, but the titan
c                        requires a function iand(i,j).
c***references
c
c***routines called
c
c***end prologue       %m%
c
      implicit integer(a-z)
      integer and
c
      and = iand(i,j) 
c
      return
      end
