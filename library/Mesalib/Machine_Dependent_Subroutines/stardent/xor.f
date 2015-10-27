*deck @(#)xor.f	5.1  11/6/94
      function xor(i,j)
c***begin prologue     %m%
c***date written       900130  
c***revision date      11/6/94      
c
c***keywords           bits, bitwise or, exlusive or, packing 
c***author             martin, richard (lanl) 
c***source             @(#)xor.f	5.1   11/6/94
c***purpose            perform a bitwise exclusive or on the two arguments. 
c***description
c                      xor(i,j) = i.xor.j
c                        this is supposedly an ansi generic name, but the titan
c                        requires a function ieor(i,j).
c***references
c
c***routines called
c
c***end prologue       %m%
c
      implicit integer(a-z)
      integer xor
c
      xor = ieor(i,j) 
c
      return
      end
