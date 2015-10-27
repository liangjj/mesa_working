*deck @(#)itobyt.f	5.1  11/6/94
      function itobyt(n)
c***begin prologue     %m%
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             
c***source             @(#)itobyt.f	5.1   11/6/94
c***purpose            
c***description
c                      returns the number of bytes in an integer.
c    
c
c***references
c
c***routines called
c
c***end prologue       %m%
c
      implicit integer(a-z)
      integer itobyt
c
c
      itobyt=4*n
c
c
      return 
      end
