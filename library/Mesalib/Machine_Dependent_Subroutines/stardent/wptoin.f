*deck @(#)wptoin.f	5.1  11/6/94
      function wptoin(i)
c***begin prologue     %m%
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             
c***source             @(#)wptoin.f	5.1   11/6/94
c***purpose            returns number of integers in a working precision(real*8) 
c                        length.  
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       %m%
c
      implicit integer(a-z)
      integer wptoin
c
      wptoin=2*i
c
c
      return
      end
