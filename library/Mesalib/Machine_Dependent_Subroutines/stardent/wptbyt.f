*deck @(#)wptbyt.f	5.1  11/6/94
      function wptbyt(n)
c***begin prologue     %m%
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             
c***source             @(#)wptbyt.f	5.1   11/6/94
c***purpose            returns number of bytes in a working precision(real*8) 
c                        number.  
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
      integer wptbyt
c
      wptbyt=8*n
c
c
      return
      end
