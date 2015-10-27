*deck @(#)convgd.f	5.1  11/6/94
      subroutine convgd(a,b,result)
c***begin prologue     convgd.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             binkley, et al. (g82)
c***source             @(#)convgd.f	5.1   11/6/94
c***purpose            tests convergence
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       convgd.f
      implicit none
c     --- input variables -----
      real*8 a,b
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
      character*(*) result
c     --- scratch arrays ---
c     --- local variables ---
c
c     --- if a.lt.b, set result='yes', otherwise 'no'.
      if(a.ge.b) then
         result='no'
      else
         result='yes'
      endif
c
c
      return
      end
