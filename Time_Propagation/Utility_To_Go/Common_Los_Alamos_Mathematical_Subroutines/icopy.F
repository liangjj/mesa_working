*deck @(#)icopy.f	5.1  11/6/94
      subroutine icopy(ivec,jvec,n)
c***begin prologue     icopy.f
c***date written       yymmdd  
c***revision date      7/16/93      
c
c***keywords           
c***author             unknown
c***source             @(#)icopy.f	5.1   11/6/94
c***purpose            copies integer vector ivec(n) -> jvec(n) 
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       icopy.f
      implicit none
c     --- input variables -----
      integer n
c     --- input arrays (unmodified) ---
      integer ivec(n)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      integer jvec(n)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer i
c
c
      do 100 i=1,n
         jvec(i)=ivec(i)
  100 continue
c
c
      return
      end
