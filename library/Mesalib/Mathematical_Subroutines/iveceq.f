*deck @(#)iveceq.f	5.1  11/6/94
      function iveceq(a,b,n)
c***begin prologue     iveceq.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             unknown
c***source             @(#)iveceq.f	5.1   11/6/94
c***purpose            compares two integer vectors
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       iveceq.f
      implicit none
c     --- input variables -----
      integer n
c     --- input arrays (unmodified) ---
      integer a(n),b(n)
c     --- input arrays (scratch) ---
c     --- function returned ---
      logical iveceq
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer i
c
      iveceq=.true.
      do 1 i=1,n
         if (a(i).ne.b(i)) then
            iveceq=.false.
            return
         endif
    1 continue
c
c
      return
      end
