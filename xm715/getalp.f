*deck %W%  %G%
      subroutine getalp(ai,aj,ak,al,index,alpha,nv,len,camcb,ca,cb,
     $                  nij,nkl)
c***begin prologue     getalp.f
c***date written       851113  
c***revision date      11/6/94      
c
c***keywords           
c***author             saxe, paul(lanl) 
c***source             %W%   %G%
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       getalp.f
      implicit none
c     --- input variables -----
      integer nv,len,nij,nkl
c     --- input arrays (unmodified) ---
      integer index(len,6)
      real*8 ai(*),aj(*),ak(*),al(*)
      real*8 ca(nij,3),cb(nkl,3)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 alpha(nv,4),camcb(nv,3)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer i,coord
c
c
      do 1 i=1,nv
         alpha(i,1)=ai(index(i,1))
    1 continue
      do 2 i=1,nv
         alpha(i,2)=aj(index(i,2))
    2 continue
      do 3 i=1,nv
         alpha(i,3)=ak(index(i,3))
    3 continue
      do 4 i=1,nv
         alpha(i,4)=al(index(i,4))
    4 continue
c
      do 6 coord=1,3
         do 5 i=1,nv
            camcb(i,coord)=ca(index(i,5),coord)-cb(index(i,6),coord)
    5    continue
    6 continue
c
c
      return
      end
