*deck @(#)pakgrd.f	5.2 11/28/95
      subroutine pakgrd(grid,wts,vwts,mxgrd,nugrd,nuwts,nuvwts,
     $                  newgrd,natoms)
c***begin prologue     pakgrd.f
c***date written       940201   (yymmdd) 
c***revision date      11/28/95      
c
c***keywords           
c***author             martin, richard(lanl)
c***source             @(#)pakgrd.f	5.2   11/28/95
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       pakgrd.f
      implicit none
c     --- input variables -----
      integer mxgrd,newgrd,natoms
c     --- input arrays (unmodified) ---
      real*8 grid(mxgrd,3,natoms),wts(mxgrd,natoms),vwts(mxgrd,natoms)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 nugrd(newgrd,3,natoms),nuwts(newgrd,natoms)
      real*8 nuvwts(newgrd,natoms)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer at,j
c
c
      do 120 at=1,natoms
         call vmove(nuwts(1,at),wts(1,at),newgrd)
         call vmove(nuvwts(1,at),vwts(1,at),newgrd)
         do 110 j=1,3
            call vmove(nugrd(1,j,at),grid(1,j,at),newgrd)
  110    continue
  120 continue
c
c
      return
      end
