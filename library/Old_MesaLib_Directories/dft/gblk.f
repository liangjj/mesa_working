*deck @(#)gblk.f	5.1 11/6/94
      subroutine gblk(natoms,mxgbsiz,mxgblk,ngrid,ngb,gblksiz)
c***begin prologue     gblk.f
c***date written       940205  
c***revision date      11/6/94      
c
c***keywords           
c***author             
c***source             @(#)gblk.f	5.1   11/6/94
c***purpose            determines the number and sizes of each
c                      atomic grid block
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       gblk.f
      implicit none
c     --- input variables -----
      integer natoms,mxgbsiz,mxgblk
c     --- input arrays (unmodified) ---
      integer ngrid(natoms)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      integer ngb(natoms),gblksiz(mxgblk,natoms)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer iat,igb
      integer mod,left
      integer inp,iout
c
      common/io/inp,iout
c
      do 120 iat=1,natoms
         ngb(iat)=ngrid(iat)/mxgbsiz
         if(mod(ngrid(iat),mxgbsiz).ne.0) ngb(iat)=ngb(iat)+1
         if (ngb(iat).gt. mxgblk) then 
            write(iout,*)' atom ',iat,' needs ',ngrid(iat),' points'
            write(iout,*)' mxgbsiz=',mxgbsiz,' needs ',ngb(iat),
     $           ' blocks'
            call lnkerr(
     $           'm611: too many grid blocks needed for one atom.'//
     $     '      use scf=maxgblk=<n> to increase above default of 5')
         else
            left=ngrid(iat)
            do 110 igb=1,ngb(iat)
               if(left.le.mxgbsiz) then
                  gblksiz(igb,iat)=left
               else
                  gblksiz(igb,iat)=mxgbsiz
                  left=left-mxgbsiz
               endif
  110       continue 
         endif
  120 continue
c
c
      return
      end
