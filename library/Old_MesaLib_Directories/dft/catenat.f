*deck @(#)catenat.f	5.2 2/5/95
      subroutine catenat(nrad,mxgrd,mxleb,nleb,rpts,rwt,angpts,
     $                  angwt,grid,gridwts,ptrad,atcen)
c***begin prologue     catenat.f
c***date written       940601  
c***revision date      2/5/95      
c
c***keywords           
c***author             martin, r.l. 
c***source             @(#)catenat.f	5.2   2/5/95
c***purpose            concatenates radial and angular points into grid.
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       catenat.f
      implicit none
c     --- input variables -----
      integer nrad,nleb,mxleb,mxgrd
c     --- input arrays (unmodified) ---
      real*8 rpts(nrad),rwt(nrad)
      real*8 angpts(mxleb,3),angwt(mxleb)
      real*8 atcen(3)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 grid(mxgrd,3),gridwts(mxgrd)
      integer ptrad(0:nrad)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer radshel,lebpt,gpt
      integer inp,iout
c
      common/io/inp,iout
c
c     --- concatenate the radial and angular points into grid
      gpt=0
      ptrad(0)=1
      do 10 radshel=1,nrad
         ptrad(radshel)=ptrad(radshel-1)+nleb
         do 20 lebpt=1,nleb
            gpt=gpt+1
            grid(gpt,1)=angpts(lebpt,1)*rpts(radshel)+atcen(1)
            grid(gpt,2)=angpts(lebpt,2)*rpts(radshel)+atcen(2)
            grid(gpt,3)=angpts(lebpt,3)*rpts(radshel)+atcen(3)
            gridwts(gpt)=angwt(lebpt)*rwt(radshel)
   20    continue 
   10 continue 
c
c
      if(gpt.gt.mxgrd) then
         write(iout,*) 'grid points larger than allowed'
         write(iout,*) 'gpt=',gpt,' mxgrd=',mxgrd
         call lnkerr('catenat')
      endif
c
c
      return
      end
