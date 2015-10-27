*deck %W%  %G%
      subroutine updorb(delx,cphf,cmo,nat3,nbf)
c***begin prologue     %M%
c***date written       yymmdd  
c***revision date      %G%      
c
c***keywords           
c***author             page, michael(nrl)
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
c***end prologue       %M%
      implicit none
c     --- input variables -----
      integer nat3,nbf
c     --- input arrays (unmodified) ---
      real*8 delx(nat3)
c     --- input arrays (scratch) ---
      real*8 cphf(nbf**2)
c     --- output arrays ---
      real*8 cmo(nbf**2)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer i
c
      common /io/ inp,iout
c
      call iosys('read real "mcscf vector" from rwf',
     $           nbf**2,cmo,0,' ')
c
c     --- loop over nuclear degrees of freedom
c         for each degree of freedom, read in the cphf solution
c         and add it to the mo's scaled by the nuclear displacement
      do 100 i=1,nat3
         call iosys(
     $     'read real "cphf solutions" from rwf without rewinding',
     $      nbf**2,cphf,0,' ')
         call vwxs(cmo,cmo,cphf,delx(i),1,nbf**2)
  100 continue
      call iosys('write real "mcscf vector" to rwf',
     $           nbf**2,cmo,0,' ')
c
c
      return
      end
