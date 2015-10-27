*deck %W%  %G%
      subroutine readfx(frcnst,n)
c***begin prologue     %M%
c***date written       yymmdd  
c***revision date      %G%      
c
c***keywords           
c***author             page, michael(nrl)
c***source             %W%   %G%
c***purpose            reads force constant matrix from chk file.
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
      integer n
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 frcnst(n)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      character*8 namchk
c
      common /io/inp,iout
c
      call iosys('read character "checkpoint filename" from rwf',
     $            0,0,0,namchk)
      call iosys('open chk as old',0,0,0,namchk)
      call iosys('read real "cartesian second derivatives"from chk',
     $            n,frcnst,0,' ')
      call iosys('close chk',namchk,0,0,' ')
c
c
      return
      end
