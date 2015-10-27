*deck @(#)readfx.f	5.1  11/6/94
      subroutine readfx(frcnst,n)
c***begin prologue     readfx.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             page, michael(nrl)
c***source             @(#)readfx.f	5.1   11/6/94
c***purpose            reads force constant matrix from chk file.
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       readfx.f
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
