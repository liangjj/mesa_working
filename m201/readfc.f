*deck @(#)readfc.f	5.1  11/6/94
      subroutine readfc(frcnst,n,file)
c***begin prologue     readfc.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             
c***source             @(#)readfc.f	5.1   11/6/94
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       readfc.f
      implicit none
c     --- input variables -----
      integer n
c     --- input arrays (unmodified) ---
      character*(*) file
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 frcnst(n)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      character*8 namchk
      character*4 filenm
      integer inp,iout
c
      common /io/inp,iout
c
c     --- open up the appropriate file
      filenm=file
      if(file.eq.'rwf') then
         call iosys('read real force_constants from rwf',n,frcnst,0,' ')
      else if (file.eq.'chk') then
         call iosys('read character "checkpoint filename" from rwf',
     $              0,0,0,namchk)
         call iosys('open chk as old',0,0,0,namchk)
         call iosys('read real force_constants from chk',n,frcnst,0,' ')
         call iosys('close chk',namchk,0,0,' ')
      endif
c
c
      return
      end
