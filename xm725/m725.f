*deck %W%  %G%
      program m725
c***begin prologue     m725.f
c***date written       yymmdd  
c***revision date      2/5/95      
c
c***keywords           
c***author             
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
c***end prologue       m725.f
      implicit none
c     --- input variables -----
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer ia,ioff,intoff,wpadti,maxcor
      real*8 z(1)
c
      common // ia(1)
      common /memory/ ioff
      equivalence (ia(1),z(1))
      integer nodeid
c
c
      call pbeginf
C      call niceftn(10)
      if (nodeid() .eq. 0) call drum
c      call llog
      call getscm(0,z,maxcor,'m725',0)
      intoff=wpadti(ioff)
      call pm725(z(ioff),ia(intoff))
      call pend
c
c
      stop
      end
