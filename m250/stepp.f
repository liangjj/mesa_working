*deck @(#)stepp.f	5.1  11/6/94
      subroutine stepp(gc,delx,cref,trvec,mass,natoms,
     $                 nat3,projg,sstep,phase)
c***begin prologue     stepp.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             page, michael(nrl)
c***source             @(#)stepp.f	5.1   11/6/94
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       stepp.f
      implicit none
c     --- input variables -----
      integer natoms,nat3
      logical projg
      real*8 sstep,phase
c     --- input arrays (unmodified) ---
      real*8 mass(natoms)
      real*8 gc(nat3)
      real*8 trvec(nat3)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 delx(nat3),cref(nat3)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer i,d2ec
      real*8 gcmnrm,zero
c
      parameter (zero=0.0d+00)
c
      common /io/ inp,iout
c
      if(projg) then
c
c        --- take a step along the negative gradient
         gcmnrm=zero
         do 10 i=1,nat3
            gcmnrm=gcmnrm+gc(i)**2
   10    continue
         write(iout,1000)
 1000    format(/' stepping along the gradient')
         gcmnrm=sqrt(gcmnrm)
         do 20 i=1,nat3
            delx(i)=-(sstep/gcmnrm)*gc(i)
   20    continue
      else
c
c        --- step along the eigenvector
         write(iout,2000)
 2000    format(/' stepping along the eigenvector ')
         write(iout,3030) (trvec(i),i=1,nat3)
 3030    format(/' eigenvector : '/,(15x,f8.4))
c
         do 30 i=1,nat3
            delx(i)=phase*sstep*trvec(i)
   30    continue
         d2ec=1
         call iosys('write integer d2e_coord to rwf',1,d2ec,0,' ')
      endif
c
c     --- record coordinates
      write(iout,3010) (cref(i),i=1,nat3)
 3010 format(//' unperturbed coordinates: ',/,(3f16.6))
c
      write(iout,3000) (delx(i),i=1,nat3)
 3000 format(//' step: ',/,(3f16.6))
      do 50 i=1,nat3
         cref(i)=cref(i)+delx(i)
   50 continue
c
      write(iout,3020) (cref(i),i=1,nat3)
 3020 format(//' new coordinates: ',/,(3f16.6))
c
c     --- store the new coordinates
      call iosys('write real coordinates to rwf',
     $            nat3,cref,0,' ')
c
c
      return
      end
