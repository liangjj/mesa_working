*deck @(#)getgrd.f	5.1  11/6/94
      subroutine getgrd(gcm,mass,natoms,natoms3)
c***begin prologue     getgrd.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             page, michael(nrl)
c***source             @(#)getgrd.f	5.1   11/6/94
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       getgrd.f
      implicit none
c     --- input variables -----
      integer natoms,natoms3
c     --- input arrays (unmodified) ---
      real*8 mass(natoms)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 gcm(natoms3)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer atom
c
c     --- read first derivatives
         call iosys('read real "cartesian first derivatives" from rwf',
     $               natoms3,gcm,0,' ')
         do 10 atom=1,natoms
            gcm(3*atom-2)=gcm(3*atom-2)/sqrt(mass(atom))
            gcm(3*atom-1)=gcm(3*atom-1)/sqrt(mass(atom))
            gcm(3*atom  )=gcm(3*atom  )/sqrt(mass(atom))
   10    continue
c
c
         return
         end
