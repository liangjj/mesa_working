*deck @(#)masswt.f	5.1  11/6/94
      subroutine masswt(fc,cmass,natoms3)
c
c***begin prologue     masswt
c***date written       850601  
c***revision date      11/6/94      
c
c***keywords           
c***author             page, michael 
c***source             @(#)masswt.f	5.1   11/6/94
c***purpose            to mass weight the force constant matrix. 
c***description
c     
c    
c
c***references
c
c***routines called
c     %m%
c       start here
c
c***end prologue       masswt
c
      implicit integer(a-z)
      real*8 fc(natoms3,natoms3),cmass(natoms3)
c
      do 40 j=1,natoms3
         do 50 k=1,natoms3
            fc(j,k)=fc(j,k)/(sqrt(cmass(j))*sqrt(cmass(k)))
   50    continue
   40 continue
c
      return
      end
