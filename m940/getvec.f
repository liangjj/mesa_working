*deck getvec.f
c***begin prologue     getvec
c***date written       960615   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           hamiltonian
c***author             schneider, barry (nsf)
c***source             
c***purpose            read in p space vectors
c***                   
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       getvec
      subroutine getvec(pvec,temp,n,nwks,ntot)
c
      implicit integer (a-z)
c
      real*8 pvec, temp
      character*4 itoc
      character*3 ans
      character*80 title
      dimension pvec(nwks,*), temp(n)
c
      common /io/ inp,iout
c
      call iosys('does "p-space vectors" exist on hpart',0,0,0,ans)
      offset=n-nwks+1
      if(ans.eq.'yes') then
         do 10 i=1,ntot
            call iosys('read real "p-space vectors" from hpart '//
     1                 'without rewinding',n,temp,0,' ')
            call copy(temp(offset),pvec(1,i),nwks)
 10      continue   
      else
         call rzero(pvec,ntot*nwks)
         do 20 i=1,ntot
            pvec(i,i)=1.d0
 20      continue
      endif   
c
      title='new vector set'
      call prntrm(title,pvec,nwks,ntot,nwks,ntot,iout)
      return
      end




