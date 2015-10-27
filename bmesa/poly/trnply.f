*deck trnply.f
c***begin prologue     trnply
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           surface functions
c***author             schneider, barry (nsf)
c***source             math
c***purpose            
c***                   
c***references         
c
c***routines called    
c***end prologue       trnply
      subroutine trnply(p,ham,srf,scr,npts,n,nrt)
      implicit integer (a-z)
      real*8 p, ham, srf, scr
      character*80 title
      dimension p(npts,n), ham(n,nrt), srf(nrt), scr(n)
      common/io/inp, iout
      do 10 i=1,n
         scr(i)=p(npts,i)
   10 continue
      call ebtc(srf,ham,scr,nrt,n,1)          
c      title='surface functions'
c      call prntrm(title,srf,nrt,1,nrt,1,iout)
      return
      end       
