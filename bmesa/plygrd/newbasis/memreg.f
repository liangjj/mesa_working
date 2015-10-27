*deck memreg.f
c***begin prologue     memreg
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            memory outlay for a particular region.
c***                   
c***references         
c
c***routines called    
c***end prologue       memreg
 1    subroutine memreg(pptwt,ppoly,ngrid,npt,ngot,maxgr,maxd)
      implicit integer (a-z)
      real*8 ptwt, poly
      dimension npt(ngrid)
      common/io/inp, iout
      pointer (pptwt,ptwt(1))
      pointer (ppoly,poly(1))
      maxgr=max(maxgr,ngrid)
      words=1
      do 10 grd=1,ngrid
         maxd=max(maxd,npt(grd))
         q=words
         wt=q+npt(grd)
         words=wt+npt(grd)
 10   continue
      words=wpadti(words)
      call memory(words,pptwt,ngot(1),'grid',0)
      words=1   
      do 20 grdi=1,ngrid
         do 30 grdj=1,ngrid
            p=words
            dp=p+npt(grdi)*npt(grdj)
            ddp=dp+npt(grdi)*npt(grdj)
            words=ddp+npt(grdi)*npt(grdj)
 30      continue
 20   continue   
      words=wpadti(words)
      call memory(words,ppoly,ngot(2),'basis',0)
      return
      end       
