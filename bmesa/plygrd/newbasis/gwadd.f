*deck gwadd.f
c***begin prologue     gwadd
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            memory outlay for a particular grid.
c***                   
c***references         
c
c***routines called    
c***end prologue       gwadd
 1    subroutine gwadd(q,wt,eigc,wtc,p,dp,ddp,pn,dpn,ddpn,
     1                 npt,nsubg,maxgrd,words)
      implicit integer (a-z)
      dimension q(nsubg), wt(nsubg), eigc(nsubg), wtc(nsubg), npt(nsubg)
      dimension p(maxgrd,maxgrd), dp(maxgrd,maxgrd), ddp(maxgrd,maxgrd)
      dimension pn(maxgrd,maxgrd), dpn(maxgrd,maxgrd)
      dimension ddpn(maxgrd,maxgrd)
      common/io/inp, iout
      ibeg=words
      do 10 i=1,nsubg
         q(i)=ibeg
         wt(i)=q(i)+npt(i)
         eigc(i)=wt(i)+npt(i)
         wtc(i)=eigc(i)+npt(i)
         ibeg=wtc(i)+npt(i)
 10   continue
      do 20 i=1,nsubg   
         do 30 j=1,nsubg
            p(j,i)=ibeg
            dp(j,i)=p(j,i)+npt(i)*npt(j)
            ddp(j,i)=dp(j,i)+npt(i)*npt(j)
            pn(j,i)=ddp(j,i)+npt(i)*npt(j)
            dpn(j,i)=pn(j,i)+npt(i)*npt(j)
            ddpn(j,i)=dpn(j,i)+npt(i)*npt(j)
            ibeg=ddpn(j,i)+npt(i)*npt(j)
 30      continue
 20   continue   
      words=ibeg
      return
      end       
