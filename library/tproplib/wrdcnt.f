*deck wrdcnt.f
c***begin prologue     wrdcnt
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            memory outlay for a particular dimension and grid set.
c***                   
c***references         
c
c***routines called    
c***end prologue       wrdcnt
      subroutine wrdcnt(q,wt,eigc,wtc,p,dp,ddp,pn,dpn,ddpn,ham,eig,v,u,
     1                  tpr,dim,words)
      implicit integer (a-z)
      common/io/inp, iout
      q=words
      wt=q+dim
      eigc=wt+dim
      wtc=eigc+dim
      p=wtc+dim
      dp=p+dim*dim
      ddp=dp+dim*dim
      pn=ddp+dim*dim
      dpn=pn+dim*dim
      ddpn=dpn+dim*dim
      ham=ddpn+dim*dim
      eig=ham+dim*dim
      v=eig+dim
      u=v+dim*dim
      tpr=u+dim*dim
      words=tpr+dim*dim
      return
      end       
