*deck wrdadd.f
c***begin prologue     wrdadd
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
c***end prologue       wrdadd
 1    subroutine wrdadd(q,wt,eigc,wtc,p,dp,ddp,pn,dpn,ddpn,ham,eig,v,u,
     1                  tpr,dim,words)
      implicit integer (a-z)
      common/io/inp, iout
      ibeg=words
      q=ibeg
      wt=q+dim
      eigc=wt+dim
      wtc=eigc+dim
      p=wtc+dim
      dp=p+dim*dim
      ddp=dp+dim*dim
      pn=ddp+dim*dim
      dpn=pn+dim*dim
      ddpn=dpn+dim*dim
      ibeg=ddpn+dim*dim
      words=ibeg
      ham=words
      eig=ham+dim*dim
      v=eig+dim
      u=v+dim*dim
      tpr=u+dim*dim
      words=tpr+dim*dim
      return
      end       
