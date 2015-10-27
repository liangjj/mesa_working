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
 1    subroutine wrdcnt(q,wt,eigc,wtc,p,dp,ddp,pn,dpn,ddpn,ham,eig,v,u,
     1                  tpr,tpra,mattyp,dim,ngrid,words)
      implicit integer (a-z)
      character*(*) mattyp
      dimension q(ngrid), wt(ngrid), eigc(ngrid), wtc(ngrid)
      dimension p(ngrid), dp(ngrid), ddp(ngrid), pn(ngrid), dpn(ngrid)
      dimension ddpn(ngrid)
      common/io/inp, iout
      ibeg=words
      do 10 i=1,ngrid
         q(i)=ibeg
         wt(i)=q(i)+dim
         eigc(i)=wt(i)+dim
         wtc(i)=eigc(i)+dim
         p(i)=wtc(i)+dim
         dp(i)=p(i)+dim*dim
         ddp(i)=dp(i)+dim*dim
         pn(i)=ddp(i)+dim*dim
         dpn(i)=pn(i)+dim*dim
         ddpn(i)=dpn(i)+dim*dim
         ibeg=ddpn(i)+dim*dim
 10   continue
       words=ibeg
       if(mattyp.eq.'complex'.or.mattyp.eq.'real-unsymmetric') then
          ham=words
          eig=ham+2*dim*dim
          v=eig+2*dim
          u=v+2*dim*dim
          tpr=u+2*dim*dim
          tpra=tpr+2*dim*dim
          words=tpra+2*dim*dim
        elseif(mattyp.eq.'real-symmetric') then               
          ham=words
          eig=ham+dim*dim
          v=eig+dim
          u=v+dim*dim
          tpr=u+dim*dim
          tpra=tpr+dim*dim
          words=tpra+dim*dim
        elseif(mattyp.eq.'real-unsymmetric') then
          ham=words
          eig=ham+dim*dim
          v=eig+2*dim
          u=v+dim*dim
          tpr=u+dim*dim
          tpra=tpr+dim*dim
          words=tpra+dim*dim 
        endif                     
      return
      end       
