*deck dropfn.f
c***begin prologue     dropfn
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           polynomials
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***                   
c***                                                          
c***references         
c
c***routines called    
c***end prologue       dropfn
      subroutine dropfn(p,dp,ddp,pn,dpn,ddpn,q,wt,eigc,wtc,
     1                  lftbc,rtbc,nin,nout)
      implicit integer (a-z)
      real*8 p, dp, ddp, pn, dpn, ddpn, q, wt, eigc, wtc
      dimension p(nin,*), dp(nin,*), ddp(nin,*)
      dimension pn(nout,*), dpn(nout,*), ddpn(nout,*)
      dimension q(nin), wt(nin), eigc(nout), wtc(nout)
      common/io/inp, iout
      nn=nin
      nup=nin
      if(rtbc.eq.0) then
         nup=nin-1
      endif
      if(lftbc.eq.0) then
         one=1
         write(iout,1) one
         do 10 i=2,nup
            eigc(i-1)=q(i)
            wtc(i-1)=wt(i)
            call copy(p(2,i),pn(1,i-1),nout)
            call copy(dp(2,i),dpn(1,i-1),nout)
            call copy(ddp(2,i),ddpn(1,i-1),nout)
 10      continue   
         nn=nn-1
      endif
      if(rtbc.eq.0) then
         last=nin
         write(iout,1) last
         nn=nn-1
      endif
      nout=nn
      call copy(pn,p,nout*nout)
      call copy(dpn,dp,nout*nout)
      call copy(ddpn,ddp,nout*nout)
      return
 1    format(/,5x,'dropping function = ',i3,' from basis set')
      end       
