*deck plymat.f
c***begin prologue     plymat
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            lobatto matrix elements.
c***                   
c***references         
c
c***routines called    
c***end prologue       plymat
      subroutine plymat(q,wt,p,dp,ddp,ov,fder,sder,fdfd,bloch,
     1                  n,coord,prn)
      implicit integer (a-z)
      real*8 q, wt, p, dp, ddp, ov, fder, sder, fdfd, bloch
      dimension q(n), wt(n), p(n,n), dp(n,n), ddp(n,n), ov(n,n)
      dimension fder(n,n), sder(n,n), fdfd(n,n), bloch(n,n)
      character*80 title
      character*(*) coord
      character*16 typarg
      logical prn
      common/io/inp, iout
      typarg='linear'
      if(coord.eq.'cylindrical') then
         typarg='quadratic'
      endif
      call fafb(p,p,wt,ov,n)
      call fadfb(p,dp,q,wt,fder,typarg,n) 
      call faddfb(p,dp,ddp,q,wt,sder,typarg,n)  
      call dfadfb(dp,dp,q,wt,fdfd,typarg,n)
      call tbloch(p,dp,q,bloch,typarg,n)  
      if(prn) then
         title='overlap matrix'
         call prntrm(title,ov,n,n,n,n,iout)
         title='psi*dpsi matrix'
         call prntrm(title,fder,n,n,n,n,iout)
         title='psi*ddpsi matrix'
         call prntrm(title,sder,n,n,n,n,iout)
         title='dpsi*dpsi matrix'
         call prntrm(title,fdfd,n,n,n,n,iout)
         title='bloch matrix'
         call prntrm(title,bloch,n,n,n,n,iout)
         call maket(p,dp,ddp,q,wt,sder,bloch,fdfd,typarg,n)
         title='symmetrized kinetic energy matrix'
         call prntrm(title,fdfd,n,n,n,n,iout)
      endif
      return
      end       






