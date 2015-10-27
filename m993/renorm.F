*deck renorm.f 
c***begin prologue     renorm
c***date written       000702   (yymmdd)
c***revision date               (yymmdd)
c***keywords           dvr
c***                   
c***author             schneider, b. i.(nsf)
c***source             renorm
c***purpose            renormalize coordinate functions
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       renorm
      subroutine renorm(p,dp,ddp,wt,n)
c
      implicit integer (a-z)
      real*8 p, dp, ddp, wt, fac
      character*80 title
      dimension p(n,n), dp(n,n), ddp(n,n), wt(n)
      common/io/inp, iout
      do 10 i=1,n
         fac=sqrt ( 1.d0/(p(i,i)*p(i,i)*wt(i)) )
         call vscale(p(1,i),p(1,i),fac,n)
         call vscale(dp(1,i),dp(1,i),fac,n)
         call vscale(ddp(1,i),ddp(1,i),fac,n)
 10   continue
      title='renormalized functions'   
      call prntrm(title,p,n,n,n,n,iout)
      title='renormalized first derivatives'   
      call prntrm(title,dp,n,n,n,n,iout)
      title='renormalized second derivatives'   
      call prntrm(title,ddp,n,n,n,n,iout)
      return
      end


















