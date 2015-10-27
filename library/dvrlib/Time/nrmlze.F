*deck nrmlze.f
c***begin prologue     nrmlze
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            normalize lobatto polynomials and their 
c**                    first derivatives.
c***                   
c***references         
c
c***routines called    
c***end prologue       nrmlze
      subroutine nrmlze(p,dp,wt,n)
      implicit integer (a-z)
      real*8 wt, p, dp, nrm
      dimension wt(n), p(n,n), dp(n,n)
      common/io/inp, iout
      do 10 i=1,n
         nrm=1.d0/sqrt(wt(i))
         call vscale(p(1,i),p(1,i),nrm,n)
         call vscale(dp(1,i),dp(1,i),nrm,n)
 10   continue   
      return
      end       


