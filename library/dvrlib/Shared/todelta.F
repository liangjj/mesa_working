*deck todelta.f
c***begin prologue     todelta
c***date written       022202   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           coordinate eigenfunctions
c***author             schneider, barry (nsf)
c***source
c***purpose            renormalization
c***           
c***description
c***references         
c
c***routines called    
c***end prologue       todelta
      subroutine todelta(p,dp,ddp,wt,n)
      implicit integer (a-z)
      real*8 p, dp, ddp, wt, norm
      dimension p(n,n), dp(n,n), ddp(n,n)
      dimension wt(n)
      common/io/inp, iout 
      do 10 i=1,n
         norm=sqrt(wt(i))
         do 20 j=1,n
            p(j,i)= p(j,i)*norm      
            dp(j,i)= dp(j,i)*norm      
            ddp(j,i)= ddp(j,i)*norm      
 20      continue
 10   continue   
      return
      end       
