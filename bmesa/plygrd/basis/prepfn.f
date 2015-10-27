*deck prepfn.f
c***begin prologue     prepfn
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           polynomials
c***author             schneider, barry (nsf)
c***source             
c***purpose            prepare orthonormal basis set from polynomials
c***                   
c***                                                          
c***references         
c
c***routines called    
c***end prologue       prepfn
      subroutine prepfn(p,dp,ddp,wt,ni,nj)
      implicit integer (a-z)
      real*8 p, dp, ddp, wt, fac, sumwt
      dimension p(nj,ni), dp(nj,ni), ddp(nj,ni)
      dimension wt(ni) 
      common/io/inp, iout
      sumwt=0.d0
      do 10 i=1,ni
         sumwt = sumwt + wt(i)
         fac=1.d0/sqrt(wt(i)) 
         call smul(p(1,i),p(1,i),fac,nj)
         call smul(dp(1,i),dp(1,i),fac,nj)
         call smul(ddp(1,i),ddp(1,i),fac,nj)
 10   continue
      write(iout,1) sumwt
      return
 1    format(/,5x,'sum of the weights = ',e15.8)      
      end       
