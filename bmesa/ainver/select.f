*deck select.f
c***begin prologue     select
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
c***end prologue       select
      subroutine select(p,dp,ddp,q,wt,n,efac,qdtyp)
      implicit integer (a-z)
      real*8 p, dp, ddp, q, wt, fac, efac, sumwt
      character*80 title
      character*(*) qdtyp
      dimension p(n,n), dp(n,n), ddp(n,n)
      dimension q(n), wt(n) 
      common/io/inp, iout
      sumwt=0.d0
      do 10 i=1,n
         sumwt = sumwt + wt(i)
         fac=1.d0/sqrt(wt(i)) 
         call smul(p(1,i),p(1,i),fac,n)
         call smul(dp(1,i),dp(1,i),fac,n)
         call smul(ddp(1,i),ddp(1,i),fac,n)
 10   continue
c       write(iout,1) sumwt
c       title='points'
c       call prntfm(title,q,n,1,n,1,iout)
c       title='weights'
c       call prntfm(title,wt,n,1,n,1,iout)
       if(qdtyp.eq.'hermite') then
          fac=.5d0/(efac*efac)
          do 20 i=1,n
             do 30 j=1,n
                ddp(j,i) = ddp(j,i) - 4.d0*fac*q(j)*dp(j,i)
     1                              + 4.d0*fac*fac*q(j)*q(j)*p(j,i)
                dp(j,i) = dp(j,i) - 2.d0*fac*q(j)*p(j,i)
 30          continue
 20       continue
       endif                         
      return
 1    format(/,5x,'sum of the weights = ',e15.8)      
      end       
