*deck tthet.f
c***begin prologue     tthet
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           hamiltonian, one-particle
c***author             schneider, barry (nsf)
c***source             
c***purpose            one particle hamiltonian, eigenvalues
c***                                and eigenvectors.
c***                   
c***description        kinetic energy operator in theta coordinate 
c***                   1.        d     [              d   ]
c***                 ----       ----   [ sin(theta) ----  ]
c***                sin(theta)  theta  [            theta ]
c***references         
c
c***routines called    
c***end prologue       tthet
      subroutine tthet(p,dp,ddp,q,wt,t,tmp,ilft,irt,n)
      implicit integer (a-z)
      real*8 p, dp, ddp, q, wt, t
      dimension p(n,0:n-1), dp(n,0:n-1), ddp(n,0:n-1)
      dimension q(n), wt(n), t(n,n), tmp(n,2)
      common/io/inp, iout
      do 10 i=1,n
         tmp(i,1) = sin(q(i))
         tmp(i,1) = cos(q(i))
 10   continue   
      do 20 i=1,n
         do 30 j=1,i
c            do 40 k=1,n
c               t(i,j) = t(i,j) - wt(k)*p(k,i-1)*( tmp(k,1)*ddp(k,j-1) 
c     1                                          + tmp(k,2)*dp(k,j-1) )
c 40         continue
            t(i,j) = t(i,j) - wt(i)*p(i,i-1)*( tmp(i,1)*ddp(i,j-1) 
     1                                       + tmp(i,2)*dp(i,j-1) )
            t(i,j) = t(i,j) + irt*tmp(n,1)*p(n,i-1)*dp(n,j-1) 
     1                      - ilft*tmp(1,1)*p(1,i-1)*dp(1,j-1)
            t(j,i) = t(i,j)
 30      continue   
 20   continue
      return
      end       
