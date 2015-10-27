*deck tcylin.f
c***begin prologue     tcylin
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           hamiltonian, one-particle
c***author             schneider, barry (nsf)
c***source             
c***purpose            one particle hamiltonian, eigenvalues
c***                                and eigenvectors.
c***                   
c***description        kinetic energy operator in rho coordinate 
c***                   2
c***                  d          1.
c***                 ----   + -------
c***                     2    
c***                 drho     4. * rho*rho
c***references         
c
c***routines called    
c***end prologue       tcylin
      subroutine tcylin(p,dp,ddp,q,wt,t,ilft,irt,n)
      implicit integer (a-z)
      real*8 p, dp, ddp, q, wt, t, fac
      dimension p(n,0:n-1), dp(n,0:n-1), ddp(n,0:n-1)
      dimension q(n), wt(n), t(n,n)
      common/io/inp, iout
      data fac / .25d0 /
      do 10 i=1,n
         do 20 j=1,i
c            do 30 k=1,n
c               t(i,j) = t(i,j) - wt(k)*p(k,i-1)*( q(k)*ddp(k,j-1) 
c     1                                         + dp(k,j-1) )
c 30         continue
            t(i,j) = t(i,j) - wt(i)*p(i,i-1)*( q(i)*ddp(i,j-1) 
     1                                      + dp(i,j-1) )
            t(i,j) = t(i,j) + irt*p(n,i-1)*q(n)*dp(n,j-1) 
     1                      - ilft*p(1,i-1)*q(1)*dp(1,j-1)
            t(j,i) = t(i,j)
 20      continue   
 10   continue
      return
      end       
