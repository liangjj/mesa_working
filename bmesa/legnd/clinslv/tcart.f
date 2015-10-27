*deck tcart.f
c***begin prologue     tcart
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           hamiltonian, one-particle
c***author             schneider, barry (nsf)
c***source             
c***purpose            one particle hamiltonian, eigenvalues
c***                                and eigenvectors.
c***                   
c***description        kinetic energy operator in cartesian coordinates 
c***references         
c
c***routines called    
c***end prologue       tcart
      subroutine tcart(p,dp,ddp,wt,t,ilft,irt,n,nonsym)
      implicit integer (a-z)
      real*8 p, dp, ddp, wt
      complex*16 t
      dimension p(n,0:n-1), dp(n,0:n-1), ddp(n,0:n-1)
      dimension wt(n), t(n,n)
      common/io/inp, iout
      do 10 i=1,n
         do 20 j=1,n
            do 30 k=1,n
               t(i,j) = t(i,j) - wt(k)*p(k,i-1)*ddp(k,j-1)
 30         continue
            t(i,j) = t(i,j) + irt*p(n,i-1)*dp(n,j-1) 
     1                      - ilft*p(1,i-1)*dp(1,j-1)
 20      continue   
 10   continue
      return
      end       
