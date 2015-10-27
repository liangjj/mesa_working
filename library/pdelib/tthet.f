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
      subroutine tthet(p,dp,ddp,q,wt,t,ct,tmp,ilft,irt,n,mattyp)
      implicit integer (a-z)
      character*(*) mattyp
      real*8 p, dp, ddp, q, wt, t, tmp
      complex*16 ct
      dimension p(n,0:n-1), dp(n,0:n-1), ddp(n,0:n-1)
      dimension q(n), wt(n), t(n,n), ct(n,n), tmp(n,2)
      common/io/inp, iout
      do 10 i=1,n
         tmp(i,1) = sin(q(i))
         tmp(i,1) = cos(q(i))
 10   continue
      if(mattyp.eq.'complex'.or.mattyp.eq.'real-unsymmetric') then
         do 20 i=1,n
            do 30 j=1,n
               ct(i,j) = ct(i,j) - wt(i)*p(i,i-1) * 
     1                             ( tmp(i,1)*ddp(i,j-1) 
     2                                       + tmp(i,2)*dp(i,j-1) )
               ct(i,j) = ct(i,j) + irt*tmp(n,1)*p(n,i-1)*dp(n,j-1) 
     1                           - ilft*tmp(1,1)*p(1,i-1)*dp(1,j-1)
 30         continue   
 20      continue
      else
         do 40 i=1,n
            do 50 j=1,i
               t(i,j) = t(i,j) - wt(i)*p(i,i-1)*( tmp(i,1)*ddp(i,j-1) 
     1                                         + tmp(i,2)*dp(i,j-1) )
               t(i,j) = t(i,j) + irt*tmp(n,1)*p(n,i-1)*dp(n,j-1) 
     1                         - ilft*tmp(1,1)*p(1,i-1)*dp(1,j-1)
               t(j,i) = t(i,j)
 50         continue   
 40      continue
      endif
      return
      end       
