*deck newerr.f
c***begin prologue     newerr
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           non-linear schroedinger equation, gross-pitaevski
c***author             schneider, barry (nsf)
c***source             tstdiis
c***purpose            form diis approximation to error fock matrix
c***                   
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       newerr
      subroutine newerr(psi,hpsi,sol,error,n,iter)
      implicit integer (a-z)
      real*8 err, psi, hpsi, sol, error
      dimension psi(n,*), hpsi(n,*), sol(iter)
      common/io/inp, iout
      error=0.d0
      do 10 i=1,n
         do 20 j=1,i
            err = 0.d0
            do 30 k=1,iter
               err = err +sol(k)*( psi(i,k)*hpsi(j,k) 
     1                                     -
     2                             psi(j,k)*hpsi(i,k) )
 30         continue
            error=max(error,abs(err))
 20      continue
 10   continue   
      return 
      end       


