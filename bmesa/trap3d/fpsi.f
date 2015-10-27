*deck fpsi.f
c***begin prologue     fpsi
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           non-linear schroedinger equation, gross-pitaevski
c***author             schneider, barry (nsf)
c***source             tstdiis
c***purpose            form current approximation to the value of
c***                   the fock operator on the wavefunction.  also
c***                   return the maximum value of [f,rho] which
c***                   is a measure of the error.
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       fpsi
      subroutine fpsi(f,psi,hpsi,error,n,iter,prnt)
      implicit integer (a-z)
      real*8 f, psi, hpsi
      real*8 eij, error
      character*3 itoc
      character*80 title
      logical prnt
      dimension f(n,n), psi(n), hpsi(n)
      common/io/inp, iout
c
c     calculate hamiltonian on psi
c
      call ebc(hpsi,f,psi,n,n,1)
      error=0.d0
      do 10 i=1,n
         do 20 j=1,i
            eij = hpsi(i)*psi(j) - hpsi(j)*psi(i)
            error = max(error,abs(eij))
 20      continue
 10   continue   
      return 
      end       
