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
      subroutine fpsi(hbuf,ibuf,diag,vnl,psi,hpsi,error,n,n1,n2,n3,
     1                lenbuf,nel,dim,incore,itdiag,prnt)
      implicit integer (a-z)
      real*8 hbuf, diag, vnl, psi, hpsi
      real*8 eij, error
      character*3 itoc
      character*1 it
      character*80 title
      logical prnt, itdiag, incore
      dimension hbuf(*), ibuf(*), diag(*), vnl(*), psi(n), hpsi(n)
      common/io/inp, iout
c
c     calculate hamiltonian on psi
c
      if(itdiag) then
         it=itoc(dim) 
         call honv(hbuf,ibuf,diag,psi,hpsi,n,1,lenbuf,nel,
     1             it,incore)
         call addvnl(psi,hpsi,vnl,n,1)
      else
         call ebc(hpsi,hbuf,psi,n,n,1)
         do 10 i=1,n
            hpsi(i) = hpsi(i) + vnl(i)*psi(i)
 10      continue            
      endif
      error=0.d0
      do 20 i=1,n
         do 30 j=1,i
            eij = hpsi(i)*psi(j) - hpsi(j)*psi(i)
            error = max(error,abs(eij))
 30      continue
 20   continue   
      return 
      end       

