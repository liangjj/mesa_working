*deck guess.f
c***begin prologue     guess
c***date written       960615   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           hamiltonian
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            guess vectors.
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       guess
      subroutine guess(diag,dinv,vec,energy,n,ntrial)
      implicit integer (a-z)
      real*8 diag, dinv, vec, energy, tmp
      dimension diag(n), dinv(n), vec(n,ntrial)
      common/io/inp, iout
      do 10 i=1,ntrial
         do 20 j=1,n
            vec(j,i) = dinv(j)*vec(j,i)
 20      continue   
 10   continue   
      return
      end       






