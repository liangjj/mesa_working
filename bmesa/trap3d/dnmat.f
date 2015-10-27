*deck dnmat.f
c***begin prologue     dnmat
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           non-linear schroedinger equation, gross-pitaevski
c***author             schneider, barry (nsf)
c***source             trap3d
c***purpose            form current approximation to density matrix
c***                   from current wavefunction
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       dnmat
      subroutine dnmat(rho,psiin,n)
      implicit integer (a-z)
      real*8 rho, psiin
      dimension rho(n,n), psiin(n)
      common/io/inp, iout
c
c     calculate density matrix.      
c
      do 10 i=1,n
         do 20 j=1,i
            rho(i,j) = psiin(i)*psiin(j)
            rho(j,i) = rho(i,j)
 20      continue
 10   continue            
      end       
