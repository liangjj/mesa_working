*deck relene.f
c***begin prologue     relene
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***                   
c***references         
c
c***routines called    
c***end prologue       relene
      subroutine relene(energy,eigt,nen,n)
      implicit integer (a-z)
      real*8 energy, eigt
      dimension energy(nen), eigt(n)
      common/io/inp, iout
      write(iout,1) (eigt(i),i=1,n)
      do 10 i=1,nen
         energy(i) = energy(i) + eigt(1)
 10   continue
      write(iout,2) (energy(i),i=1,nen)
 1    format(/,1x,'target energies',(/,1x,5e15.8))
 2    format(/,1x,'scattering energies',(/,1x,5e15.8))
      return      
      end       






