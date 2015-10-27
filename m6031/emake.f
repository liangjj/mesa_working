*deck emake
      subroutine emake(energy,e0,dele,nen)
      implicit integer (a-z)
      real*8 energy, e0, dele
      dimension energy(nen)
      energy(1)=e0
      do 10 ien=2,nen
         energy(ien)=energy(ien-1)+dele 
   10 continue            
      return
      end
