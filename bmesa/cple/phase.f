*deck phase
c***begin prologue     phase
c***date written       940213   (yymmdd)
c***revision date               (yymmdd)
c***keywords           s-wave, phase
c***author             schneider, barry (nsf)
c***source             
c***purpose            s-wave phase shifts
c***description        
c***references       
c
c***routines called
c***end prologue       phase
      subroutine phase(eig,flst,energy,rbox,nen,n)
      implicit integer (a-z)
      common /io/ inp, iout
      real*8 eig, flst, energy, rbox, rmat, k, sn, cn, tandel, pse
      dimension eig(n), flst(n), energy(nen)
      write(iout,1)
      do 10 i=1,nen
         rmat=0.d0
         do 20 j=1,n
            rmat=rmat+flst(j)*flst(j)/(eig(j)-energy(i))
   20    continue
         rmat=.5d0*rmat
         k=sqrt(2.d0*energy(i))
         sn=sin(k*rbox)
         cn=cos(k*rbox)
         tandel=(rmat*k*cn-sn)/(cn+rmat*k*sn)
         pse=atan(tandel)
         write(iout,2) energy(i), k, rmat, tandel, pse         
   10 continue                     
    1 format(/,1x,'      energy     ','    k value    ',
     1            '   r-matrix    ','   tan phase   ',
     2            '    phase     ')
    2 format(1x,5e15.8)                          
      return
      end
