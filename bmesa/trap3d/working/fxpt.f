*deck fxpt.f
c***begin prologue     fxpt
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
c***end prologue       fxpt
      subroutine fxpt(nfix,nreg,card)
      implicit integer (a-z)
      logical nfix
      character*(*) card
      character*2 itoc
      dimension nfix(2,nreg)
      common/io/inp, iout
      totfx=0
      do 10 i=1,nreg
         nfix(1,i)=logkey(card,'fix-left-end-point-region-'//itoc(i),
     1                    .false.,' ')
         nfix(2,i)=logkey(card,'fix-right-end-point-region-'//itoc(i),
     1                    .false.,' ')
         if(nfix(1,i)) then
            totfx=totfx+1
         endif
         if(nfix(2,i)) then
            totfx=totfx+1
         endif
 10   continue
c      write(iout,1) totfx
      return
 1    format(/,1x,'total number of fixed points = ',i2)
      end       
