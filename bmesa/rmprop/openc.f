*deck openc 
c***begin prologue     openc
c***date written       941019   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6236, link 6236, open channels
c***author             schneider, b. i.(nsf)
c***source             m6236
c***purpose            determine number of open and closed channels 
c***                   at a given total energy
c***references
c
c***routines called    iosys, util and mdutil
c***end prologue       m6236
      subroutine openc(ec,energy,nc,nopen,nclosd)
      implicit integer (a-z)
      real*8 ec, energy
      character*8 val
      dimension ec(nc)
      common /io/ inp, iout
      write (iout,1) energy
      write (iout,2)
      nopen=0
      do 10 i=1,nc
         teste=energy-ec(i)
         val='closed'
         if (teste.ge.0.d0) then
             nopen=nopen+1
             val='open'
         endif
         write (iout,3) i, ec(i), val
   10 continue
c
      nclosd=nc-nopen
      write (iout,*)
      return
 1    format(/,15x,'total energy = ',e15.8)  
 2    format(/,1x,'channel',5x,'channel energy',5x,'open/closed',/)
 3    format(3x,i3,5x,e15.8,9x,a8)
      end
