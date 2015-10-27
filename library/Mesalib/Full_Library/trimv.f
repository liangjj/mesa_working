*deck trimv
c***begin prologue     trimv
c***date written       920209   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           trim a vector
c***author             schneider, barry (lanl)
c***source             mylib 
c***                                         
c***purpose            trim components of a vector
c***                   
c***                  
c***references         none
c
c***routines called    
c***end prologue       trimv
      subroutine trimv(vin,vout,vscr,keep,nkeep,n)
      implicit integer (a-z)
      real*8 vin, vout, vscr
      dimension vin(n), vout(nkeep), vscr(nkeep), keep(nkeep)
      do 10 i=1,nkeep
         vscr(i)=vin(keep(i))
   10 continue
      call copy(vscr,vout,nkeep) 
      return
      end
