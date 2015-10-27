*deck timpts.f
c***begin prologue     timpts
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            time-dependent potential
c***                   
c***description        automat time grid
c***                   
c***references         
c
c***routines called    
c***end prologue       timpts 
      subroutine timpts(t,npts,ntreg)
      implicit integer (a-z)
      real*8 t, delt
      logical posinp
      character*1600 card
      dimension t(ntreg+1), npts(ntreg)
      common/io/inp, iout
      if ( posinp('$pts',card) ) then
           read(inp,*) t(1), delt, npts(1)
      else
           call lnkerr('not found')
      endif
      do 10 i=2,ntreg+1       
         t(i)=t(i-1)+delt
 10   continue
      do 20 i=2,ntreg
         npts(i)=npts(i-1)
 20   continue   
      return
      end       


