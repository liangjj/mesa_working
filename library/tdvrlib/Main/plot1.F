*deck plot1.f
c***begin prologue     plot1
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           plot, time-dependent, wavefunction
c***author             schneider, barry (nsf)
c***source             
c***purpose            plot wavefunction.
c***                   
c***references         
c
c***routines called    
c***end prologue       plot1
      subroutine plot1(rex,rapprox,iex,iapprox,t,nt,pun)
      implicit integer (a-z)
      real*8 rex, rapprox, iex, iapprox, t
      dimension rex(nt), rapprox(nt), iex(nt), iapprox(nt)
      dimension t(nt)
      common/io/inp, iout	 
      do 10 i=1,nt,5
         write(pun,1) t(i), rex(i), rapprox(i), iex(i), iapprox(i)
 10   continue   
      return
 1    format(e15.8,1x,'&',e15.8,1x,'&',e15.8,'&',e15.8,1x,
     1                '&',e15.8,3x,'\\ \hline')
      end       
