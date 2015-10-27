*deck @(#)filff.f	1.1 9/8/91
c***begin prologue     filff
c***date written       930417   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           filff, link 6018
c***author             schneider, barry (nsf)  
c***source             m6018
c***purpose            fill diagonal free-free hamiltonian matrices
c***                   
c***references        
c
c***routines called    iosys, util and mdutil
c***end prologue       filff
      subroutine filff (ovpp,ovpm,hpvhp,opp,opm,omm,hpp,maxlm,nl)
      implicit integer(a-z)
      complex*16 ovpp, ovpm, hpvhp, opp, opm, hpp
      real*8 omm
      dimension ovpp(maxlm), ovpm(maxlm), hpvhp(maxlm)
      dimension opp(nl), opm(nl), omm(nl), hpp(nl)
      common /io/ inp,iout
      do 10 l=1,nl
         opp(l)=ovpp(l)
         opm(l)=ovpm(l)
         omm(l)=imag(opm(l))
         hpp(l)=hpvhp(l)
   10 continue
      return
      end
