c $Header$
*deck vhalf.f
c***begin prologue     vhalf
c***date written       921223   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           vhalf, link 6201, wavefunction
c***author             schneider, barry (lanl)
c***source             m6203
c***purpose            scale diagonal elements of local potential
c***                   by .5       
c***references         
c***routines called
c***end prologue       vhalf
      subroutine vhalf (vloc,ns,nr,nth,nph,nchan)
      implicit integer(a-z)
      real*8 vloc
      dimension vloc(*), nr(ns), nth(ns), nph(ns)
      ntri=nchan*(nchan+1)/2
      count=0
      do 10 is=1,ns
         npts=nr(is)*nth(is)*nph(is)
         do 20 ic=1,nchan
            idiag=ic*(ic+1)/2
            loc=count+(idiag-1)*npts
            call sscal(npts,.5d0,vloc(loc+1),1)
   20    continue
         count=count+npts*ntri
   10 continue
      return
      end
