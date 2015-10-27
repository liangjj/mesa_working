*deck filpot
c***begin prologue     filpot
c***date written       940213   (yymmdd)
c***revision date               (yymmdd)
c***keywords           potential
c***author             schneider, barry (nsf)
c***source             
c***purpose            potential on grid
c***description        
c***references       
c
c***routines called
c***end prologue       filpot
      subroutine filpot(vij,vcij,pt,range,chni,chnj,type)
      implicit integer (a-z)
      common /io/ inp, iout
      real*8 vij, vcij, pt, range, d
      character*(*) type
      if (type.eq.'square-well') then
          if (pt.le.range) then
              vij=vcij
          else
              vij=0.d0
          endif        
      elseif (type.eq.'exponential') then
              vij=vcij*exp(-pt)
      elseif (type.eq.'hazi-taylor') then
              vij=vcij*pt*pt*exp(-pt)
      elseif (type.eq.'henry') then
           if (chni.eq.chnj) then 
               d=2.d0*vcij
               d=sqrt(d)
               vij=-vcij/(pt*pt+d)**2
           else
               vij=-vcij*(1.d0-exp(-pt)**3/(pt*pt))
           endif                                          
      else
          call lnkerr('error in potential type')
      endif
      return                                               
      end
