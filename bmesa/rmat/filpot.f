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
      subroutine filpot(vij,pot,pt,range,npts,chni,chnj,type)
      implicit integer (a-z)
      common /io/ inp, iout
      real*8 pot, pt, range, vij, d
      character*(*) type
      dimension pot(npts), pt(npts)
      if (type.eq.'square-well') then
          do 10 i=1,npts
             if (pt(i).le.range) then
                 pot(i)=vij
             else
                 pot(i)=0.d0
             endif        
   10     continue
      elseif (type.eq.'exponential') then
          do 20 i=1,npts
             pot(i)=vij*exp(-pt(i))
   20     continue
      elseif (type.eq.'hazi-taylor') then
          do 30 i=1,npts
             pot(i)=vij*pt(i)*pt(i)*exp(-pt(i))
   30     continue
      elseif (type.eq.'coulomb') then
          do 40 i=1,npts
             pot(i)=vij/pt(i)
   40     continue             
      elseif (type.eq.'henry') then
           if (chni.eq.chnj) then 
               d=2.d0*vij
               d=sqrt(d)
               do 50 i=1,npts
                  pot(i)=-vij/(pt(i)*pt(i)+d)**2
   50          continue
           else
               do 60 i=1,npts
                  pot(i)=-vij*(1.d0-exp(-pt(i)))**3/(pt(i)*pt(i))
   60          continue
           endif                                          
      else
          call lnkerr('error in potential type')
      endif
      return                                               
      end
