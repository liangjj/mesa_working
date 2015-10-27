*deck necote
      subroutine necote(rleft,rright,pt,wt,n,scale)
c***begin prologue     necote
c***date written       920324   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***purpose            newton cotes weights
c***
c***description
c***            
c               
c               
c***references
c
c***routines called
c
c***end prologue       necote
c
      implicit integer (a-z)
      real*8 rleft, rright, pt, wt, stp
      dimension wt(n,n-1), pt(n)
      logical scale
      common /io/ inp, iout
      stp=(rright-rleft)/(n-1)
      pt(1)=rleft
      do 10 i=2,n
         pt(i)=pt(i-1)+stp
   10 continue
      if ( n.eq.2 ) then
           call rtwo(wt,stp)
      elseif ( n.eq.3 ) then
           call rthree(wt,stp)
      elseif ( n.eq.4 ) then
           call rfour(wt,stp)
      elseif ( n.eq.5 ) then
           call rfive(wt,stp)
      elseif ( n.eq.6 ) then
           call rsix(wt,stp)
      elseif ( n.eq.7 ) then
           call rseven(wt,stp)
      elseif ( n.eq.8 ) then
           call reight(wt,stp)
      elseif ( n.eq.9 ) then
           call rnine(wt,stp)
      else
           call lnkerr('error in quadrature order')
      endif
      if (scale) then
          do 20 i=1,n-1
             do 30 j=1,n
                wt(j,i)=wt(j,i)*pt(j)*pt(j)
   30        continue
   20     continue           
      endif
      return
      end















