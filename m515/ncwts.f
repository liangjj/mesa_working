*deck @(#)ncwts.f	5.1  11/28/95
      subroutine ncwts(pt,wt,n,stp,scale)
c***begin prologue     ncwts.f
c***date written       920324   (yymmdd)
c***revision date      11/6/94      (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             @(#)ncwts.f	5.1   11/28/95
c***purpose            newton cotes weights
c***
c***description
c   these are the newton-cotes weights for doing an integral from a0 to an
c      as (a0,a1) +(a1,a2) +(a2,a3) +... 
c      note that rthree, e.g., implements second-order newton-cotes, etc.
c
c               
c               
c***references
c
c***routines called
c
c***end prologue       ncwts.f
c
      implicit integer (a-z)
c
      integer n
      real*8 wt(n,n-1), pt(n)
      real*8 stp
      logical scale
c
      integer inp,iout
      common /io/ inp, iout
c
c     --- call appropriate subroutine, depending on order.
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
            write(iout,*) 'requested newton-cotes order:',n
           call lnkerr('quadrature not available')
      endif
      if (scale) then
          do 20 i=1,n-1
             do 30 j=1,n
                wt(j,i)=wt(j,i)*pt(j)*pt(j)
   30        continue
   20     continue           
      endif
c
c
      return
      end
