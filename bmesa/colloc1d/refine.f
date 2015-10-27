*deck refine.f
c***begin prologue     refine
c***date written       910128   (yymmdd)
c***revision date               (yymmdd)
c***keywords           potential
c***author             schneider, barry(lanl)
c***source             @(#)m6020
c***purpose            refine grids
c***
c***
c***references
c
c***routines called    util
c***end prologue
      subroutine refine(type,quad,nin,nout,nwts)
      implicit integer (a-z)
      character*(*) type, quad
      common /io/ inp, iout
      data title /1/
      save title
      if(title.eq.1) then
         if(type.eq.'collocation') then
            write(iout,1)
         elseif(type.eq.'finite-difference') then
            write(iout,2)
         elseif(type.eq.'fitting') then
            write(iout,3)
         else
            call lnkerr('error in call to refine')
         endif
         title=0
      endif
      nout=nin+nin-1
      nwts=nout
      if(quad.eq.'newton-cotes'.or.quad.eq.'fixed') then
         nwts=nwts*(nwts-1)
      endif
      return
 1    format(/,6x,'   grid     collocation points')
 2    format(/,6x,'   grid     finite difference points')
 3    format(/,6x,'   grid     matching points')
      end







