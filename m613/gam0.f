*deck @(#)gam0.f	5.1 11/6/94
c***begin prologue     gam0.f
c***date written       880528   (yymmdd)
c***revision date      @(#)gam0.f	5.1
c***keywords           gam0, gaussian integrals
c***author             schneider, barry (lanl)
c***source             @(#)gam0.f	5.1 11/6/94
c***purpose            starting value for recursion for gaussian integrals
c***description        computes sqrt(pi)/2. * erf(t) to begin recursion 
c***                   for even values of n in basint
c***                   
c
c***references         bis notes
c
c***routines called
c***end prologue       gam0.f
      function gam0(t)
      implicit integer (a-z)
      common/io/inp,iout
c      character *5 card
      real *8 t, gam0, erf, sqpi2
      common/card/card
      data sqpi2/.88622692545275801364d+00/
c      if (card.eq.'print') then
c      write (iout,1)
c    1 format(/,5x,'gam0 called')
c      endif
      gam0=sqpi2*erf(t)
      return
      end
