*deck drvint.f
c***begin prologue     drvint
c***date written       921223   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           drvint, link 2702, y(l,m), projections
c***author             schneider, barry (nsf)
c***source             m6202
c***purpose            driver for integral equation
c***references         none
c
c***routines called
c***end prologue       drvint
      subroutine drvint (psi,flm,j,y,wt,old,nrad,mval,lndex,nl,nm,ltop,
     1                   ldim,nshell,type)
      implicit integer (a-z)
      real*8 psi, flm, j, y, wt, old
      character*(*) type
      dimension psi(*), flm(*), j(*), y(*), nrad(nshell)
      dimension mval(nm), nl(nm), lndex(*), wt(*), old(*)
      common /io/ inp, iout
      locpsi=1
      locflm=1
      locj=1
      locy=1
      locind=1
      locwt=1
      call rzero(old,ldim)
      do 10 ns=1,nshell
c----------------------------------------------------------------------c
c         for each m, compute the forward radial integral              c
c--------------------------------------------------------------------- c     
         do 20 mu=1,nm
            if (type.eq.'gauss') then
                call fgauss(flm(locflm),psi(locpsi),old,j(locj),
     1                      y(locy),wt(locwt),lndex(locind),nl(mu),
     2                      ltop,nrad(ns))
                locflm=locflm+nl(mu)*nrad(ns)
                locpsi=locpsi+nl(mu)*nrad(ns)
                locind=locind+nl(mu)
            elseif (type.eq.'newton-cotes') then
                call fnewt(flm(locflm),psi(locpsi),j(locj),y(locy),
     1                     wt(locwt),lndex(locind),nl(mu),
     2                     ltop,nrad(ns))
            else
                    call lnkerr('quadrature type error')
            endif
   20    continue              
         locwt=locwt+nrad(ns)
         locj=locj+nrad(ns)*(ltop+1)
         locy=locy+nrad(ns)*(ltop+1)
      return
      end
