*deck drvint.f
c***begin prologue     drvint
c***date written       930922   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           drvint, m6201 y(l,m), projections
c***author             schneider, barry (nsf)
c***source             m6201
c***purpose            driver for integral equation
c***references         none
c
c***routines called
c***end prologue       drvint
      subroutine drvint (flm,psilm,rad,j,y,wt,scr,nrad,l,m,nshell,nr,
     1                   nints,type,prnt)
      implicit integer (a-z)
      real*8 psilm, flm, j, y, wt, scr, rad
      character*(*) type
      character*80 title
      character*3 itoc
      logical prnt
      dimension flm(nrad,*), psilm(nints,*), j(nrad,0:l), y(nrad,0:l)
      dimension wt(*), nr(nshell), scr(*), rad(nrad)
      common /io/ inp, iout
      locpf=1
c----------------------------------------------------------------------c
c         for each m, compute the radial integral                      c
c--------------------------------------------------------------------- c     
      locold=1
      do 10 mu=0,m
         no=2
         if ( mu.eq.0) then
              no=1
         endif
         do 20 count=1,no     
            if (type.eq.'legendre') then
                call fbgint(flm(1,locpf),psilm(1,locpf),j,y,wt,
     1                      scr,l,mu,nrad,prnt)
            elseif (type.eq.'newton-cotes') then
                call fbnint(flm(1,locpf),psilm(1,locpf),rad,j,y,wt,
     1                      scr,l,mu,nrad,nshell,nr,nints,prnt)
            else
                call lnkerr('quadrature type error')
            endif
            if (prnt) then
                lv=l-mu+1
                title='radial partial waves for m = '//itoc(mu)
                call prntfm(title,psilm(1,locpf),nints,lv,nints,lv,iout)
            endif                  
            locpf=locpf+nrad*(l-mu+1)
            locold=locold+l-mu+1
   20    continue
   10 continue
      return
      end
