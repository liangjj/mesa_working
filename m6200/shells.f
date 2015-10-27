*deck shells.f
c***begin prologue     shells
c***date written       9308021   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           shells
c***author             schneider, barry (nsf)
c***source             m6200
c***purpose            calculate angular and radial quadrature points
c***                   and weights
c***references         
c
c***routines called    gaussq ( math )
c***end prologue       shells
      subroutine shells(type,icen,cen,r,nr,eta,rpt,thpt,phpt,wtr,wtsum,
     1                  wtth,wtph,angleb,work,wtleb,sthet,sphi,cphi,
     2                  yuk,scr,grid,wt,spgrid,ns,nrmax,nthet,nphi,
     3                  ngrid,nwts,ncen,ncplus,prnt,nleb,nang,nonsep,
     4                  yukon,nodisk)
      implicit real*8 (a-h,o-z) 
      integer thefac, phifac
      logical prnt, nonsep, yukon, nodisk
      character*30 rtyp, thtyp, phtyp, str
      character*3 itoc, chra
      character*(*) type
      dimension rpt(nrmax,ns), thpt(nthet), phpt(nphi)
      dimension wtr(*), wtth(nthet), wtph(nphi), wt(nwts), wtleb(*)
      dimension wtsum(*), nr(ns), r(ns+1), delthe(2), delphi(2)
      dimension cen(3,ncplus), sthet(nthet), sphi(nphi), cphi(nphi)
      dimension grid(3,ngrid), scr(*), dummy(2), yuk(*), eta(ncplus)
      dimension spgrid(3,ngrid), angleb(*), work(*)
      common/io/inp,iout
      data pi / 3.141592653589793d0  /
      if (type.eq.'atom') then
          write(iout,*)
          chra=itoc(icen)
          len=length(chra)
          write(iout,*) 'generating grid for atom-'//chra
          str='atom-'//chra
      else
          str='scattering center'
          write(iout,*)
          write(iout,*) 'generating grid for scattering center'
      endif                      
      call iosys('read character "radial quadrature type '//
     1            str//'" from lamdat',-1,0,0,rtyp)
      call iosys('read real "theta range '//
     1            str//'" from lamdat',2,delthe,0,' ') 
      call iosys('read real "phi range '//
     1            str//'" from lamdat',2,delphi,0,' ') 
      nedges=ns+1
      call iosys ('read real "radial edges '//str//
     1            '" from lamdat',nedges,r,0,' ')
      call iosys ('create real "radial points '//str//
     1            '" on lamdat',nrmax*ns,0,0,' ')
      call iosys ('create real "scaled radial weights '//
     1             str//'" on lamdat',nrmax*nrmax*ns,
     2             0,0,' ')
      call iosys ('create real "unscaled radial weights '//
     1             str//'" on lamdat',nrmax*nrmax*ns,
     2             0,0,' ')
      if (rtyp.eq.'newton-cotes') then
          call iosys ('create real "summed radial weights '//
     1                 str//'" on lamdat',nrmax*ns,0,0,' ')
      endif
      if (.not.nodisk) then     
           call iosys ('create real "yukawa potential '//
     1                  str//' coordinates" on lamdat',ngrid,
     2                  0,0,' ')
           call iosys ('create real "atomic grid '//str//
     1                 '" on lamdat',3*ngrid,0,0,' ')
           call iosys ('create real "spherical atomic grid '//str//
     1                 '" on lamdat',3*ngrid,0,0,' ')
           call iosys ('create real "unscaled atomic weights '//
     1                  str//'" on lamdat',nwts,0,0,' ')
       endif     
      if (nonsep) then 
c      
c          in a lebedev quadrature the number of theta and phi points are
c          identical since the quadrature is non-separable. the combined
c          angular weight is stored in wtth.      
c
          call lebdev(angleb,angleb,thpt,sthet,phpt,sphi,cphi,wtleb,
     1                nleb,nang,str) 
c
      else
          call iosys('read character "theta quadrature type '//
     1                str//'" from lamdat',-1,0,0,thtyp)
          call iosys('read character "phi quadrature type '//
     1                str//'" from lamdat',-1,0,0,phtyp)         
          difth=delthe(2)-delthe(1)
          difph=delphi(2)-delphi(1)
c     multiplicative factors for theta and phi integration if symmetry
c     allows integrals to be calculated on subdomains.
          thefac=pi/difth
          phifac=2.d0*pi/difph
          write(iout,1) thefac, phifac
          tmp=cos(delthe(1))
          ampt=tmp
          abpt=cos(delthe(2))
          ampt=(ampt-abpt)*.5d0
          abpt=(abpt+tmp)*.5d0
          call gaussq(thtyp,nthet,0.d+00,0.d+00,0,dummy,scr,
     1                thpt,wtth)
          do 10 i=1,nthet
             thpt(i)=ampt*thpt(i)+abpt
             sthet(i)=sqrt((1.d0-thpt(i)*thpt(i)))
             wtth(i)=thefac*ampt*wtth(i)
   10     continue
          ampp=(delphi(2)-delphi(1))*.5d0
          abpp=(delphi(2)+delphi(1))*.5d0
          call gaussq(phtyp,nphi,0.d+00,0.d+00,0,dummy,scr,
     1                phpt,wtph)
          do 20 i=1,nphi
             phpt(i)=ampp*phpt(i)+abpp
             sphi(i)=sin(phpt(i))
             cphi(i)=cos(phpt(i))
             wtph(i)=phifac*ampp*wtph(i)
   20     continue
          call iosys ('write real "theta points '//str//
     1                '" to lamdat',nthet,thpt,0,' ')
          call iosys ('write real "theta weights '//str//
     1                '" to lamdat',nthet,wtth,0,' ')
          call iosys ('write real "phi points '//str//
     1                '" to lamdat',nphi,phpt,0,' ')
          call iosys ('write real "phi weights '//str//
     1                '" to lamdat',nphi,wtph,0,' ')
          nang=nthet*nphi   
      endif
      if (prnt) then
          write(iout,*)
          write(iout,*) 'cos(theta) points and weights'
          write(iout,2) (thpt(ii),ii=1,nthet)
          write(iout,2) (wtth(ii),ii=1,nthet)
          write(iout,*)
          write(iout,*) 'phi points and weights'
          write(iout,2) (phpt(ii),ii=1,nphi)
          write(iout,2) (wtph(ii),ii=1,nphi)
      endif
c   
c         do the radial quadrature points. they are divided into shells.
c
c
      ngcnt=0
      nwtcnt=0
      aprvol=0.d0
      yukint=0.d0          
      call rzero(yuk,ngrid)
      icount=0
      do 30 i=1,ns
         if (r(i).eq.0.d0) then
             r(i)=1.d-10
         endif    
         locwtr=nrmax*(i-1)+1
c     get radial quadrature points and weights for this shell
         if (rtyp.eq.'legendre') then
             ampr=(r(i+1)-r(i))*.5d0
             abpr=(r(i+1)+r(i))*.5d0
             call gaussq(rtyp,nr(i),0.d+00,0.d+00,0,dummy,scr,
     1                   rpt(1,i),wtr(locwtr))
             nwtsr=nr(i)
             icount=icount+nr(i)
             loccnt=locwtr
              do 40 j=1,nr(i)
                 numr=numr+1
                 rpt(j,i)=ampr*rpt(j,i)+abpr
                 wtr(loccnt)=ampr*wtr(loccnt)
                 loccnt=loccnt+1
   40         continue
         elseif(rtyp.eq.'newton-cotes') then
                locwtr=nrmax*(nrmax-1)*(i-1)+1
                call necote(r(i),r(i+1),rpt(1,i),wtr(locwtr),nr(i),
     1                      .false.)
                nwtsr=(nr(i)-1)*nr(i)
                icount=icount+nr(i)
         else
                call lnkerr('error in radial quadrature')
         endif
         call iosys ('write real "radial points '//str//
     1               '" to lamdat without rewinding',nr(i),
     2                  rpt(1,i),0,' ')
         call iosys ('write real "unscaled radial weights '//
     1                str//'" to lamdat without '//
     2                'rewinding',nwtsr,wtr(locwtr),0,' ')
         if(rtyp.eq.'newton-cotes') then
            call sumncw(wtr(locwtr),wtsum,nr(i))
            call iosys ('write real "summed radial weights '//
     1                   str//'" to lamdat without rewinding',
     2                   nr(i),wtsum,0,' ')
         endif
c      the procedure is different if the angular quadrature is separable or
c      non-separable in  (theta,phi)
         call mkgr(grid(1,ngcnt+1),spgrid(1,ngcnt+1),rpt(1,i),thpt,
     1             sthet,sphi,cphi,cen(1,icen),nr(i),nthet,
     2             nphi,nang,nonsep)
         if (yukon) then
             call yukawa(yuk(ngcnt+1),grid(1,ngcnt+1),eta,cen,
     1                   nr(i),nthet,nphi,ncen,nang,nonsep)
         endif 
         call mkwt(wtr(locwtr),wtr(locwtr),wtth,wtph,wtleb,
     1             wt(nwtcnt+1),work,nr(i),nthet,nphi,rtyp,nang,nonsep)
         if (yukon) then
             call mkyunt(rpt(1,i),wt(nwtcnt+1),yuk(ngcnt+1),aprvol,
     1                   yukint,work,nr(i),nthet,nphi,rtyp,nang,
     2                   nonsep)
         endif
         ngcnt=ngcnt+nr(i)*nang
         nwtcnt=nwtcnt+nwtsr*nang
         call scalwt(rpt(1,i),wtr(locwtr),wtr(locwtr),work,nr(i),rtyp)
         call iosys ('write real "scaled radial weights '//
     1                str//'" to lamdat without '//
     2                'rewinding',nwtsr,wtr(locwtr),0,' ')           
         if (prnt) then
             write(iout,*) '     data on radial points for shell = ',i
             write(iout,*)
             write(iout,*) 'radial points and weights'
             write(iout,2) (rpt(ii,i),ii=1,nr(i))
             write(iout,2) (wtr(ii),ii=locwtr,locwtr+nwtsr-1)
         endif
   30 continue
      if (ngcnt.ne.ngrid) then
          call lnkerr('error in grid point count')
      endif
      if (nwtcnt.ne.nwts) then
          call lnkerr('error in grid weight count')
      endif
      call iosys ('write integer "total number of radial '//
     1            'points '//str//'" to lamdat',1,icount,0,' ')
      if (.not.nodisk) then
           call iosys ('write real "yukawa potential '//
     1                 str//' coordinates" to lamdat',ngcnt,
     2                 yuk,0,' ')
           call iosys ('write real "atomic grid '//
     1                  str//'" to lamdat',3*ngcnt,grid,0,' ')
           call iosys ('write real "spherical atomic grid '//
     1                  str//'" to lamdat',3*ngcnt,spgrid,0,' ')
           call iosys ('write real "unscaled atomic weights '//
     1                  str//'" to lamdat',nwts,wt,0,' ')
      endif
      if (yukon) then
          yukext=0.d0
          do 80 nc=1,ncen
             yukext=yukext+4.d+00*pi/eta(nc)**2
   80     continue
          exvol=4.d0*pi*r(ns+1)**3/3.d0
          write(iout,3) exvol
          write(iout,4) aprvol       
          write(iout,5) yukext
          write(iout,6) yukint
      endif
      return
    1 format(' multiplicative factors for theta and phi  :',2f10.5)
    2 format(/4(2x,e15.8))
    3 format(/,5x,'exact       value of volume =',e15.8)
    4 format(/,5x,'approximate value of volume =',e15.8)
    5 format(/,5x,'exact     value of yukawa integral = ',e15.8)
    6 format(/,5x,'numerical value of yukawa integral = ',e15.8)
      end
