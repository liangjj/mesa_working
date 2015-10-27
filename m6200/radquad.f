*deck radquad.f
c***begin prologue     radquad
c***date written       930802   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           radquad, link m6200
c***author             schneider, barry (nsf)
c***source             m6200
c***purpose            read in radial quadrature information.
c
c***routines called
c***end prologue       radquad
      subroutine radquad(card,str,nr,r,numshl,igrid,iwt,nrmax,maxgrd,
     1                   maxshl,ntrad,nang)
      implicit integer (a-z)
      character*(*) card, str
      logical shinfo
      character*30 rtyp, chrkey
      real*8 fpkey, r, rleft, rright, delta 
      dimension nr(numshl), r(numshl+1)
      common /io/ inp, iout
      len=length(str)      
      shinfo=logkey(card,'input-shell-information',.false.,' ')
      if (shinfo) then
          nshell=intkey(card,'no-radial-shells',1,' ')
          nedges=nshell+1
          call fparr(card,'radial-boundaries',r,nedges,' ')
          call intarr(card,'radial-quadrature-orders',nr,nshell,' ')
          if (r(1).eq.0.d0) then
              r(1)=1.d-10
          endif
      else
          rleft=1.d-10
          rright=fpkey(card,'maximum-radial-point',10.d0,' ')
          delta=fpkey(card,'interval-size',1.d0,' ')
          nshell=rright/delta
          nrdef=intkey(card,'default-quadrature-size',5,' ')
          do 10 ns=1,nshell
             nr(ns)=nrdef
   10     continue
          r(1)=rleft
          nedges=nshell+1
          do 20 ns=2,nedges
             r(ns)=r(ns-1)+delta
   20     continue
      endif
      maxshl=max(maxshl,nshell)
      write(iout,100) (nr(j),j=1,nshell)
      write(iout,200) (r(i),i=1,nedges)   
      write(iout,300) nshell
      rtyp=chrkey(card,'type-radial-quadrature','legendre',' ')
      igrid=0
      icrad=0
      iwt=0
      do 30 i=1,nshell
         igrid=igrid+nr(i)*nang
         nrmax=max(nrmax,nr(i))
         icrad=icrad+nr(i)
         iwt=iwt+nr(i)*(nr(i)-1)*nang
   30 continue
      if (rtyp.ne.'newton-cotes') then
          iwt=igrid
      endif    
      ntrad=max(icrad,ntrad)
      maxgrd=max(igrid,maxgrd)
      write(iout,400) igrid
      call iosys ('write integer "number of shells '//
     1             str//'" to lamdat',1,nshell,0,' ')
      call iosys ('write integer "number radial points per shell '//
     1            str//'" to lamdat',nshell,nr,0,' ')
      call iosys ('write real "radial edges '//str//'" to lamdat',
     1             nedges,r,0,' ')
      call iosys('write character "radial quadrature type '//
     1            str//'" to lamdat',0,0,0,rtyp)
      call iosys ('write integer "total number of radial points '//
     1             str//'" to lamdat',1,icrad,0,' ') 
      call iosys ('write integer "total number of points '//
     1             str//'" to lamdat',1,igrid,0,' ')                  
      return
  100 format(' r-quadrature orders:',(/20i3))
  200 format(' radial boundaries of shells',/,(2x,6f12.5))
  300 format(2x,i4,' radial shells'/)
  400 format('total number of grid points',1x,i5)
      end

