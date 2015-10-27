c $Header: satshl.f,v 1.1 92/12/31 15:23:26 bis Exp $
*deck @(#)satshl.f	1.1 9/9/91
      subroutine satshl(icen)
      parameter (numshl=50)
      implicit real*8 (a-h,o-z)
      common/spheri/ nshell, nr(numshl), nthell(numshl), nphell(numshl),
     1               ntheta(numshl,numshl), nphi(numshl,numshl)
      common/spherr/ r(numshl),theta(numshl,numshl),phi(numshl,numshl)
      common /quadpt/ nsub(numshl), nth(numshl), nph(numshl)
      common/quadtp/ itype, jtype, ktype
      common /io/ inp, iout
      character*20 cpass
      character*1600 card
      character*20 itype, jtype, ktype
      character*3 itoc, chrv
      chrv=itoc(icen)
      len=length(chrv)
      write(iout,900) icen
      call posinp('$quadrature-atom-'//chrv(1:len),cpass)
      call cardin(card)
      nshell=intkey(card,'no-radial-shells',1,' ')
      write(iout,100) nshell
      call iosys ('write integer "number of shells for atom-'//
     1             chrv//'" to lamdat',1,nshell,0,' ')
      nedges=nshell+1
      call fparr(card,'radial-boundaries',r,nedges,' ')
      write(iout,200) (r(i),i=1,nedges)
      call intarr(card,'radial-quadrature-orders',nr,nshell,' ')
      write(iout,300) (nr(j),j=1,nshell)
      call intarr(card,'no-theta-divisions-per-shell',nthell,
     1            nshell,' ')
      call intarr(card,'no-phi-divisions-per-shell',nphell,
     1            nshell,' ')
      call iosys ('write integer "number radial points per shell for '//
     1            'atom-'//chrv//'" to lamdat',nshell,nr,0,' ')
      igrid=0
      do 10 i=1,nshell
         nthp1=nthell(i)+1
         nphp1=nphell(i)+1
         call posinp('$atom-'//chrv(1:len)//'-shell-'//itoc(i),cpass)
         call cardin(card)
         call fparr(card,'theta-boundaries',theta(1,i),nthp1,' ')
         call fparr(card,'phi-boundaries',phi(1,i),nphp1,' ')
         call intarr(card,'theta-quadrature-orders',ntheta(1,i),
     1               nthell(i),' ')
         call intarr(card,'phi-quadrature-orders',nphi(1,i),
     1               nphell(i),' ')
         nth(i)=0
         do 20 j=1,nthell(i)
            nth(i)=nth(i)+ntheta(j,i)
   20    continue
         nph(i)=0 
         do 30 j=1,nphell(i)
            nph(i)=nph(i)+nphi(j,i)
   30    continue
         nsub(i)=0
         do 40 j=1,nthell(i)
            do 50 k=1,nphell(i)
               nsub(i)=nsub(i)+nr(i)*ntheta(j,i)*nphi(k,i) 
               igrid=igrid+nr(i)*ntheta(j,i)*nphi(k,i)
   50       continue
   40    continue
         write(iout,400) i, (theta(j,i),j=1,nthp1)
         write(iout,500) (phi(j,i),j=1,nphp1)
         write(iout,600) (ntheta(j,i),j=1,nthell(i))
         write(iout,700) (nphi(j,i),j=1,nphell(i))
   10 continue
      call iosys ('write integer "number theta points per shell for '//
     1            'atom-'//chrv//'" to lamdat',nshell,nth,0,' ')
      call iosys ('write integer "number phi points per shell for '//
     1            'atom-'//chrv//'" to lamdat',nshell,nph,0,' ')
      call iosys ('write integer "total number of points per shell '//
     1            'for atom-'//chrv//'" to lamdat',nshell,nsub,
     2             0,' ')
      call iosys ('write integer "total number of points for atom-'//
     1             chrv//'" to lamdat',1,igrid,0,' ')
      write(iout,800) igrid
  100 format(2x,i4,' radial shells'/)
  200 format(' radial boundaries of shells',/,(2x,6f12.5))
  300 format(' r-quadrature orders:',30i3)
  400 format(//' for shell',i3/' theta edges:',8f12.4/(13x,8f12.4))
  500 format(/' phi edges:  ',8f12.4/(13x,8f12.4))
  600 format(' theta orders:',30i3)
  700 format(' phi orders:  ',30i3)
  800 format(/,5x,'no. quadrature points',1x,i6)
  900 format(/,5x,'quadrature data for atom = ',i2)
      return
      end







