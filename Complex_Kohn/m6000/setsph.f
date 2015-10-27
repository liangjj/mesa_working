*deck @(#)setsph.f	1.1 9/9/91
      subroutine setsph(ops,igrid,isave)
      parameter (numshl=50)
      implicit real*8 (a-h,o-z)
      logical iplot
      common /centsi/ncent,iplot
      common/spheri/ nshell,nr(numshl),nthell(numshl), nphell(numshl),
     1               ntheta(numshl,numshl), nphi(numshl,numshl)
      common/spherr/ r(numshl),theta(numshl,numshl),phi(numshl,numshl)
      common/quadtp/ itype, jtype, ktype
      common /io/ inp, iout
      character *10 cpass
      character *1600 card
      character *(*) ops
      character *20 itype, jtype, ktype, chrkey
      character *3 itoc
      itype=chrkey(ops,'type-radial-quadrature','legendre',' ')
      jtype=chrkey(ops,'type-theta-quadrature','legendre',' ')
      ktype=chrkey(ops,'type-phi-quadrature','simpson',' ')
      call posinp('$shells',cpass)
      call cardin(card)
      nshell=intkey(card,'no-radial-shells',1,' ')
      write(iout,431) nshell
431   format(2x,i4,' radial shells'/)
      nedges=nshell+1
      call fparr(card,'radial-boundaries',r,nedges,' ')
      write(iout,432) (r(i),i=1,nedges)
  432 format(' radial boundaries of shells',/,(2x,6f12.5))
      call intarr(card,'radial-quadrature-orders',nr,nshell,' ')
      write(iout,104)(nr(j),j=1,nshell)
104   format(' r-quadrature orders:',30i3)
      call intarr(card,'no-theta-divisions-per-shell',nthell,
     1            nshell,' ')
      call intarr(card,'no-phi-divisions-per-shell',nphell,
     1            nshell,' ')
      igrid=0
      do 10 i=1,nshell
      nthp1=nthell(i)+1
      nphp1=nphell(i)+1
      call posinp('$shell-'//itoc(i),cpass)
      call cardin(card)
      call fparr(card,'theta-boundaries',theta(1,i),nthp1,' ')
      call fparr(card,'phi-boundaries',phi(1,i),nphp1,' ')
      call intarr(card,'theta-quadrature-orders',ntheta(1,i),nthell(i),
     1            ' ')
      call intarr(card,'phi-quadrature-orders',nphi(1,i),nphell(i),' ')
      do 500 j=1,nthell(i)
      do 600 k=1,nphell(i)
      igrid=igrid+nr(i)*ntheta(j,i)*nphi(k,i)
  600 continue
  500 continue
      write(iout,100)i,(theta(j,i),j=1,nthp1)
      write(iout,101)(phi(j,i),j=1,nphp1)
      write(iout,102)(ntheta(j,i),j=1,nthell(i))
      write(iout,103)(nphi(j,i),j=1,nphell(i))
10    continue
      write(iout,700) igrid
      isave=0
      if (iplot) then
      isave=igrid/4
      endif
100   format(//' for shell',i3/' theta edges:',8f12.4/(13x,8f12.4))
101   format(/' phi edges:  ',8f12.4/(13x,8f12.4))
102   format(' theta orders:',30i3)
103   format(' phi orders:  ',30i3)
700   format(/,5x,'no. quadrature points',1x,i6)
      return
      end
