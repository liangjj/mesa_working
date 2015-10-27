*deck gdat.f
c***begin prologue     gdat
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            grid input.
c***                   
c***references         
c
c***routines called    
c***end prologue       gdat
      subroutine gdat(card,qtyp,qdtyp,xl,xr,psil,psir,
     1                dpsil,dpsir,npt,nmax,nfix,ngrid)
      implicit integer (a-z)
      character*(*) qtyp, qdtyp, card
      character*2 itoc
      character*80 chrkey, title
      logical logkey, fix, nfix
      real*8 fpkey, xl, xr
      dimension nfix(2)
      common/io/inp, iout
      write(iout,1) ngrid
      qtyp=chrkey(card,'coordinate-type','x',' ')
      qdtyp=chrkey(card,'quadrature-type','legendre',' ')
      xl=fpkey(card,'xl',0.d0,' ')
      xr=fpkey(card,'xr',1.d0,' ')
      psil=intkey(card,'f-xl',0,' ')
      psir=intkey(card,'f-xr',0,' ')
      dpsil=intkey(card,'df-xl',1,' ')
      dpsir=intkey(card,'df-xr',1,' ')
      fix=logkey(card,'fix-end-points',.false.,' ')
      nfix(1)=.false.
      nfix(2)=.false.
      if(fix) then
         nfix(1)=logkey(card,'fix-left-end-point',.false.,' ')
         nfix(2)=logkey(card,'fix-right-end-point',.false.,' ')
      endif
      if(qdtyp.eq.'chebyshev-1'.or.qdtyp.eq.'chebyshev-2') then
         nfix(1)=.false.
         nfix(2)=.false.
      endif
      npt=intkey(card,'number-of-points',0,' ')
      nmax=npt
      if(psil.eq.0.and.nfix(1)) then
         nmax=nmax-1
      endif
      if(psir.eq.0.and.nfix(2)) then
         nmax=nmax-1
      endif
      write(iout,2) qdtyp, xl, xr, psil, psir, dpsil, dpsir, npt, nmax
      return
 1    format(/,5x,'data grid = ',i2)
 2    format(/,5x,
     1             'quadrature type                = ',a32,/,5x,     
     2             'left boundary                  = ',e15.8,/,5x,
     3             'right boundary                 = ',e15.8,/,5x,
     4             'left-boundary-condition        = ',i2,/,5x,
     5             'right-boundary-condition       = ',i2,/,5x,
     6             'derivative at left boundary    = ',i2,/,5x,
     7             'derivative at right boundary   = ',i2,/,5x,
     8             'number of points               = ',i4,/,5x,
     9             'number of basis functions      = ',i4)
      end






