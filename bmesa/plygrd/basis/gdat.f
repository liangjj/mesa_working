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
     1                dpsil,dpsir,nsubg,npt,nfix,maxd)
      implicit integer (a-z)
      character*(*) qtyp, qdtyp, card
      character*80 chrkey
      logical logkey, fix, nfix
      real*8 fpkey, xl, xr
      dimension nfix(2), npt(*)
      common/io/inp, iout
      nsubg=intkey(card,'number-of-subgrids',1,' ')
      write(iout,1) qtyp
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
      call intarr(card,'number-of-points-per-subgrid',npt,nsubg,' ')
      write(iout,2) qdtyp, xl, xr, psil, psir, dpsil, dpsir, nsubg
      write(iout,3) (npt(i), i=1,nsubg)
      do 10 i=1,nsubg
         maxd=max(maxd,npt(i))
 10   continue   
      return
 1    format(/,5x,'data coordinate type = ',a24)
 2    format(/,5x,
     1             'quadrature type                = ',a32,/,5x,     
     2             'left boundary                  = ',e15.8,/,5x,
     3             'right boundary                 = ',e15.8,/,5x,
     4             'left-boundary-condition        = ',i2,/,5x,
     5             'right-boundary-condition       = ',i2,/,5x,
     6             'derivative at left boundary    = ',i2,/,5x,
     7             'derivative at right boundary   = ',i2,/,5x,
     8             'number of sub-grids            = ',i4)
 3    format(/,5x,'number of points in each subgrid = ',10(i3,1x))
      end






