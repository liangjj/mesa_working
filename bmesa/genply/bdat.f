*deck bdat.f
c***begin prologue     bdat
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            basis function input.
c***                   
c***references         
c
c***routines called    
c***end prologue       bdat
      subroutine bdat(frmdsk,card,qtyp,qdtyp,xl,xr,psil,psir,
     1                dpsil,dpsir,npt,nmax,nfix,ngrid)
      implicit integer (a-z)
      character*(*) qtyp, qdtyp, card
      character*2 itoc
      character*80 chrkey, title
      logical logkey, fix, nfix, frmdsk
      real*8 fpkey, xl, xr
      dimension nmax(ngrid), npt(ngrid)
      dimension nfix(2)
      common/io/inp, iout
      if(.not.frmdsk) then
         write(iout,1)
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
         do 10 i=1,ngrid
            npt(i)=intkey(card,'number-of-points-grid-'//itoc(i),0,' ')
            nmax(i)=npt(i)
            if(psil.eq.0.and.nfix(1)) then
               nmax(i)=nmax(i)-1
            endif
            if(psir.eq.0.and.nfix(2)) then
               nmax(i)=nmax(i)-1
            endif
 10      continue   
      else
         write(iout,2)
         do 20 i=1,ngrid
            call iosys ('read character "coordinate type '//
     1                  'grid = '//itoc(i)//'" from rwf',
     2                   -1,0,0,qtyp)
            call iosys('read integer "number of points '//
     1                 'grid = '//itoc(i)//'" from rwf',1,
     2                  npt(i),0,' ') 
            call iosys('read integer "number of functions '//
     1                 'grid = '//itoc(i)//'" from rwf',1,
     2                  nmax(i),0,' ')                          
            call iosys('read real "left end point grid = '
     1                 //itoc(i)//'" from rwf',1,xl,0,' ')
            call iosys('read real "right end point grid = '
     1                 //itoc(i)//'" from rwf',1,xr,0,' ')
            call iosys('read integer "left derivative grid = '//
     1                  itoc(i)//'" from rwf',1,dpsil,0,' ')
            call iosys('read integer "right derivative grid = '//
     1                  itoc(i)//'" from rwf',1,dpsir,0,' ')
            call iosys('read integer "left boundary condition '//
     1                 'grid = '//itoc(i)//'" from rwf',1,psil,0,' ')
            call iosys('read integer "right boundary condition '//
     1                 'grid = '//itoc(i)//'" from rwf',1,psir,0,' ')
 20      continue   
      endif
      title='information for this coordinate'
      write(iout,3) title, qdtyp, xl, xr, psil, psir, dpsil, dpsir
      do 30 i=1,ngrid
         title='data for grid '//itoc(i)
         write(iout,4) title, npt(i), nmax(i)       
 30   continue   
      return
 1    format(/,5x,'grid data from input')
 2    format(/,5x,'grid data from disk')
 3    format(/,15x,a60,/,5x,/,5x,
     1             'quadrature type                = ',a32,/,5x,     
     2             'left boundary                  = ',e15.8,/,5x,
     3             'right boundary                 = ',e15.8,/,5x,
     4             'left-boundary-condition        = ',i2,/,5x,
     5             'right-boundary-condition       = ',i2,/,5x,
     6             'derivative at left boundary    = ',i2,/,5x,
     7             'derivative at right boundary   = ',i2)
 4    format(/,15x,a60,/,5x,/,5x,
     1             'number of points          = ',i4,/,5x,
     2             'number of basis functions = ',i4)
      end






