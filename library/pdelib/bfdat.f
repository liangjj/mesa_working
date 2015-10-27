*deck bfdat.f
c***begin prologue     bfdat
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
c***end prologue       bfdat
      subroutine bfdat(frmdsk,card,qtyp,qdtyp,left,right,
     1                 dleft,dright,lftbc,rtbc,npt,nmax,nfix,ngrid)
      implicit integer (a-z)
      character*(*) qtyp, qdtyp, card
      character*2 itoc
      character*80 chrkey, title
      logical logkey, fix, nfix, frmdsk
      real*8 fpkey, left, right, dleft, dright
      dimension nmax(ngrid), npt(ngrid)
      dimension nfix(2)
      common/io/inp, iout
      if(.not.frmdsk) then
         write(iout,1)
         qtyp=chrkey(card,'coordinate-type','x',' ')
         qdtyp=chrkey(card,'quadrature-type','legendre',' ')
         left=fpkey(card,'left-end-point',0.d0,' ')
         right=fpkey(card,'right-end-point',1.d0,' ')
         dleft=fpkey(card,'left-derivative',0.d0,' ')
         dright=fpkey(card,'right-derivative',0.d0,' ')
         lftbc=intkey(card,'left-boundary-condition',0,' ')
         rtbc=intkey(card,'right-boundary-condition',0,' ')
         fix=logkey(card,'fix-end-points',.false.,' ')
         nfix(1)=.false.
         nfix(2)=.false.
         if(fix) then
            nfix(1)=logkey(card,'fix-left-end-point',.false.,' ')
            nfix(2)=logkey(card,'fix-right-end-point',.false.,' ')
         endif
         do 10 i=1,ngrid
            npt(i)=intkey(card,'number-of-points-grid-'//itoc(i),0,' ')
            nmax(i)=npt(i)
            if(lftbc.eq.0) then
               nmax(i)=nmax(i)-1
            endif
            if(rtbc.eq.0) then
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
     1                 //itoc(i)//'" from rwf',1,left,0,' ')
            call iosys('read real "right end point grid = '
     1                 //itoc(i)//'" from rwf',1,right,0,' ')
            call iosys('read real "left derivative grid = '//
     1                  itoc(i)//'" from rwf',1,dleft,0,' ')
            call iosys('read real "right derivative grid = '//
     1                  itoc(i)//'" from rwf',1,dright,0,' ')
            call iosys('read integer "left boundary condition '//
     1                 'grid = '//itoc(i)//'" from rwf',1,lftbc,0,' ')
            call iosys('read integer "right boundary condition '//
     1                 'grid = '//itoc(i)//'" from rwf',1,rtbc,0,' ')
 20      continue   
      endif
      title='information for this coordinate'
      write(iout,3) title, left, right, dleft, dright, lftbc, rtbc
      do 30 i=1,ngrid
         title='data for grid '//itoc(i)
         write(iout,4) title, npt(i), nmax(i)       
 30   continue   
      return
 1    format(/,5x,'grid data from input')
 2    format(/,5x,'grid data from disk')
 3    format(/,15x,a60,/,5x,/,5x,     
     1             'left boundary                  = ',e15.8,/,5x,
     2             'right boundary                 = ',e15.8,/,5x,
     3             'derivative at left boundary    = ',e15.8,/,5x,
     4             'derivative at right boundary   = ',e15.8,/,5x,
     5             'left-boundary-condition        = ',i2,/,5x,
     6             'right-boundary-condition       = ',i2)
 4    format(/,15x,a60,/,5x,/,5x,
     1             'number of points          = ',i4,/,5x,
     2             'number of basis functions = ',i4)
      end






