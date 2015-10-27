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
      subroutine bfdat(card,qtyp,qdtyp,left,right,
     1                 dleft,dright,lftbc,rtbc,npt,nmax,nfix)
      implicit integer (a-z)
      character*(*) qtyp, qdtyp, card
      character*2 itoc
      character*80 chrkey, title
      logical logkey, fix, nfix
      real*8 fpkey, left, right, dleft, dright
      dimension nfix(2)
      common/io/inp, iout
      write(iout,1)
      qtyp=chrkey(card,'coordinate-type','x',' ')
      qdtyp=chrkey(card,'quadrature-type','legendre',' ')
      qdtyp='legendre'
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
      npt=intkey(card,'number-of-points',0,' ')
      nmax=npt
      if(lftbc.eq.0) then
         nmax=nmax-1
      endif
      if(rtbc.eq.0) then
         nmax=nmax-1
      endif
      title='information for this coordinate'
      write(iout,2) title, left, right, dleft, dright, lftbc, rtbc
      write(iout,3) title, npt, nmax       
      return
 1    format(/,5x,'grid data from input')
 2    format(/,15x,a60,/,5x,/,5x,     
     1             'left boundary                  = ',e15.8,/,5x,
     2             'right boundary                 = ',e15.8,/,5x,
     3             'derivative at left boundary    = ',e15.8,/,5x,
     4             'derivative at right boundary   = ',e15.8,/,5x,
     5             'left-boundary-condition        = ',i2,/,5x,
     6             'right-boundary-condition       = ',i2)
 3    format(/,15x,a60,/,5x,/,5x,
     1             'number of points          = ',i4,/,5x,
     2             'number of basis functions = ',i4)
      end






