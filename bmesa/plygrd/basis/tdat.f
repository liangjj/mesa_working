*deck tdat.f
c***begin prologue     tdat
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
c***end prologue       tdat
      subroutine tdat(card,qtyp,tl,tr,nsubg,npt,tgrd,maxd)
      implicit integer (a-z)
      character*(*) card, qtyp
      character*80 chrkey
      real*8 fpkey, tl, tr
      dimension npt(*)
      common/io/inp, iout
      write(iout,1) tgrd
      qtyp=chrkey(card,'quadrature-type','legendre',' ')
      nsubg=intkey(card,'number-of-subgrids',1,' ')
      tl=fpkey(card,'tl',0.d0,' ')
      tr=fpkey(card,'tr',1.d0,' ')
      call intarr(card,'number-of-points-per-subgrid',npt,nsubg,' ')
      write(iout,2) qtyp, tl, tr, nsubg
      write(iout,3) (npt(i), i=1,nsubg)
      do 10 i=1,nsubg
         maxd=max(maxd,npt(i))
 10   continue   
      return
 1    format(/,5x,'data for time variable region = ',i3)
 2    format(/,5x,
     1             'quadrature type                     = ',a24,/,5x,
     2             'left time boundary                  = ',e15.8,/,5x,
     3             'right time boundary                 = ',e15.8,/,5x,
     4             'number of time sub-grids            = ',i4)
 3    format(/,5x,'number of points in each subgrid = ',10(i3,1x))
      end






