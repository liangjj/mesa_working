*deck @(#)wrttpe.f	1.1 9/7/91
c***begin prologue     wrttpe
c***date written       890404   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           wrttpe, link 6001, orbital file
c***author             schneider, barry (lanl)
c***source             m6001
c***purpose            orbital file setup
c***description        sets up iosys file for numerical orbitals
c*** 
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       wrttpe
      subroutine wrttpe(file,grdtyp,ncen,ncon,nreg,npnts,pntbuf,
     1                  nolst,aosym,group)
      implicit integer (a-z)
      parameter (dimpr=300 , dimcen=10)
      real *8 rloc, charge
      character * (*) file, aosym, group
      character *11 titl
      character *13 grdtyp
      logical logky
      dimension aosym(dimpr)
      common /rloc/ rloc(3,dimcen) , charge(dimcen)
      common/logprt/logky(2)
      common/ io/ inp, iout
      titl='"con array"'
      call iosys ('write character "file title" to rwf',0,0,0,titl)
      if (logky(2)) then
          call iosys ('open orbs as new',npnts*20,0,0,file)
      else
          call iosys ('open orbs as new on ssd',npnts*20,0,0,file)
          write (iout,500)
      endif
      call iosys ('write integer "orb dimen" to orbs',1,dimpr,0,' ')
      call iosys ('write integer "centr dimen" to orbs',1,dimcen,0,' ')
      call iosys ('write integer "no. centers" to orbs',1,ncen,0,' ')
      call iosys ('write real centers to orbs',3*dimcen,rloc,0,' ')
      call iosys ('write real "nuclear charges" to orbs',ncen,
     1             charge,0,' ') 
      call iosys ('write character group to orbs',0,0,0,group)
      call iosys ('write integer "no. regions" to orbs',1,nreg,0,' ')
      call iosys ('write character "grid type" to orbs',0,0,0,grdtyp)
      call iosys ('write integer "no. grid pts" to orbs',1,npnts,0,' ')
      call iosys ('write integer "point buffer" to orbs',1,pntbuf,
     1            0,' ')
      call iosys ('write integer "no. cont" to orbs',1,ncon,0,' ')
      call iosys ('write character symmetry to orbs',0,0,0,aosym)
      call iosys ('write integer "final pts" to orbs',1,nolst,0,' ')
      return
  500 format (/,5x,'orbital file to ssd')
      end
