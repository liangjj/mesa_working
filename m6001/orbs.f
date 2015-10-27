*deck @(#)orbs.f	1.1 9/7/91
c***begin prologue     m6001
c***date written       890404   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           m6001, link 6001, orbital decomposition
c***author             schneider, barry (lanl)
c***source             m6001
c***purpose            driver for numerical orbital tabulation
c***description        calculates contracted ao's on physical grid
c***                   for integral calculation
c*** 
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       m6001
      program orbs
      implicit integer (a-z)
      parameter (dimpr=300 , dimcen=10)
      logical logkey, logky, grdtyp, posinp
      character *4096 ops
      character *1600 card
      character *8 chrkey, cpass
      character *128 filkne, filgrd, filorb
      character *3 itoc, ans, aosym, group, ctype
      character *16 grdnam
      real *8 z, alf, cont, rloc, charge, anorm, pi, fpkey
      common a(1)
      dimension z(1)
      equivalence (z,a)
      common /memory / ioff
      common /io/ inp,iout
      common/chrpss/ filgrd, filorb, aosym(dimpr)
      common /aosi/ npr, ncon, nxyzc(dimpr,4), nprc(dimpr)
      common /aosr/ alf(dimpr), cont(dimpr)
      common /rloc/ rloc(3,dimcen), charge(dimcen)
      common/logprt/logky(2)
      common/ctype/ctype(0:3,0:3,0:3)
      dimension anorm(dimpr), iholdc(dimpr,2)
      equivalence (anorm,iholdc)
      data pi /3.14159265358979323846d+00/
      call drum
      call iosys ('read character options from rwf',-1,0,0,ops)
      call iosys ('read integer maxsiz from rwf',1,memmax,0,' ')
      call iosys ('read character "kohn data filename" from rwf',-1,0,0,
     1             filkne)
      call iosys ('open kohndt as old',0,0,0,filkne)
      logky(1)=logkey(ops,'print=m6001=orbitals',.false.,' ')
      logky(2)=logkey(ops,'no-ssd',.false.,' ')
      write (iout,2)
c---------------------------------------------------------------------c
c            initialize x y and z values for cartesian aos            c
c---------------------------------------------------------------------c
      call nlmxyz
c---------------------------------------------------------------------c
c                read in data                                         c
c---------------------------------------------------------------------c
      grdtyp=.false.
      grdnam='"trns grid"'
      if( posinp ('$orbs',cpass) ) then
          call cardin (card)
          grdtyp=logkey(card,'untransformed-grid',.false.,' ')
          grdnam='"trns grid"'
      endif
      if (grdtyp) then
          grdnam='"untrns grid"'
          write (iout,500)
      else
          write(iout,600)
          endif
          call iosys ('read character "grid filename" from rwf',-1,0,0,
     1                 filgrd)
c----------------------------------------------------------------------c
c               open grid file and get no. of points                   c
c               need four arrays of dimension pntbuf                   c
c----------------------------------------------------------------------c
      call iosys ('open grid as old',0,0,0,filgrd)
      call iosys ('read integer "no. grid pts" from grid',1,npnts,0,
     1            ' ')
      call iosys ('read integer "point buffer" from grid',1,
     1            pntbuf,0,' ')
      call iosys ('read character "orbital filename" from rwf',-1,0,0,
     1             filorb)
      call iosys ('read "number of atoms" from kohndt',1,ncen,0,' ')
      call iosys ('read real "nuclear charges" from kohndt',ncen,charge,
     1             0,' ')
      do 4 i=1,ncen
         call iosys ('read real "x-y-z atom-'//itoc(i)//'" from kohndt',
     1               3,rloc(1,i),0,' ')
    4 continue
      if ( posinp('$geom',cpass)) then
           call cardin(card)
           ncen=intkey(card,'no-centers',1,' ')
           do 3 i=1,ncen
              charge(i)= fpkey(card,'charge-atom-'//itoc(i),1.d+00,' ')
              call fparr(card,'x-y-z-atom-'//itoc(i),rloc(1,i),3,' ')
    3      continue
      endif
      if( posinp('$basis',cpass) ) then
          call cardin(card)
          group=chrkey(card,'group','c2v',' ')
          call locase(group,group)
          npr=intkey(card,'no-prim',1,' ')
          ncon=intkey(card,'no-contracted',1,' ')
          call intarr (card,'no-prim/cont',nprc,ncon,' ')
      endif
c----------------------------------------------------------------------c
c           calculate number of passes of grid file needed             c
c----------------------------------------------------------------------c
      nreg=npnts/pntbuf
      nolst=npnts-nreg*pntbuf
      if (nolst.ne.0) then
          nreg=nreg+1
      else
          nolst=pntbuf
      endif
c----------------------------------------------------------------------c
c                 read in the basis set information                    c
c----------------------------------------------------------------------c
      call rdbsis(ntot,nprmx,aosym)
      write (iout,440)
      do 80 i=1,ncen
         write (iout,450) (rloc(j,i),j=1,3), charge(i)
   80 continue
      write (iout,460) npr,ncon
      write (iout,470) (nprc(i),i=1,ncon)
      if (logkey(ops,'print=m6001=basis',.false.,' ')) then
          call wrtbas (ntot)
      endif
      if (ntot.ne.npr) call lnkerr('no. primitives error')
c----------------------------------------------------------------------c
c               get memory                                             c
c               at least enough for one contracted                     c
c                     function                                         c
c               store x, y, z, wt and function                         c
c----------------------------------------------------------------------c        
      minsze=7*pntbuf
      if (minsze.gt.memmax) call lnkerr('cannot get enough '//
     1                                  'memory:quit')
c----------------------------------------------------------------------c
c         make buffer a multiple of number of contracted functions     c
c----------------------------------------------------------------------c
      nwleft=memmax-6*pntbuf
      nfmax=nwleft/pntbuf
      nfmax=min(nfmax,ncon)
      write (iout,350) nfmax,pntbuf
      nwords=wpadti((6+nfmax)*pntbuf)
      call iosys ('write integer maxsiz to rwf',1,nwords,0,' ')
      call getscm (nwords,z,ngot,'m6001',0)
      grid=ioff
      pre=grid+4*pntbuf
      rsq=pre+pntbuf
      fmo=rsq+pntbuf
      write (iout,340) nwords
      call wrttpe(filorb,grdnam,ncen,ncon,nreg,npnts,pntbuf,nolst,
     1            aosym,group)
      call rzero (anorm,ncon)
      size=ncon*npnts
      call iosys ('create real "con array" on orbs',size,0,0,' ')
      call iosys ('create real '//grdnam//' on orbs',4*npnts,0,0,' ')
      nwds=0
      nwrite=0
      do 200 ireg=1,nreg
         noptrg=pntbuf
         if (ireg.eq.nreg) then
             noptrg=nolst
         endif
c----------------------------------------------------------------------c
c           read in this block of grid points                          c
c----------------------------------------------------------------------c
         call iosys ('read real '//grdnam//' from grid without '//
     1               'rewinding',4*noptrg,z(grid),0,' ')   
         call iosys ('write real '//grdnam//' to orbs without '//
     1               'rewinding',4*noptrg,z(grid),0,' ')
         call makecn (z(fmo),z(grid),z(pre),z(rsq),anorm,ncon,nfmax,
     1                noptrg,nwrite,nwds,ireg,logky(1))
  200 continue
      if (nwds.ne.size) call lnkerr('error in words to orbs')
      write (iout,140) nwrite
      write (iout,150) nwds
      call iosys ('rewind all on grid read-and-write',0,0,0,' ')
      call iosys ('rewind all on orbs read-and-write',0,0,0,' ')
      call iosys ('close grid',0,0,0,' ')
      call iosys ('close orbs',0,0,0,' ')
      write (iout,360) (anorm(ii),ii=1,ncon)
      call iosys ('write integer maxsiz to rwf',1,memmax,0,' ')
      call chainx(0)
      stop
    1 format(a80)
    2 format(//,20x,'***** m6001:numerical orbital tabulation program **
     1***')
  140 format (//,5x,'no. writes to orbs',2x,i4)
  150 format (//,5x,'no. words written to orbs',2x,i8)
  340 format (//,5x,'memory adjusted by',2x,i8,2x,'words to calculate',
     1' contracted aos')
  350 format (//,5x,'maximum no. contracted functions stored',1x,i4,/
     1        5x,'point array buffer size',1x,i7)
  360 format (/,5x,'normalization integrals',(/,5x,5e15.8))
  440 format (//,10x,'position and charge of atomic centers',//,7x,
     1               'x',17x,'y',17x,'z',17x,'charge')
  450 format (/,1x,e15.8,3x,e15.8,3x,e15.8,3x,e15.8)
  460 format (//,5x,'data on atomic orbitals',//,5x,'no. prim.',
     1 2x,i4,5x,'no. of cont.',2x,i4)
  470 format (//,5x,'no. primitives / contracted function',(/,10x,20(1x
     1 ,i2)))
  500 format(/,5x,'untransformed grid')
  600 format(/,5x,'transformed grid')
      end
