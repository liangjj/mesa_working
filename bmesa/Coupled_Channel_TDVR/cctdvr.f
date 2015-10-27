*deck cctdvr.f 
c***begin prologue     cctdvr
c***date written       960718   (yymmdd)
c***revision date               (yymmdd)
c***keywords           time, dvr, orthogonal polynomial
c***                   bose-einstein
c***author             schneider, b. i.(nsf)
c***source             tprop
c***purpose            driver for solution of coupled, time dependent
c***                   schroedinger equation.  a dvr basis in space
c***                   and time is used to expand the wavefunction.
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       cctdvr
      program cctdvr
c
      implicit integer (a-z)
      parameter ( maxreg=100, mxch=50 )
      character*4096 ops
      character*2 itoc
      character*8 prtflg
      character*8 type
      character*80 cpass, chrkey, precon, qtyp
      character*24 coord, key, typpot, system
      character*800 card
      character*24 i0stat
      character*128 filbec, filham
      character*3 timkey, trapon
      logical prn, dollar, logkey
      logical itsolv, itdiag
      logical spac
      real*8 edge, fpkey 
      real*8 thresh, cnverg, eps
      real*8 energy 
      integer*8 pham, pvt
      integer*8 pjunk
      dimension pham(4), n(4), nphy(4), ngrds(4), edge(maxreg+1)
      dimension coord(4), key(4), prn(30)
      dimension type(2), typpot(4), nwham(4), nwpot(4)
      common/io/inp, iout      
      common/punch/pun
c
      call drum
      write(iout,*)
      call iosys ('read character options from rwf',-1,0,0,ops)
      call iosys ('read character "print flag" from rwf',-1,0,0,prtflg)
      pnch=logkey(ops,'punch=on',.false.,' ')
c
c         set spatial dimensionality of problem and program options 
c
      spdim=intkey(ops,'number-of-space-dimensions',1,' ')
      ntreg=intkey(ops,'number-of-time-regions',1,' ') 
      nc=intkey(ops,'number-of-channels',1,' ')
      dim=spdim+1
      timkey=chrkey(ops,'m6295=time','on',' ')
      system=chrkey(ops,'coordinate-system','cartesian',' ')
      pnch=logkey(ops,'punch=on',.false.,' ')
      itsolv=logkey(ops,'iterative-linear-system-solve',.false.,' ')
      type(1)='unknown'
      type(2)=chrkey(ops,'open-ham','unknown',' ')
      write(iout,1) spdim
c
      write(iout,2)
      call iosys ('read character "bec filename" from rwf',
     1             -1,0,0,filbec)
      call iosys ('open bec as '//type(1),0,0,0,filbec)
      call iosys('rewind all on bec read-and-write',0,0,0,' ')
      spac=logkey(ops,'no-spatial-hamiltonian',.false.,' ')
      i0stat=chrkey(ops,'driver','state-vector',' ')
      state=intkey(ops,'initial-state',0,' ')
      write(iout,3) spac, i0stat, state, ntreg
      state=state+1  
c
c     get all of the one-dimensional matrices needed to construct
c     the spatial part of the hamiltonian and associated quantities.
c
      n3d=1
      do 10 i=1,spdim
         coord(i)=chrkey(ops,'dimension-'//itoc(i),'x',' ')
         call basis(ops,cpass,card,i,pham(i),ngrds(i),coord(i),
     1              key(i),typpot(i),n(i),nphy(i),nwham(i),prn)
         n3d=n3d*nphy(i)
 10   continue
c
c             set linear system procedures
c      
      if(itsolv) then
         if( dollar('$gmres',card,cpass,inp) ) then
             call lindat(card,cpass,cnverg,thresh,eps,precon,nblck,
     1                   prn(14),filham,type(2))
         endif
      endif
c
c         Begin Time Calculation
c
      do 80 tim=1,ntreg
c
         coord(dim)=itoc(tim)
         len=length(coord(dim))
         coord(dim)='t'//coord(dim)(1:len)
c
c        do the same for the time coordinate as the spatial basis
c
         call basis(ops,cpass,card,dim,pham(dim),ngrds(dim),
     1              coord(dim),key(dim),typpot(dim),n(dim),
     2              nphy(dim),nwham(dim),prn)
c
         call psit(pham,pvt,spdim,dim,nc,nphy,key,typpot,spac,
     1             tim,i0stat,state,system,itsolv,
     2             cnverg,thresh,eps,precon,nblck,nwpot(dim),
     3             card,prn(12),hdiag,pnch)
 80   continue   
      call chainx(0)               
      stop
 1    format(/,20x,'time-dependent basis function code',//,20x,
     1             'number of spatial dimensions = ',i1)      
 2    format(/,15x,'calculation = solve time-dependent schrodinger'
     1             ' equation')
 3    format(/,5x,'time-dependent data',/,5x,
     1            'no spatial hamiltonian   = ',l1,/,5x,
     2            'driver                   = ',a24,/,5x,
     3            'state                    = ',i2,/,5x,
     4            'number of time intervals = ',i3)
      end









