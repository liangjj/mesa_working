*deck @(#)pm1902.f	5.1 11/6/94
      program m1902
c***begin prologue     pm1902.f
c***date written       860826   (yymmdd)
c***revision date      11/6/94
c   28 may, 1994       rlm at lanl
c      adding electric field gradients
c***keywords           m1902, link 1902, one-electron, integrals, property
c***author             saxe, paul and martin, richard    (lanl)
c***source             @(#)pm1902.f	5.1 11/6/94
c***purpose            computes 1-e property integrals over a general 
c                      contraction scheme.
c***description
c       e0             compute the electric monople integrals.
c       e1             compute the electric dipole integrals.
c       e2             compute the electric quadrupole integrals.
c       e3             compute the electric octupole integrals.
c       e4             compute the electric hexadecapole integrals.
c       v0             compute the electrostatic potential integrals.
c       v1             compute the electric field integrals.
c       v2             compute the electric field gradient integrals.
c       fermi          compute the fermi contact integrals (dirac delta).
c       del            compute the momentum integrals.
c       mv             compute the mass-velocity integrals.
c
c***references
c
c***routines called
c
c***end prologue       pm1902.f
      implicit none
c     --- input variables -----
      integer as, ax, ay, az
      real*8 s, x, y, z
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer maxnbf,inp,iout, idum
      integer nprps(0:4)
      integer lmult,nderiv,nat,nbf,nnp,ngot,maxcor
      integer cgrid,ngrid,zan,c,top,nprim,ncont,ntypes,nbtype
      integer lenxyz,mxcont,maxl,ptprim,noprim,ptcont,nocont
      integer start,nocart,nobf,minmom,maxmom,mintyp
      integer nx,ny,nz,cont,ex,nmat,need,left,lens, maxatm
      integer i
      integer wpadti,wptoin,iadtwp
      parameter (maxnbf=2000)
      parameter (maxatm=2000)
      logical dolp,fermif,mvint,delint,logkey
      logical positn,grid,inau
      character*4096 ops
      character*4096 refops
      character*80 card
      character*16 bflabl(maxnbf)
      real*8 toang
      dimension need(10), ngot(10)
c
      common/io/inp,iout
c
      data nprps/1,3,6,10,15/
      save nprps
      pointer (ps,s(1)), (ps,as(1))
      pointer (px,x(1)), (px,ax(1))
      pointer (py,y(1)), (py,ay(1))
      pointer (pz,z(1)), (pz,az(1))
c
      call drum
      refops='hf=tcscf,ci=(nroots),mcscf,casscf,'//
     $     'properties=(e0,e1,e2,e3,e4,v0,v1,fermi,del,mv),'//
     $     'populations=(full,orbital,bond,charge),'//
     $     'transition,print=(populations=(del,mv,fermi,v0,v1,e0,'//
     $     'e1x,e1y,e1z,e2xx,e2yy,e2zz,e2xy,e2xz,e2yz))'
c
c
c     --- collect the options string.
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     --- determine the properties we have been requested to do.
c         first check the electric multipoles.
      lmult=-1
      if(logkey(ops,'properties=e0',.false.,refops)) lmult=0
      if(logkey(ops,'properties=e1',.false.,refops)) lmult=1
      if(logkey(ops,'properties=e2',.false.,refops)) lmult=2
      if(logkey(ops,'properties=e3',.false.,refops)) lmult=3
      if(logkey(ops,'properties=e4',.false.,refops)) lmult=4
c
c     --- now look for electrostatic properties.
      nderiv=-1
      if(logkey(ops,'properties=v0',.false.,refops)) nderiv=0
      if(logkey(ops,'properties=v1',.false.,refops)) nderiv=1
      if(logkey(ops,'properties=v2',.false.,refops)) nderiv=2
c
c     --- check for the fermi contact.
      fermif=.false.
      fermif=logkey(ops,'properties=fermi',.false.,refops)
c
c     --- check for the mass-velocity integrals.
      mvint=logkey(ops,'properties=mv',.false.,refops)
c
c     --- check for the momentum integrals.
      delint=logkey(ops,'properties=del',.false.,refops)
c
c     --- get the basis function labels.
      call iosys('read character "basis function labels" from rwf',
     $          -1,0,0,bflabl)
c
c     --- get lengths, dimensions, etc. needed for core allocation --
      call iosys('read integer "number of atoms" from rwf',1,nat,0,' ')
      call iosys('read integer "number of basis functions" from rwf',
     $           1,nbf,0,' ')
      nnp=(nbf+1)*nbf/2
c
c     --- if there is a $grid section in the input stream, then read the 
c         coordinates of the points, else evaluate the potential, etc.
c         at the nuclear centers.  note that when there are point charges
c         in the calculation, we don't want to evaluate properties at each
c         of their centers. the "ngrid" variable is used to distinguish
c         "real" atomic centers.
c
c         first find how much core we have to work with.
c      call getscm(0,z(1),ngot,'m1902: main',0)
c      maxcor=iadtwp(ngot)
c
c     the test points will be in z(cgrid).
      cgrid=1
      need(1)=wpadti(cgrid+3*maxatm)
      call getmem(need(1),py,ngot(1),'tmpgrid',0)
      grid=positn('$grid',card,inp)
      if(grid) then
         inau=logkey(ops,'geom=inau',.false.,' ')
         call iosys('read real angstrom/bohr from rwf',1,toang,0,' ')
         call gcoord(ngrid,y(cgrid),inau,toang) 
      else
         call iosys('read integer '//
     $              '"number of atoms with basis functions"'//
     $              ' from rwf',1,ngrid,0,' ')
      endif
c
c     --- allocate core for the atomic numbers and coordinates
c
      zan=cgrid+3*ngrid
      c=zan+nat
      need(2)=wpadti(c+3*nat)
      call getmem(need(2),px,ngot(2),'grid',0)
      call copy(y,x,3*ngrid)
      call getmem(-ngot(1),py,idum,'tmpgrid',idum)
      call iosys('read real "nuclear charges" from rwf',
     $           -1,x(zan),0,' ')
      call iosys('read real coordinates from rwf',-1,x(c),0,' ')
c
c     --- if grid points have not been read, then evaluate at the
c         atomic centers.
      if(.not.grid) then
         call scopy(3*ngrid,x(c),1,x(cgrid),1)
      endif
c
c     at this point we have x(cgrid),x(zan), and x(c) saved and ready
c     for further use.
c
c     --- retrieve basis set pointers.
c
      call bfpntr(nat,nprim,ncont,ntypes,nbtype,ptprim,noprim,ptcont,
     $            nocont,start,nocart,nobf,minmom,maxmom,mintyp,
     $            nx,ny,nz,cont,ex,need(3))
      call getmem(need(3),py,ngot(3),'basis',0)
c
c     the pointers come back and point to the variables,
c     ay(ptprim), ay(noprim), ay(nocont), ay(start), ay(nocart), ay(nobf),
c     ay(minmom), ay(maxmom), ay(mintyp), ay(nx), ay(ny), ay(nz),
c     y(cont), and y(ex)
c
      call bf(ptprim,noprim,ptcont,nocont,start,
     $        nocart,nobf,minmom,maxmom,mintyp,nx,ny,nz,
     $        cont,ex,y,ay)
c
c     --- main scratch space
c      s=iadtwp(top)
c
c     --- calculate one-electron multipole integrals.
      if(lmult.ge.0) then
c        determine room needed to hold the electric multipole blocks.
         nmat=0
         do 10 i=0,lmult
            nmat=nmat+nprps(i)
   10    continue
         need(4)=wptoin(nnp*nmat)
c
c        get the memory for the s matrix.
c
         call getmem(need(4),ps,ngot(4),'s',0)
         call mltpol(x(c),y(ex),y(cont),s,ay(ptprim),
     $               ay(noprim),ay(nocont),ay(ptcont),nat,nprim,
     $               ntypes,nbtype,nnp,ncont,ay(start),nbf,x(zan),
     $               ay(nocart),ay(nobf),ay(maxmom),ay(mintyp),
     $               ay(nx),ay(ny),ay(nz),ay(minmom),lmult,nmat,
     $               ops,refops,bflabl)
         call getmem(-ngot(4),ps,idum,'s',idum)
      endif
c
c     --- calculate electrostatic properties if desired.
      if(nderiv.ge.0) then
c
c
c        allocate core for the the electrostatic potential,etc.
         if(nderiv.eq.0) then
            lens=1
            need(4)=wptoin(nnp*ngrid)
            call getmem(need(4),ps,ngot(4),'s',0)            
         else if(nderiv.eq.1) then
            lens=lens+1+3
            need=wptoin((1+3)*nnp*ngrid)
            call getmem(need(4),ps,ngot(4),'s',0)            
         else if(nderiv.eq.2) then
            lens=lens+1+3+6
            need=wptoin((1+3+6)*nnp*ngrid)
            call getmem(need(4),ps,ngot(4),'s',0)
         endif
         dolp=.false.
         call efld(x(c),y(ex),y(cont),ay(ptprim),ay(noprim),
     $             ay(nocont),ay(ptcont),nat,nprim,
     $             ntypes,nbtype,nnp,ncont,ay(start),nbf,x(zan),
     $             ay(nocart),ay(nobf),ay(maxmom),ay(mintyp),
     $             ay(nx),ay(ny),ay(nz),ay(minmom),ops,refops,s,
     $             nderiv,dolp,bflabl,ngrid,x(cgrid),lens)
         call getmem(-ngot(4),ps,idum,'s',idum)
      endif
c
c     --- calculate one-electron fermi contact integrals.
c
      if(fermif) then
c
c        determine room needed to hold the fermi contact blocks.
c
         need=wptoin(nnp*ngrid)
         call getmem(need(4),ps,ngot(4),'s',0)
         call dirac(x(c),y(ex),y(cont),s,ay(ptprim),ay(noprim),
     $            ay(nocont),ay(ptcont),nat,nprim,
     $            ntypes,nbtype,nnp,ncont,ay(start),nbf,
     $            ay(nocart),ay(nobf),ay(maxmom),ay(mintyp),
     $            ay(nx),ay(ny),ay(nz),ay(minmom),ops,refops,
     $            bflabl,ngrid,x(cgrid))
         call getmem(-ngot(4),ps,idum,'s',idum)
      endif
c
c     --- calculate one-electron momentum related integrals.
      if(mvint.or.delint) then
c        determine room needed to hold the momentum blocks.
         need=wptoin(nnp*3)
         call getmem(need(4),ps,ngot(4),'s',0)
         call del(x(c),y(ex),y(cont),s,ay(ptprim),ay(noprim),
     $            ay(nocont),ay(ptcont),nat,nprim,
     $            ntypes,nbtype,nnp,ncont,ay(start),nbf,
     $            ay(nocart),ay(nobf),ay(maxmom),ay(mintyp),
     $            ay(nx),ay(ny),ay(nz),ay(minmom),ops,refops,bflabl)
         call getmem(-ngot(4),ps,idum,'s',idum)
      endif
c
c     --- exit gracefully.
      call chainx(0)
c
c
      stop
      end
