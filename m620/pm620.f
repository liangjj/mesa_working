*deck @(#)pm620.f	5.1  11/6/94
      subroutine pm620(z,a)
c***begin prologue     pm620.f
c***date written       880712
c***revision date      11/6/94
c   12 february, 1994  rlm at lanl
c      making compatible with mesa; removing restrictions on size
c
c***keywords
c***author             martin, richard (lanl)
c***source             @(#)pm620.f	5.1   11/6/94
c***purpose            driver routine for cox-williams code
c***description
c
c    this code is a modification of pdm88, a code written by d.e. williams
c    and available from QCPE.
c
c    find least squares optimum values for net atomic charges (monopoles),
c    dipole components, and/or quadrupole components from the electric
c    potential.
c    alternatively, optimum bond dipoles may be found.  if desired, the
c    bond dipoles may be restricted to lie along the bond direction.
c    atomic lone pair dipoles, restricted in direction, may be added
c    to the restricted bond dipole model.
c-----------------------------------------------------------------------
c    input variables:
c       nat  number of atoms 
c
c       nbond  =0  for atomic multipole calculations
c              =1  for unrestricted bond dipole calculations
c      the input atom list is converted to a list of bond centers, and
c      the three components of each bond dipole are found.
c      a bond site is generated for every two atoms less than 1.6a apart
c      except no bond centers between h's are considered
c      a check is made for exactly one bond to a given hydrogen
c      x-h bond lengths are assumed in the range 0.7-1.2a
c              =2  for restricted bond dipole calculations
c      as with nbond=1, except the dipoles are restricted to the bond 
c      direction when nbond=2 additional atomic dipoles can also be included
c      in the model of either trigonal or tetrahedral types.  in the trigonal
c      type the atomic dipole is directed along the bisector of the bond.
c      in the tetrahedral type the dipole is directed perpendicular to
c      a plane formed by three bonds of normalized distance.
c
c      nsave   normally zero, except for a restricted bond
c              dipole calculation, where it indicates
c              how many (if any) atom save lines are to be read.
c              nsave.ne.0 allows the possibility of combining a restricted 
c      lone pair dipole calculation with a bond dipole calculation.
c      example: if nbond=2 and nsave=2, two extra lines are read which
c      specify which atoms carry the lone pair dipole and whether these
c      atoms are of trigonal or tetrahedral type.
c
c      prnt  print flag 
c-----------------------------------------------------------------------
c  lines 3. multipole selection integers, read if nbond=0 (i10)
c      nat lines are read only if nbond=0.  the ten integers specify
c    monopoles q, three dipole components ux, uy, uz, and 6 quadrupole
c    moments qxx, qyy, qzz, qxy, qxz, and qzz for each of the nat atoms.
c    0=do not select, 1=select
c    note that qxx should never be selected, since it is treated as a
c    dependent variable.
c      if nbond is 1 or 2, dipoles or restricted dipoles are 
c    automatically selected.
c-----------------------------------------------------------------------
c  lines 4. atomic lone pair dipole information, read only if nbond=2
c           and nsave.gt.0
c      if itype=1, the lone pair dipole bisects a trigonal set of atoms.
c      if itype=2, the lone pair dipole trisects a pyramidal set of atoms.
c    col.  1- 5   itype       1 for trigonal, 2 for tetrahedral
c    col.  6-10      i1       number of vertex atom (input list sequence)
c    col. 11-15      i2       number of bonded atom
c    col. 16-20      i3       number of bonded atom
c    col. 21-25      i4       number of bonded atom (tetrahedral only)
c-----------------------------------------------------------------------
c  units conversion factors
c  1 au (energy)   = 2625.5211  kj/mol
c  1 au (dipole)   = 2.541765   debye
c  1 ea (dipole)   = 4.803241   debye
c
c  variable            definition
c
c  ident(i,j)  i=1 has label number of first bond defining atom
c              i=2 has label number of second bond defining atom
c  isel(i,l)   selection integers for zm (0=nosel, 1=sel)
c  jz(i)       nuclear charge, used for shorten routine
c                should be 1 for hydrogen, 0 for non-nuclear sites
c  nbegin      normally 1, to indicate beginning of atom list or bond list
c                if greater than 1, indicates the beginning of a bond list
c                following a saved atom list
c  label(j)    symbol for element (h,c,n,o,f,si,p,s,cl) followed
c                by sequence number in input
c  site(m,i)   x,y,z (cartesian) for the ith site, m=1,3
c  vo(m,j)     x,y,z,vo for the jth grid point, m=1,4
c                converted to cartesian during input
c  zl(i,j)     direction cosines of bond j
c  zm(i,l)     lth multipole (variable) of the ith site
c  zmag(i)     magnitude of dipole on atom i
c  amat(jk)    ls matrix (upper triangle with initial work vector),
c                converted to inverse and correlation matrix later
c  xsdev(i,l)  esd's of zm
c
c
c***references
c
c***routines called
c
c***end prologue       pm620.f
      implicit none
c     --- input variables -----
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
      integer a(*)
      real*8 z(*)
c     --- local variables ---
      integer maxatm,nat,nbond,nsave,nsite,npoints,nq,nbegin
      integer i,nbf,nnp,ndmat,ngot,ioff
      integer inp,iout
      integer intkey
      integer vnuc,d,dtot,cgrid,top,poten,scr
      integer zm,zmag,avec,amat,der1,sdev,xsdev,site,tsite
      integer zl,ident,jz,ndip,isel,vo
      integer wpadti,iadtwp,need
      parameter (maxatm=2000)
      character*4096 ops
      character*5 label(maxatm)
      character*16 chdum(2*maxatm)
      logical logkey,cinz
      logical prnt
      real*8 sdot,dipole(3)
      real*8 two
      real*8 toang,jph,jpcal,avog,etoesu,slight,emass,fines
      real*8 toev,tokcal
c
      parameter (two=2.0d+00)
c
      common/io/inp,iout
c
 1000 format(1x,'m620: cox-williams potential fit:')
 1010 format(x,3i5)
 1020 format(x,a5,i5,3f10.5)
 1030 format(5i5)
 1040 format(x,5i5)
c
c     --- recover the options ---
      call iosys('read character options from rwf',-1,0,0,ops)
c
      write(iout,1000)
c
c     --- get some constants
      call iosys('read real angstrom/bohr from rwf',1,toang,0,' ')
      call iosys('read real j/hartree from rwf',1,jph,0,' ')
      call iosys('read real j/cal from rwf',1,jpcal,0,' ')
      call iosys('read real avogadro from rwf',1,avog,0,' ')
      call iosys('read real esu/e- from rwf',1,etoesu,0,' ')
      call iosys('read real lightspeed from rwf',1,slight,0,' ')
      call iosys('read real "electron mass" from rwf',1,emass,0,' ')
      call iosys('read real "fine-structure" from rwf',
     $            1,fines,0,' ')
      toev=(1.0d-05)*fines*fines*emass*slight*slight*slight/etoesu
      tokcal=jph*avog/(1.0d03*jpcal)
c
c     --- set up some input parameters
      call iosys('read integer "number of atoms" from rwf',1,nat,0,' ')
      nbond=intkey(ops,'cox-williams=nbond',0,' ')
      if(nbond.lt.0.or.nbond.gt.2) then
         write(iout,*) 'nbond=',nbond
         call lnkerr('bad value of nbond found for cox-williams')
      endif
      if(nbond.eq.2) then
         nsave=intkey(ops,'cox-williams=nsave',0,' ')
      endif
c     --- set print flag
      prnt=logkey(ops,'cox-williams=print-potential',.false.,' ')
      call iosys('read integer "number of solvent surface points"'
     $         //' from rwf',1,npoints,0,' ')
c
c     --- get the atomic site labels and transfer them into the local
c         label array ---
      cinz=.not.logkey(ops,'geom=coord',.false.,' ')
      if(cinz) then
         call iosys('read character z-names from rwf',-1,0,0,chdum(1))
      else
         call iosys('read character "z-names w/o dummies" from rwf',
     $              -1,0,0,chdum(1))
      endif
      do 10 i=1,nat
         label(i)=chdum(i)
   10 continue
c
c     --- allocate some core
      call getscm(0,z,ngot,'m620',0)
      site=1
      zl=site+3*nat
      ident=wpadti(zl+3*nat)
      isel=ident+2*nat
      top=isel+10*nat
      jz=top
      ndip=jz+nat
      tsite=iadtwp(ndip+5*nsave)
      need=wpadti(tsite+3*nat)
      if(need.gt.ngot) then
         call lnkerr('m620: need more core')
      endif
c     --- read the atomic coordinates and convert to angstroms ---
      call iosys('read real coordinates from rwf',-1,z(site),0,' ')
      call sscal(3*nat,toang,z(site),1)
c     --- get the atomic numbers
      call iosys('read integer "atomic numbers" from rwf',
     $            nat,a(jz),0,' ')
c
c     --- set the site selections flags,etc.
      call vary(ops,nbond,nat,a(ndip),nsave,nsite,nq,nbegin,a(jz),
     $          z(site),z(zl),z(tsite),a(ident),label,a(isel))
c
c     --- output site listing
      call sitout(z(site),z(zl),nsite,a(ident),nbegin,label,nbond)
c
c     --- determine the potential on the grid points --- 
      call iosys('read integer "number of basis functions" from rwf',
     $           1,nbf,0,' ')
      nnp=nbf*(nbf+1)/2
      call iosys('read integer "number of hf density matrices"'
     $         //' from rwf', 1,ndmat,0,' ')
c     --- allocate some core
      vo=iadtwp(top)
      vnuc=vo+4*npoints
      top=wpadti(vnuc)
      d=vnuc+npoints
      scr=d+ndmat*nnp
      dtot=scr+nbf*nbf
      poten=dtot+nbf*nbf
      cgrid=poten+nnp
      need=wpadti(cgrid+3*npoints)
      if(need.gt.ngot) then
         call lnkerr('m620: need more core')
      endif
c      
      call iosys('read real "solvent surface coordinates"'
     $         //' from rwf',npoints*3,z(cgrid),0,' ')
c     --- convert to angstroms
      do 250 i=1,npoints
         z(vo+4*(i-1))  =z(cgrid+3*(i-1))*toang
         z(vo+4*(i-1)+1)=z(cgrid+3*(i-1)+1)*toang
         z(vo+4*(i-1)+2)=z(cgrid+3*(i-1)+2)*toang
  250 continue
c     --- read in the density matrices, surface electrostatic potential
c         matrices, and generate potential 
c         first density matrix is the closed shell piece, second one open.
c         the closed shell piece must be doubled.
      call iosys('read real "hf density matrix" from rwf',
     $            ndmat*nnp,z(d),0,' ')
      call trtosq(z(dtot),z(d),nbf,nnp)
      call vadd(z(dtot),z(dtot),z(dtot),nbf*nbf)
      if(ndmat.eq.2) then
         call trtosq(z(scr),z(d+nnp),nbf,nnp)
         call vadd(z(dtot),z(dtot),z(scr),nbf*nbf)
      endif
      call iosys('read real "solvent surface nuclear electrostatic'
     $         //' potentials" from rwf',npoints,z(vnuc),0,' ')
      ioff=0
      do 260 i=1,npoints
         call iosys('read real "solvent surface electrostatic'
     $            //' potential matrices" from rwf',
     $            nnp,z(poten),ioff,' ')
         call trtosq(z(scr),z(poten),nbf,nnp)
         z(vo+4*(i-1)+3)=z(vnuc+i-1)+sdot(nbf*nbf,z(scr),1,z(dtot),1)
         ioff=ioff+nnp
  260 continue
c
c     --- convert potentials from au to kj/mol
      do 280 i=1,npoints
         z(vo+4*(i-1)+3)=(jph*avog/1000.0d+00)*z(vo+4*(i-1)+3)
  280 continue
c
c     --- allocate core needed for solver
      zm=iadtwp(top)
      zmag=zm+10*nsite
      avec=zmag+nsite
      amat=avec+nq
      der1=amat+nq+nq*(nq+1)/2
      sdev=der1+nq
      xsdev=sdev+nq
      need=xsdev+10*nsite
      if(need.gt.ngot) then
         call lnkerr('m620: need more core')
      endif
c     --- solve linear equations to get fit to the potential
      call pdm88(z(site),z(vo),z(zm),a(isel),z(amat),z(avec),z(der1),
     $           z(zmag),z(zl),a(ident),dipole,z(sdev),z(xsdev),label,
     $           nsite,nq,npoints,nbond,nsave,prnt)
c
c     --- store the monopoles, dipoles, quadrupoles on the rwf
      if(nbond.eq.0) then
         call iosys('write real "cox-williams monopoles" to rwf',
     $              nsite,z(zm),0,' ')
         call iosys('write real "cox-williams dipoles" to rwf',
     $              3*nsite,z(zm+nsite),0,' ')
         call iosys('write real "cox-williams quadrupoles" to rwf',
     4              6*nsite,z(zm+4*nsite),0,' ')
      endif
c
c
      return
      end
