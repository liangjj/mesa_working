*deck %W% %G%
      subroutine pm511(z,a)
c
c***begin prologue     %M%
c***date written       930511  (yymmdd)
c***revision date      %G%
c   january 31,1995    rlm at lanl
c      moving grid generation to lower levels (kmtxguts,kmopen)
c   march 13, 1994     rlm at lanl
c      fixing bug associated with number of lebedev points for l=13
c***keywords           m511, link 511, density-functional, dft
c***author             martin, richard (lanl)
c***source             %W% %G%
c***purpose            solves the kohn-sham equations.
c***description
c
c
c***references
c                      p.hohenberg and w.kohn,phys.rev.b 136,864(1964).
c                      w.kohn and l.j.sham,phys.rev.a 140,1133(1965).
c
c***routines called
c***end prologue       %M%
c
      implicit none
c     --- input variables -----
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
      integer a(*)
      real*8 z(*)
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer i
      integer nbf,nnp,multip,natoms,nae,nbe,nirrep
      integer nprim,ncont,ntypes,nbtype
      integer ptprim,noprim,ptcont,nocont,start
      integer nocart,nobf,minmom,maxmom,mintyp
      integer nx,ny,nz
      integer cont,ex,top,ztop
      integer inp,iout
      integer ncoul,nexch,nshell,ndmat
      integer mxiter,stdiis,finish,mxds
      integer intkey
      integer ian,lenxyz,ngrid,mxgrd,pmxgrd
      integer s,t,v,h,smhalf,f,mo
      integer c,zan,grid,wts,nradial,nomega
      integer rmax,lmax,atom,biggrd,maxr,pmaxr
      integer ptr,dsptr,shlnbf,shlmin,shlmax
      integer numso,lambda,orbsym,occsym,salc
      integer u,d,dlast,eigval,kept,temp,t1,t2,t3,t4,t5
      integer triang,jmat,jlast,kmat,focksv,error,diissv,energs
      integer alpha,beta,fcoef
      integer maxnbf,maxatm,maxrep,mxcont,maxl
      integer canget,values,wish,ntriang,valuesi
      integer left,kleft
      integer ptrad
      integer vwts,rnuc,amu,pwtx,rr,radii,akl
      integer dftcore,mxgbsiz,ngb,
     $        dfttmp,dftabs
      integer npf,nnprim,pstart,prtoao,dpr,jprim
      integer iadtwp,wpadti,itobyt
      integer ntotal,nonzer
      integer imaxl,bigl,charge
      integer angsiz,minesz,grdsiz
      integer shlblks,nnshl,qint,dijmax
      integer vorder,vlmax,vradial,vncrule,vomega
      
      parameter (maxatm=2000,maxnbf=2000,maxrep=14)
c
      character*8 prtflg
      character*16 bflabl(maxnbf)
      character*16 names(maxatm)
      character*16 grdtyp(maxatm)
      character*16 pgrdtyp(maxatm)
      character*8 symlabl(maxnbf)
      character*8 lirrep(maxrep)
      character*128 namint
      character*4096 ops
      character*8 calc
      character*6 exchf
      character*4 corrf
c
      real*8 fpkey,level
      real*8 zero,ten
      parameter (zero=0.0d+00,ten=1.0d+01)
      real*8 dencut,dmcut,defcut,defbf,defkm,bfcut,kmcut,toosmall
      parameter (defcut=1.0d-16,defbf=1.0d-10,defkm=1.0d-10,
     $           toosmall=1.0d-50)
c
      logical prnt,logkey,dodiis,pulay,page,usesym,gscf
      logical slater,becke,vwn,lyp,cnull
      logical dirj,dograd,poisj
      logical adjust
c
      common /io/     inp,iout
c
      data prnt/.true./
      data usesym/.false./
      save prnt,usesym
c
 1000 format(1x,'m511:')
 1010 format(5x,'memory available(bytes)',12x,i9)
 1011 format(5x,'diis convergence:                   ',1pe8.1)
 1012 format(5x,'energy convergence:                 ',1pe8.1)
 1015 format(5x,'maximum block size',17x,i9)
 1016 format(5x,'exchange-correlation functional: ',a6,'-',a4)
 1017 format(8x,'density-matrix cutoff:           ',1pe8.1)
 1018 format(8x,'density cutoff:                  ',1pe8.1)
 1019 format(8x,'basis function cutoff:           ',1pe8.1,
     $      /8x,'k-matrix cutoff:                 ',1pe8.1)
 1020 format(8x,'voronoi cells are size adjusted')
 1021 format(5x,'all integrals held in core.')
 1030 format(5x,'# integral triangles in core',i4)
 1040 format(5x,'need(bytes)',18x,i9)
 1050 format(5x,'level-shift used:',f5.2)
 1060 format(8x,'grid size; atom',i3,':',a8,8x,i6,2x,i3,' blocks')
 1070 format(5x,80a1)
c
c     --- get max core available in integers
      call getscm(0,z,canget,'m511:',0)
c
c     ----- recover the options string -----
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     ----- has printing been turned off externally? -----
      call iosys('read character "print flag" from rwf',-1,0,0,prtflg)
      if(prtflg.eq.'minimum') prnt=.false.
c
c     ----- check on options and set defaults -----
      mxiter=intkey(ops,'scf=cycles',30,' ')
      pulay=logkey(ops,'scf=pulay',.true.,' ')
      page=logkey(ops,'scf=page',.false.,' ')
      if(page) pulay=.false.
      gscf=logkey(ops,'scf=gscf',.false.,' ')
      dodiis=logkey(ops,'scf=diis',.true.,' ')
      stdiis=intkey(ops,'scf=startdiis',-4,' ')
      mxds=intkey(ops,'scf=maxdiis',5,' ')
      mxgrd=intkey(ops,'scf=mxgrid',8000,' ')
      if(.not.dodiis) mxds=0
      finish=intkey(ops,'scf=convergence',8,' ')
      level=fpkey(ops,'scf=level-shift',zero,' ')
      usesym=.not.logkey(ops,'sym=off',.false.,' ')
      dirj=logkey(ops,'scf=directj',.false.,' ')
      poisj=logkey(ops,'scf=poissonj',.false.,' ')
c     --- poisson parameters
      vorder=intkey(ops,'poisson=spline-order',3,' ')
      vorder=vorder+1
      vlmax=intkey(ops,'poisson=lmax',4,' ')
      vradial=intkey(ops,'poisson=nradial',51,' ')
      vncrule=intkey(ops,'poisson=rule',5,' ')
      if (dirj .and. poisj) call lnkerr(
     $     'incompatible options scf=directj and scf=poissonj')
c     --- quadrature options ---
      rmax=intkey(ops,'scf=radgrid',51,' ')
      lmax=intkey(ops,'scf=lebord',23,' ')
      dmcut=fpkey(ops,'scf=denmat-cutoff',defcut,' ')
      dencut=fpkey(ops,'scf=density-cutoff',toosmall,' ')
      bfcut=fpkey(ops,'scf=bf-cutoff',defbf,' ')
      kmcut=fpkey(ops,'scf=km-cutoff',defkm,' ')
      minesz=intkey(ops,'scf=minesz',100,' ')
      adjust=logkey(ops,'scf=adjustcell',.false.,' ')
c
c     --- pick a functional set (exchange,correlation)
      slater=logkey(ops,'scf=exchf=slater',.false.,' ')
      becke=logkey(ops,'scf=exchf=becke',.false.,' ')
      lyp=logkey(ops,'scf=corrf=lyp',.false.,' ')
      vwn=logkey(ops,'scf=corrf=vwn',.false.,' ')
      cnull=logkey(ops,'scf=corrf=null',.false.,' ')
      if (becke .and. slater)
     $     call lnkerr('m511: two correlation functionals chosen')
      if (lyp .and. vwn)
     $     call lnkerr('m511: two exchange functionals chosen')
c     have to pick at least one, make slater-vwn default
      if (.not.(becke .or. slater)) slater=.true.
      if (.not.(vwn.or.lyp.or.cnull)) vwn=.true.
      if(slater) then
         exchf='slater'
      else if(becke) then
         exchf='becke '
      endif
      if(vwn) then
         corrf='vwn'
      else if(lyp) then
         corrf='lyp'
      else
         corrf='null'
      endif
c
c     --- print some basic information
      if(prnt) then
         write(iout,1000)
         write(iout,1010) canget*itobyt(1)
         write(iout,1011) ten**(-finish)
         write(iout,1012) ten**(-finish-1)
         write(iout,1016) exchf,corrf
      endif
c
c     ----- get some basic information -----
      call iosys('read integer "number of atoms" from rwf',
     $           1,natoms,0,' ')
      call iosys('read integer "number of basis functions" from rwf',
     $           1,nbf,0,' ')
      nnp=(nbf+1)*nbf/2
      call iosys('read integer "number of primitive functions"'
     $           //' from rwf',1,npf,0,' ')
      nnprim=(npf+1)*npf/2
      call iosys('read integer "spin multiplicity" from rwf',
     $           1,multip,0,' ')
      call iosys('read integer "number of alpha electrons" from rwf',
     $           1,nae,0,' ')
      call iosys('read integer "number of beta electrons" from rwf',
     $           1,nbe,0,' ')
      if (usesym)
     $call iosys('read integer "number of irreducible representations"'
     $           //' from rwf',1,nirrep,0,' ')
c
      call iosys('read character "atomic grid name" from rwf',
     $           -1,0,0,grdtyp)
      call iosys('read character "z-names w/o dummies" from rwf',
     $           -1,0,0,names)
c
      if (nbf.gt.maxnbf) then
         call lnkerr('character core not long enough')
      end if
c
c     ----- set up some parameters depending on multip -----
      if (gscf.or.nbe.eq.0) then
         calc='general'
         if(nbe.eq.0) then
            nshell=1
         else
            nshell=2
         end if
c        get the general scf information from input
         nshell=intkey(ops,'scf=nshell',nshell,' ')
         ncoul=intkey(ops,'scf=ncoul',nshell,' ')
         nexch=intkey(ops,'scf=nexch',nshell,' ')
         ndmat=intkey(ops,'scf=ndmat',nshell,' ')
c        only the occupied shells are required in input,
c        add one for the virtuals
         nshell=nshell+1
      else
         if(multip.eq.1) then
            calc='closed'
            nshell=2
            ncoul=1
            nexch=1
            ndmat=1
         else
            calc='open'
            if(poisj) then
               call lnkerr('open shell poisson not ready')
            endif
            nshell=3
            ncoul=2
            nexch=2
            ndmat=2
         endif
      end if
c
      call iosys('write integer "number of shells" to rwf',1,nshell,
     $            0,' ')
      call iosys('write integer "number of hf density matrices" '//
     $            'to rwf',1,ndmat,0,' ')
c
c     --- retrieve symmetry and shell information ---
c     allocate core
      shlnbf=1
      shlmin=shlnbf+nshell
      shlmax=shlmin+nshell
      numso=shlmax+nshell
      lambda=numso+nirrep
      orbsym=lambda+nirrep
      occsym=orbsym+nbf
      alpha=iadtwp(occsym+nirrep*nshell)
      beta=alpha+nshell*nshell
      fcoef=beta+nshell*nshell
      salc=fcoef+nshell
      top=wpadti(salc+nbf*nbf)
      if (top .gt. canget)
     $     call lnkerr('m511: not enough core for setup')
      call setup(a(shlnbf),a(shlmin),a(shlmax),nshell,nbf,nae,nbe,
     $           ops,calc,a(occsym),a(numso),a(lambda),
     $           lirrep,nirrep,usesym,a(orbsym),z(salc),
     $           z(alpha),z(beta),z(fcoef))
c
c     --- get atomic numbers and coordinates
      ian=top
      c=iadtwp(ian+natoms)
      grdsiz=wpadti(c+3*natoms)
      top=grdsiz+natoms
      ztop=iadtwp(top)
      call iosys('read integer "atomic numbers" from rwf',
     $            -1,a(ian),0,' ')
      call iosys('read real coordinates from rwf',-1,z(c),0,' ')
c
c     --- generate grid points and weights ---
      if (lmax.lt.3 .or. lmax.gt.29)
     $     call lnkerr('m511: invalid lebedev order requested')
      nomega=angsiz(lmax)
c
c     --- mxgrd is defaulted to a typical standard grid size.  
      mxgrd=max(mxgrd,(rmax-1)*nomega)
      left=canget-top
c
c     --- generate the grid once, finding the actual grid sizes 
c         and the largest one. this is the grid used in kmatrix.
c         appropriate pointers are returned and mxgrd is redefined.
      call getgrid(z,a,left,natoms,a(ian),
     $             z(c),ops,rmax,lmax,grdtyp,adjust,
     $             mxgrd,ngrid,grid,wts,vwts,ptrad,maxr,top)
      call iosys('write real "external grid" to rwf',
     $           3*mxgrd*natoms,z(grid),0,' ')
      call iosys('write real "external grid weights" to rwf',
     $           mxgrd*natoms,z(wts),0,' ')
      call iosys('write integer "external atomic grid size" to rwf',
     $           natoms,a(ngrid),0,' ')
      call icopy(a(ngrid),a(grdsiz),natoms)
c
c     --- if we're doing poisson calculation, generate the
c         poisson grid.
      if(poisj) then
         do 10 i=1,natoms
            pgrdtyp(i)='general'
   10    continue
         top=wpadti(ztop)
         vomega=angsiz(2*vlmax+1)
         pmxgrd=max(mxgrd,vradial*vomega)
c        note that i bump the number of radial points by one so as to
c        be consistent with gofish2 -- when this is at it should be
c        this should be changed.
         call getgrid(z,a,left,natoms,a(ian),
     $                z(c),ops,vradial+1,2*vlmax+1,pgrdtyp,adjust,
     $                pmxgrd,ngrid,grid,wts,vwts,ptrad,pmaxr,top)
         call iosys('write real "poisson grid" to rwf',
     $              3*pmxgrd*natoms,z(grid),0,' ')
         call iosys('write real "poisson grid weights" to rwf',
     $              pmxgrd*natoms,z(wts),0,' ')
         call iosys('write integer "poisson atomic grid size" to rwf',
     $              natoms,a(ngrid),0,' ')
c         call iosys('write real "poisson voronoi weights" to rwf',
c     $              pmxgrd*natoms,z(vwts),0,' ')
         call iosys('write integer "poisson radial pointers"'
     $              //' to rwf',pmaxr*natoms,a(ptrad),0,' ')
      endif
c
c     --- now reallocate core based on what we know about grid sizes.
c     mxgrd is the largest atomic grid size for k-matrix, while
c     pmxgrd is the largest poisson grid size.
c     should pass mxgrd, pmxgrd to each routine independently
c     probably a mistake to allocate core here.
      biggrd=max(mxgrd,pmxgrd)
      ptrad=wpadti(ztop)
      grid=ptrad+max(maxr,pmaxr)+1
      wts=grid+3*biggrd
      vwts=wts+biggrd
      rnuc=vwts+biggrd
      amu=rnuc+natoms*natoms
      pwtx=amu+natoms*natoms
      rr=pwtx+natoms
      radii=rr+natoms
      akl=radii+natoms
      top=wpadti(akl+natoms*natoms)
      if(top.gt.canget) then
         write(iout,*) 'top,canget',top,canget
         call lnkerr('m511: not enough core to do grid') 
      endif
c
c     --- retrieve basis set information; returns pointers as well.
      call iosys('length of exponents on rwf',nprim,0,0,' ')
      call iosys('length of "contraction coefficients" on rwf',
     $            ncont,0,0,' ')
      call iosys('length of "number of pure functions" on rwf',
     $           ntypes,0,0,' ')
      call iosys('read integer "number basis types" from rwf',
     $           1,nbtype,0,' ')
      call iosys('length of "power of x" on rwf',lenxyz,0,0,' ')
      call iosys('read integer maxcont from rwf',1,mxcont,0,' ')
      call iosys('read integer maxl from rwf',1,maxl,0,' ')
      call iosys('read integer "number of shell blocks" from rwf',
     $            1,shlblks,0,' ')
      nnshl=shlblks*(shlblks+1)/2
      call basis(natoms,nbf,nprim,ncont,ntypes,nbtype,lenxyz,
     $           mxcont,maxl,ptprim,noprim,ptcont,nocont,start,
     $           nocart,nobf,minmom,maxmom,mintyp,nx,ny,nz,
     $           cont,ex,top,z,a)
      pstart=top
      prtoao=iadtwp(pstart+ntypes*natoms)
      imaxl=wpadti(prtoao+npf*nbf)
      charge=iadtwp(imaxl+natoms)
      top=wpadti(charge+2*natoms)
c
c     --- get additional basis set information
      call iosys('read integer "pointer to first primitive" from rwf',
     $           -1,a(pstart),0,' ')
      call iosys('read real t(prim,cont) from rwf',
     $            npf*nbf,z(prtoao),0,' ')
c
c     --- allocate core for kohn-sham solver ---
      s=iadtwp(top)
      t=s+nnp
      v=t+nnp
      h=v+nnp
      smhalf=h+nnp
      f=smhalf+nnp
      zan=f+nnp
      u=zan+natoms
      d=u+nbf**2
      dlast=d+nnp*ndmat
      dpr=dlast+nnp*ndmat
      jprim=dpr+nnprim*ndmat
      qint=jprim+nnprim*ndmat
      dijmax=qint+nnshl
      mo=dijmax+nnshl
      eigval=mo+nbf**2
      ptr=wpadti(eigval+nbf)
      dsptr=ptr+mxiter
      kept=dsptr+mxds
      temp=kept+nshell
      t1=iadtwp(temp+2*nbf)
      t2=t1+max(nbf**2,nnp*ncoul)
      t3=t2+max(nbf**2,nnp*nexch)
      t4=t3+nbf**2
      t5=t4+nbf**2*nexch
      triang=max(t5+nbf**2,t1+mxds**2+mxds)
      jmat=triang+nnp
      jlast=jmat+nnp*ncoul
      kmat=jlast+nnp*ncoul
      error=kmat+nnp*nexch
      focksv=error+nbf**2*mxds
      energs=focksv+nnp*mxds
      diissv=energs+mxiter
      ntotal=wpadti(diissv+mxiter)
      nonzer=ntotal+mxiter
      valuesi=nonzer+mxiter
      values=iadtwp(valuesi)
      left=iadtwp(canget)-values
c
c     --- finally we know how much we'll have left for kmatrix.
c     first, see if the gradient of the basis functions is needed
      if(becke.or.lyp) then
         dograd=.true.
c        increment maxl by 1 for gradients
         bigl=max(maxl+1,2)
      else
         dograd=.false.
         bigl=max(maxl,2)
      endif
c     it needs some arrays that are independent of the grid.
      if(calc.eq.'closed') then
c        (nnp for dtmp, nbf for phibar,
c         nbf for phigrad, minesz*nbf for scratch)
         dftabs=nnp+2*nbf+minesz*nbf
      else if (calc.eq.'open') then      
c        there are 2*nnp for dtmp.
c        and the others for phibar,gradbar(3),minesz*(nbf+1) scratch
         dftabs=2*nnp+4*nbf+(minesz+1)*nbf
      endif
      kleft=left-dftabs
c     --- the remaining core is used to hold arrays which depend
c         on the grid size. we are going to determine the largest
c         grid block we can handle with the available memory to 
c         ensure we can run in the amount of space provided.
c         this means we need to know how much space we need
c         PER GRID POINT. this is given by dftcore
c
      if(calc.eq.'closed') then
c        dengrida
         dftcore=1
c        fout
         dftcore=dftcore+5 
         if(lyp) then
c           fout needs an extra one
            dftcore=dftcore+1
         endif
c        becke/lyp needs even more
         if(becke.or.lyp) then
c           dengrada,ga,queue
            dftcore=dftcore+3+1+3
         endif
      else if(calc.eq.'open') then
c        dengrida,dengridb
         dftcore=2
c        fout
         dftcore=dftcore+5 
         if(lyp) then
            dftcore=dftcore+1
         endif
c        becke/lyp needs even more
         if(becke.or.lyp) then
c           dengrada,dengradb,ga,gb,queue
            dftcore=dftcore+3+3+1+1+3
            if(lyp) then
c              gab
               dftcore=dftcore+1
            endif
         endif
      endif
c     both open/closed need
c     phi,grad,scr,scr2,tea,tmpgwt,nzptrs(treat this as if were a real*8
c                                         array for now)
      dftcore=dftcore+nbf+3*nbf+1+1+nbf+1+1
c
c     --- this takes us up to the beginning of itch in kmatrix.
c         the functionals (except slater) need some scratch space
      dfttmp=0
      if (becke) dfttmp=5
      if (vwn) then
         dfttmp=max(dfttmp,16)
      else if (lyp) then
         dfttmp=max(dfttmp,12)
      endif 
c     --- need some scratch too for the direct basgrd calls 
c         it is overlapped with the functionals, however.
      dfttmp=max(dfttmp,1+2*mxcont+3*bigl)
c
c     ta da..... this is how much room we need PER GRID POINT
      dftcore=dftcore+dfttmp
c
c     --- finally, determine the maximum block size which will fit
      mxgbsiz=kleft/dftcore
      
c
c     let user override, to avoid using all the memory there is and 
c     choking the machine, but only to DECREASE it, not to increase it.
      mxgbsiz=min(intkey(ops,'scf=mxgbsiz',mxgbsiz,' '),mxgbsiz)
      mxgbsiz=min(mxgrd,mxgbsiz)
      if(prnt) then
         write(iout,1015) mxgbsiz
         do 20 atom=1,natoms
            ngb=a(grdsiz+atom-1)/mxgbsiz
            if(mod(a(grdsiz+atom-1),mxgbsiz).eq.0) then
c              just fits.
            else
               ngb=ngb+1
            endif
            write(iout,1060) atom,grdtyp(atom),a(grdsiz+atom-1),ngb
   20    continue
      endif
c
c     --- divvy up remaining core if doing conventional scf
      if (.not.dirj .and. .not. poisj) then
c        --- open the integral file ---
         call iosys('read character "integral filename" from rwf',
     $        0,0,0,namint)
         call iosys('open ints as old',0,0,0,namint)
c        --- want to use the remaining core to hold as many two-electron
c            integrals as possible.  the basic unit is a triangle of length
c            nnp.  first find out how much core is available, and 
c            then determine the number, ntriang, of triangles
c            which we can hold in core. we would like to get them
c            all in.
         wish=wpadti(values+nnp**2)
         if (wish.le.canget) then
            ntriang=nnp
         else
            ntriang=left/nnp
         end if
      end if
c
c      
      if(prnt) then
         write(iout,1017) dencut
         write(iout,1018) dmcut
         write(iout,1019) bfcut,kmcut
         if(adjust) write(iout,1020)
         if(.not.dirj .and. .not. poisj) then
            if(ntriang.eq.nnp) then
               write(iout,1021)
            else
               if(ntriang.ge.1) then
                  write(iout,1030) ntriang
               else
                  write(iout,1040) wpadti(values+nnp)*itobyt(1)
                  call lnkerr(' m511: need more core')
               endif
            endif
         endif
         if (level.ne.zero) write (iout,1050) level
      endif
c
c     --- read in s, t and v one-electron integrals ---
      call iosys('read real "overlap integrals" from rwf',
     $           nnp,z(s),0,' ')
      call iosys('read real "kinetic integrals" from rwf',
     $           nnp,z(t),0,' ')
      call iosys('read real "potential integrals" from rwf',
     $           nnp,z(v),0,' ')
c
c     --- recall s**(-1/2) ---
      call iosys('read real "overlap to minus half" from rwf',
     $           nnp,z(smhalf),0,' ')
c     --- do the scf iterations ---
      call ksham(z(s),z(t),z(v),z(h),z(smhalf),z(f),z(u),z(eigval),
     $          z(t1),z(t2),z(t3),z(triang),z(values),a(ian),z(c),
     $          z(mo),z(wts),nradial,nomega,natoms,nbf,nnp,z(d),
     $          z(dlast),z(t4),z(t5),z(error),z(focksv),mxiter,finish,
     $          stdiis,ntriang,z(t1),nshell,a(shlmin),a(shlmax),
     $          z(alpha),z(beta),dodiis,pulay,page,calc,
     $          ncoul,nexch,ndmat,z(energs),z(diissv),prnt,ops,
     $          mxds,a(ptr),a(dsptr),level,z(fcoef),
     $          z(jmat),z(jlast),z(kmat),bflabl,usesym,symlabl,z(salc),
     $          a(numso),a(lambda),lirrep,nirrep,a(orbsym),a(occsym),
     $          a(temp),a(kept),left,mxgrd,pmxgrd,a(grdsiz),dmcut,
     $          dencut,bfcut,kmcut,a(ptprim),a(noprim),nbtype,
     $          z(ex),a(nx),a(ny),a(nz),lenxyz,
     $          a(nocart),a(mintyp),a(maxmom),mxcont,nprim,
     $          a(pstart),z(prtoao),z(dpr),z(jprim),npf,nnprim,
     $          dirj,poisj,a(ntotal),a(nonzer),
     $          z(cont),a(nocont),a(ptcont),ntypes,ncont,a(start),
     $          a(nobf),a(minmom),z(grid),z(charge),a(imaxl),bigl,
     $          dograd,slater,becke,vwn,lyp,
     $          mxgbsiz,a(valuesi),names,
     $          grdtyp,
     $          z(vwts),z(rnuc),z(amu),z(pwtx),z(rr),z(radii),z(akl),
     $          a(ptrad),rmax,lmax,adjust,minesz,nnshl,z(qint),
     $          z(dijmax))
c
c     --- form the lagrangian ---
c         actually, for now just form the energy weighted density matrix
c         the scf vector is still in z(mo) and the eigenvalues are 
c         in z(eigval)
      if (calc .eq. 'closed') then
         call nrgwden(z(mo),z(eigval),z(d),nbf,nnp,nbe)
      endif
      
c     --- and exit gracefully ---
      call chainx(0)
c
c
      return
      end
