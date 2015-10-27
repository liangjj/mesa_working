*deck @(#)pm990.f	1.2  7/30/91
      subroutine pm990(z,a)
c
c***begin prologue     m990
c***date written       yymmdd   (yymmdd)
c***revision date      890119   (yymmdd)
c
c  19 january  1989   bhl at llnl
c     only call the routine enforce if ci=no=enforce
c
c  02 february 1988   bhl at brl
c     read mcscf unconverged vectors if mcscf failed
c
c  11 august 1987     pws at lanl
c     1. fixed bug causing frozen core density to be 0 for diagonal
c     densities when also running transition density matrices.
c     2. made cutoff an input option, and set so if small all
c     configurations can be searched.
c
c  10 august 1987     pws at lanl
c     fixing up so don't try more than the number-of-configurations
c     roots.
c
c  25 february 1987   pws at lanl
c      adding transition density capabilities.
c
c   2 february 1987   modified by rlm at lanl to handle average natural
c                     orbitals with associated weights.
c   28 august 1986    modified by pws at lanl to not do no's for
c        mcscf case, since they have already been done in mcscf
c
c
c***keywords           one-particle density matrix, most important
c                      configurations, ci analysis
c***author             saxe, paul (lanl)
c***source             @(#)pm990.f	1.2   7/30/91
c
c***purpose            to print the most important configurations in a
c                     ci wavefunction and form the one-particle density
c                     matrix.
c
c***description
c       m990 will analyze the ci vectors from one or more
c       roots of a ci, as well as form and print the approximate
c       natural orbitals for each root. currently the following options
c       are recognized:
c
c       nroots=n         the number of roots to analyze.
c       nimportant=n     the number of important configurations to
c                        print, with a default of 15.
c       starting_row=n   the row in the drt to use as the top of the
c                        graph. this is used to pick out a smaller ci
c                        from a larger list.
c       print_nos        print the natural orbitals.
c       print_occ        print only the occupancies of the natural
c                        orbitals.
c       average=(n,m..)  do average natural orbitals from the densities
c                        of roots n,m,...
c       transition=(n,m) calculate one-particle transition density matrix
c                        between roots n and m.
c
c          at the minute, this link requires the following constants
c       from the iosys routine:
c
c       symmetries in ci,  norbs,  nbasis,  nrows,  nwks,  nlevs,
c       orbfrm
c
c       and requires the following files on the read-write file:
c
c       arc,  weight,  orbtbf,  iout,  bfsym,  nlwks,  b,  orbsym,
c       "xform vector",  "ci root n"
c
c       where n is the numbers of the ci roots to analyze. this link
c       places the following files on the read-write file:
c
c       "mo 1pdm n", "no vector n", "no occ n"
c
c       if average natural orbitals are evaluated, they are denoted
c       as vector 0.
c
c
c***references
c
c***routines called    (none)
c
c***end prologue       m990
c
      implicit integer (a-z)
c
      real*8 one
      real*8 z(*),cutoff
      real*8 fpkey
      character*4096 ops
      character*3 symlab(8)
      character refops*80, chrkey*16
      logical doavg,avgall,logkey
      logical iskey
      integer a(*)
c
      common /io/ inp,iout
c
c
      parameter (one=1.0d+00)
c
      data symlab /'a','b','c','d','e','f','g','h'/
      data refops/'no average print'/
c
 1001 format (//,t5,'transition density between roots',i4,' and',i4)
c     ----- open the read-write file, input and output -----
c
c
c     ----- recover the options string -----
c
      write(iout,*)'  m990: wfn and density codes'
c
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     ----- check for mcscf running an scf, in which case exit -----
c
      if (logkey(ops,'mcscf=hf',.false.,' ')) go to 5000
c
c     ----- get the constants for dividing up core -----
c
      call iosys('read integer "symmetries in ci" from rwf',
     $           1,nsym,0,' ')
      call iosys('read integer norbs from rwf',1,norbs,0,' ')
      call iosys('read integer "number of basis functions"'//
     $        'from rwf',1,nbf,0,' ')
      call iosys('read integer "number of drt functions" from rwf',
     $     1,nbfdrt,0,' ')
      call iosys('read integer nrows from rwf',1,nrows,0,' ')
      call iosys('read integer nwks from rwf',1,nwks,0,' ')
      call iosys('read integer nlevs from rwf',1,nlevs,0,' ')
      call iosys('read integer orbfrm from rwf',1,orbfrm,0,' ')
      nnp=norbs*(norbs+1)/2
      nnpbf=nbf*(nbf+1)/2
c
c     ----- check for other options -----
c
         nroots=intkey(ops,'ci=nroots',1,' ')
         nimprt=intkey(ops,'ci=nimportant',15,' ')
         srow=intkey(ops,'ci=starting_row',1,' ')
         cutoff=fpkey(ops,'ci=coefficient=cutoff',0.01d+00,' ')
         if(.not.iskey(ops,'ci=no=average',refops)) then
c           the average keyword does not exist.
            doavg=.false.
            avgall=.false.
         else if(chrkey(ops,'ci=no=average',' ',refops).eq.' ') then
c           the average keyword exists as a standalone.
            doavg=.true.
            avgall=.true.
         else
c           the average keyword exists as a replacement string.
            doavg=.true.
            avgall=.false.
         endif
c
c     ----- the number of important configurations must be no larger
c              than the number of configurations!
c
      nimprt=nwks
      nroots=min(nroots,nwks)
      nimprt=min(nimprt,nwks)
      nlarge=min(10000,nwks)
      if (cutoff.lt.0.0001d+00) nlarge=nwks
c
c     ----- check for transition density matrix -----
c


c
c     ----- allocate core for drt arrays, vectors, etc. --
c
      walk=1
      brkdn=walk+nimprt
      arc=brkdn+norbs
      weight=arc+4*nrows
      orbtbf=weight+4*nrows
      ref=orbtbf+norbs
      alpha=ref+norbs
      beta=alpha+norbs
      bfsym=beta+norbs
      symnum=bfsym+nbfdrt
      nlwks=symnum+nbfdrt
      orbsym=nlwks+nrows
      b=orbsym+norbs
      irowsv=b+nrows
      segsv=irowsv+nlevs
      pagesv=segsv+nlevs
      iwtsv=pagesv+nlevs
      traksv=iwtsv+nlevs
      jrowsv=traksv+nlevs
      jwtsv=jrowsv+nlevs
      hrowsv=jwtsv+nlevs
      hwt=hrowsv+nlevs
      harcsv=hwt+nlevs
      bftorb=harcsv+nlevs
      lgwk=bftorb+nbf
c
      coef=iadtwp(lgwk+nlarge)
      c=coef+nimprt
         s=c
      lgcoef=s+nwks
      acoef=lgcoef+nlarge
      bcoef=acoef+nlevs
      dave=bcoef+nlevs
         d1=dave
      natocc=d1+nnp
      eigvec=natocc+norbs
      t1=eigvec+norbs**2
      t2=t1+nbf**2
      aotono=t2+nbf**2
      occ=aotono+nbf**2
      aotomo=occ+nbf
      motono=aotomo+nbf**2
      smat=motono+nbf**2
      smhalf=smat+nnpbf
      tocc=smhalf+nnpbf
      eigval=tocc+nbf
      u=eigval+nbf
      triang=u+nbf**2
      top=wpadti(triang+nnpbf)
c
c     ----- get core, and then read in arrays -----
c
      call getscm(top,a,maxcor,'m904',0)
c
      call iosys('read integer arc from rwf',4*nrows,a(arc),0,' ')
      call iosys('read integer weight from rwf',4*nrows,a(weight),0,
     $           ' ')
      call iosys('read integer orbtbf from rwf',norbs,a(orbtbf),0,' ')
      call iosys('read integer iout from rwf',nbfdrt,a(bftorb),0,' ')
      call iosys('read integer bfsym from rwf',nbfdrt,a(bfsym),0,' ')
      call iosys('read integer nlwks from rwf',nrows,a(nlwks),0,' ')
      call iosys('read integer b from rwf',nrows,a(b),0,' ')
      call iosys('read integer orbsym from rwf',norbs,a(orbsym),0,' ')
      do 1 i=orbsym,orbsym+norbs-1
         a(i)=a(i)-1
    1 continue
c
      call rzero(z(dave),nnp)
c
c     ----- fix up bftorb in the case of an mcscf -----
c
c
c     ----- find the most important configurations in each root -----
c
c
      write (iout,2)
    2 format (/,'     configuration list ',/)
c
c
      do 99014 i=1,nwks
         z(coef-1+i)=.1
         a(walk-1+i)=i
99014 continue
c
      if (logkey(ops,'drt=p-and-q',.false.,' ')) then
         write(iout,*)' q-space configurations '
      end if
c
      srow=1
      call prwalk(ops,z(coef),a(walk),a(brkdn),a(arc),a(weight),
     $            nimprt,nrows,a(orbtbf),norbs,a(ref),a(alpha),a(beta),
     #            a(bfsym),a(symnum),nsym,nbfdrt,symlab,iout,srow)
c
      if (logkey(ops,'drt=p-and-q',.false.,' ')) then
         write(iout,*)' p-space configurations '
         srow=2
         call prwalk(ops,z(coef),a(walk),a(brkdn),a(arc),a(weight),
     $            nimprt,nrows,a(orbtbf),norbs,a(ref),a(alpha),a(beta),
     #            a(bfsym),a(symnum),nsym,nbfdrt,symlab,iout,srow)
      end if
c
c     ----- and exit with grace -----
c
 5000 continue
c
      call chainx(0)
c
c
      stop
      end
