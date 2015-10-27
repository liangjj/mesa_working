*deck @(#)m801.f	5.1  11/6/94
      program m801
c
c***begin prologue     m801
c***date written       810101   (yymmdd)
c***revision date      900410   (yymmdd)
c
c      10 april 1990 rlm at lanl
c                    modifying so that the symmetry information can be
c                    picked up from the rwf.
c
c      7 september 1987  pws at lanl
c                    adding the ability to restrict configurations
c                    with a minumim and maximum number of open shells
c                    these options are drt=nopen=(min=n,max=m)
c
c     12 july 1987   mesa 2.0 version implemented on both the
c                    sun's and cray's at los alamos. the input is
c                    changed to only orbital codes, with everything
c                    else passed in as keywords.
c
c
c    25 april 1985   lanl
c                    mesax version created on cray x-mp
c                                                           pws
c
c     7 may 1985      lanl                                  pws
c                    modified to default to cisd in all orbitals if no
c                    input can be found. also, other aspects such as number
c                    of symmetry types and number of basis functions can be
c                    defaulted.
c
c     1 january 1984 pe 3252 version adapted to gould 8705, brl/ibd.
c                    1 problem in scan with function in if statement,
c                    and changed block data names so they differ from
c                    the common block names. version 6
c                                                           pws
c
c     28 june 1983   set up interacting option of 2 to restrict
c                    valence type refernces to same symmetry.
c                                                           pws
c
c     28 june 1983   modified paldus again so that can get correct
c                    interacting space for at least 3-in-3 type
c                    references                             pws
c
c     24 june 1983   added routine ref to calculate weight (number)
c                    of reference configuration(s) for single
c                    reference, multi-reference and valence type
c                    calculations. will not work for special
c                    calculations.                          pws
c
c     22 june 1983   modified paldus to calculate interacting space
c                    for valence type calculations, at least for
c                    closed shell type reference function   pws
c
c     31 march 1983  changed orbital order to handle ci gradients
c                    with restricted core and virtual orbitals.
c                                                           pws
c
c
c***keywords           drt
c***author             saxe, paul (lanl)
c***source             @(#)pm801.f	5.1   11/6/94
c
c***purpose            to construct the drt
c
c***description
c
c
c   distinct row table program for the sdguga ci system.
c     written by paul saxe, uc berkeley, 1981-1982
c
c   this program will read a flexible input format of orbital codes,
c   rearrange the orbitals to a form suitable for the ci and generate  *
c   all the arrays needed to describe the ci calculation and the
c   integral storage. the input consists of the following cards
c
c   occupation codes: (repeat count)(key)code symmetry
c          fzc  frozen core
c          fzv  frozen virtual
c          cor  restricted core (integrals are transformed)
c          vir  restricted virtual
c          doc  doubly occupied
c          uoc  virtual
c          alp  alpha occupancy (spin increase)
c          bet  beta occupancy  (spin decrease)
c          spe  special orbitals to be defined later in input
c
c          repeat count (optional) is number of identical orbitals
c          key is % for valence excitation orbitals
c                 / for orbitals differing in different references
c          symmetry is a number (1-8) identifying symmetry type,
c                n.b. for d2h mandatory to number in cottons way.
c                if the symmetry identifiers are assigned automatically
c                (by not specifying 'nsym=integer' in the drt option
c                of the route) then this field should be left blank.
c
c          example: fzc1 2%doc1 3 %uoc1 15uoc3  (blanks ignored)
c
c   2. extra codes for references greater than first. give only
c      codes for orbitals corresponding to those flagged with a '/'
c      in section 4! for example, three-configs from two-orbitals
c      could be as follows
c            4)  /doc1 3doc1 /uoc1 5uoc1
c            5)   uoc1 doc1     alp1 bet1
c   3. if there are special codes, the program needs to know how
c      to handle them. the possibilities are:
c        a. gvb  placing this in first three positions of this card
c                causes a two-reference interacting calculation to
c                be run.
c        b. oss  for open-shell singlet
c        c. mat  to enter matrix,etc.  next cards are
c            1. no. e- in special orbitals, spin*2, and total symmetry *
c                  special orbitals. (3i5)
c            2. matrix of excitations into orbitals, dimension 4**no.  *
c               special orbitals. entered as multi-dimensional array
c               in fortran, leftmost index giving case (1-4) for walk  *
c               for first special orbital, etc. (4(4i1,1x)).
c               example, for gvb pair: 0000 0111 0111 0112
c
c   the ci program can handle arbitrary reference sets,arbitrary
c   excitation levels, references such as trtiples in a selected space *
c   and say singles outside that space. there is essentially no limit  *
c   to the total number of unpaired electrons, orbitals or
c   configurations except for computer time avaliable.
c
c
c***references
c
c***routines called    (none)
c
c***end prologue       m801
c
      implicit integer (a-z)
c
      integer numint
      parameter (maxrep=14,maxnbf=2000)
      real*8 versin
c----------------------------------------------------------------
      real*8 z
      integer a
      pointer (p,z(1)), (p,a(1))
c---------------------------------------------------------------
      character*1 multrf,valenc
      character*1 c(maxnbf)
      character*16 cjunk
      character*4096 ops
      character*3 codes
      character*18 words
      character*8 lirrep(maxrep)
      character*4 subgrp
      character chrkey*16, inpflg*16, label*16
      character*8 dsk
      logical noinp,positn,no34x,symrwf,usesym
      logical logkey,spini,symi
c
      common /drtinf/ na,nb,ns,nespec,maxb,levfrm,levval,levopn,levmul
     #,               levocc,spec,sspesh,val
      common /drtcod/ ncodes,dela(9),delb(9),delele(9)
     #,               ntypes,virtul,occupd,valocc,rescor,resvir,frozen
     #,               valvir,opensh,multi,speshl
      common /drtchr/ codes(9),words(9),multrf,valenc
      common /dimens/ nbf,nsym,norbs,nrowsp,nrows4p,nrows,nrows4
     #,               nlevs,nrefs,nrowoc,nrow4o,nwks,nwksoc,nlevoc
     #,               orbfrm,symorb,numij,ngroup,numint,nmax,nspc,nvref
     #,               nijvir
      common /tapes/  out,errout,input,drttap
      common /versin/ versin
      common /bloksz/ blksiz,absmax,maxsiz
      common /io/     inp,ioutq
c
      save cortop
c
c     drt:'version 7; 25-april-1985.'
   30 format(1x,'m801:drt generation',
     $        /,5x,a16,
     $        /,5x,'number of symmetries',16x,i2,
     $        /,5x,'number of references',15x,i3,
     $        /,5x,'virtual excitation level',12x,i2,
     $        /,5x,'valence excitation level',12x,i2,
     $        /,5x,'integral block size',7x,i12)
   31 format(5x,'3-4 external arrays',16x,' no')
   32 format(5x,'3-4 external arrays',16x,'yes')
   33 format(5x,'interacting space',17x,'spin')
   34 format(5x,'interacting space',17x,'full')
   35 format(5x,'symmetry restrict valence space',4x,'yes')
   36 format(5x,'symmetry restrict valence space',4x,' no')
c
c__________________________________________________________________
      call drum
      input=inp
      out=ioutq
      errout=ioutq
c__________________________________________________________________
c     versions with same integer part are compatible,fraction is ddmmyy
      versin=7.250485
c
c     ----- recover the options string -----
c
      call iosys('read character options from rwf',-1,0,0,ops)
      call getmem(0,p,ngot,'first',0)
      call iosys('read integer mxcore from rwf',1,top,0,' ')
      call iosys('read integer mncore from rwf',1,bottom,0,' ')
      write(out,*) 'minimum core available = ',bottom
      write(out,*) 'maximum core available = ',top
      if(chrkey(ops,'core-use=m801','maximum',' ').eq.'minimum') then
         top=bottom
      endif 
      call getmem(top,p,cortop,'m801',1)
      dsk='rwf'
      call iosys('write character "drt file name drt" to rwf',0,0,0,dsk)
c
c     ----- start the timing routines ------
c
c
c
c      do 1 i=1,cortop
c         a(i)=0
c    1 continue
c
c     check for input flag replacement.
      inpflg='$'//chrkey(ops,'drt=repnam','drt',' ')
c
c     read the options string.
      label=chrkey(ops,'drt=label',' ',' ')
      exv=intkey(ops,'drt=vexcite',2,' ')
      exref=intkey(ops,'drt=refexcite',0,' ')
      no34x=logkey(ops,'drt=no34x',.false.,' ')
      symi=logkey(ops,'drt=ref=symlimit',.false.,' ')
c..bhl      spini=logkey(ops,'drt=space=spin',.true.,' ')
      spini=.true.
       if(logkey(ops,'drt=noint',.false.,' ')) spini=.false.
c..bhl
c     out=chrkey(ops,'drt=out','out',' ')
      blksiz=intkey(ops,'drt=blksiz',1000000,' ')
      absmax=blksiz
      levfrm=intkey(ops,'drt=fermilevel',0,' ')
      nbf=intkey(ops,'drt=nbf',0,' ')
      nrefs=intkey(ops,'drt=nrefs',1,' ')
      nsym=intkey(ops,'drt=nsym',0,' ')
      if(logkey(ops,'drt=rdsym',.false.,' ')) then
         if(nsym.eq.0) nsym=1
         symrwf=.false.
      else
         symrwf=.true.
      end if
c
c
      noinp=.not.positn(inpflg,cjunk,inp)
c
c     ----- fill in defaults if needed -----
c
c     has symmetry been turned off for some reason?
      call iosys('read integer usesym from rwf',1,usesym,0,' ')
      symrwf=symrwf.and.usesym
      if (symrwf) then
         call iosys('read character "group symbol" from rwf',0,0,
     $               0,subgrp)
         call iosys('read integer "number of irreducible '//
     $              'representations" from rwf',1,nsym,0,' ')
         call grpsym(subgrp,lirrep,nsym)
      end if
      if (nsym.eq.0) nsym=1
c
c     turn off the interacting space if a multireference calculation.
      if (nrefs.gt.1) spini=.false.
c
      if (nbf.eq.0) call iosys('read integer "number of basis '//
     $     'functions" from rwf',1,nbf,0,' ')
      call iosys('write integer "number of drt functions" to rwf',
     $     1,nbf,0,' ')
      call iosys('write integer "symmetries in ci" to rwf',1,nsym,0,' ')
      call iosys('write character "drt label" to rwf',0,0,0,label)
c
c
      ver=versin
      day=(versin-ver)*100
      mo=((versin-ver)*100-day)*100
      yr=(((versin-ver)*100-day)*100-mo)*100+0.5
c     write (out,30) ver,day,mo,yr,label,idum,exv,exval,spini,symi,blksiz,,nsym,
c    #,              nrefs
      write(out,30) label,nsym,nrefs,exv,exref,blksiz
      if(no34x) then
         write(out,31)
      else
         write(out,32)
      endif
      if(spini) then
         write(out,33)
      else
         write(out,34)
      endif
      if(symi) then
         write(out,35)
      else
         write(out,36)
      endif
c  30 format ('1',79('*'),//,t30,'distinct row table program version'
c    #,        i2,/,t40,i2,'/',i2,'/',i2,//,1x,26a3,//
c    #,       '   prnt   x  valx  int  blk  out ',/,1x,6i5,//,t10
c    #,       'number of symmetries, of basis functions, of references'
c    #, /,t10,'          ----------     ---------------     ----------'
c    #,     /,t25,i1,t38,i4,t57,i4)
c
c
      ntype=1
      bfnum=ntype+ntypes
      numsym=bfnum+nbf*ntypes
      bfsym=numsym+ntypes*nsym
      bfcode=bfsym+nbf
      mosym=bfcode+nrefs*nbf
      top=mosym+nbf
      if (top.le.cortop) go to 5
      write (errout,4)
    4 format (//,' not enough core space to even read in orbital ',/
     #,          ' information -- something is very wrong',//)
      call lnkerr('not enough core for orbital information')
    5 continue
c
c     allocate character core.
      bfkey=1
      if(nbf.gt.maxnbf) call lnkerr('not enough character core in m801')
c
c
      if(symrwf) then
         call iosys('read integer "transformation vector symmetries"'
     $             //' from rwf',nbf,a(mosym),0,' ')
      end if
      call scan(a(ntype),a(bfnum),a(numsym),a(bfsym),c(bfkey),
     #          a(bfcode),a(mosym),symrwf,noinp,spini)
c
      orbtbf=top
      iout=orbtbf+norbs
      orbsym=iout+nbf
      top=orbsym+norbs
      if (top.le.cortop) go to 7
      write (errout,6)
    6 format (//,' not enough core space to reorder orbitals -- ',/
     #,          '    something is very wrong',//)
      call lnkerr('not enough core to reorder orbitals')
    7 continue
      call reordr(a(ntype),a(bfnum),a(numsym),a(bfsym),c(bfkey)
     #,           a(bfcode),a(orbtbf),a(iout),a(orbsym))
c
c     ----- set fermi-level if desired -----
c
      if(levfrm.ne.0) orbfrm=levfrm-1
c
c
c
      nspc=0
      if (spec.gt.0) nspc=4**spec
      nlevs=norbs+1
      nlevoc=nlevs-levfrm+1
      levpt=top
      levnr=levpt+nlevs
      nelecs=levnr+nlevs
      spc=nelecs+nrefs
      top=spc+nspc
      nrowsp=(cortop-top)/(9+nrefs)
      nrows4p=nrowsp*4
      nopen=top
      ap=nopen+nrowsp
      bp=ap+nrowsp
      sp=bp+nrowsp
      arcp=sp+nrowsp
      nlwksp=arcp+nrows4p
      x=nlwksp+nrowsp
      if (top.le.cortop) go to 9
      write (errout,8)
    8 format (//,' not enough core to try to make drt -- something',/
     #,          ' is quite wrong',//)
      call lnkerr('not enough core to make drt')
    9 continue
      if(spec.gt.0) call getspe(a(spc),noinp)
      call paldus(ops,a(bfsym),c(bfkey),a(bfcode),a(orbtbf)
     #,           a(ap),a(bp),a(sp),a(levpt),a(levnr),a(arcp)
     #,           a(nlwksp),a(x),a(nelecs),a(spc),
     $            exv,exref,spini,symi,a(nopen))
      b=ap+nrows
      s=b+nrows
      arc=s+nrows
      nlwks=arc+nrows4
      top=nlwks+nrows
      call squeez(a(bp),a(sp),a(arcp),a(nlwksp),a(b),a(s),a(arc)
     #,           a(nlwks))
      wght=top
      wtab=wght+nrows4
      wtw=wtab+orbfrm
      wtx=wtw+nsym*orbfrm
      wty=wtx+nsym*orbfrm
      top=wty+orbfrm
      if (top.le.cortop) go to 11
      write (errout,10) top,cortop
   10 format (//,' not enough core to run weight, need',i7,' and have'
     #,           i7,' words of core available',//)
      call lnkerr('not enough core for wieghts')
   11 continue
      call weight(a(ap),a(b),a(s),a(arc),a(levnr),a(levpt),a(nlwks)
     #,           a(wtab),a(wtw),a(wtx),a(wty),a(wght))
      symorb=nsym*norbs
      numij=norbs*(norbs+1)/2
      nijvir=orbfrm*(orbfrm+1)/2
      kadd=top
      ladd=kadd+symorb
      ijadd=ladd+symorb
      ijgrp=ijadd+numij
      inint=ijgrp+numij
      inext=inint+norbs
      jmnnxt=inext+norbs
      jmxnxt=jmnnxt+norbs
      ijxx=jmxnxt+norbs
c     ijxx=ningrp+norbs
      klxx=ijxx+numij
      nklxx=klxx+nijvir
      ijww=nklxx+nsym*orbfrm
      klww=ijww+numij
      nklww=klww+nijvir
      ningrp=nklww+nsym*orbfrm
      top=ningrp
      if (top.le.cortop) go to 15
      write (errout,14) top,cortop
   14 format (//,' not enough core to complete calculation -- need',i7
     #,          ' and have',i7,' words of core',//)
      call lnkerr('not enough core to wrap up calculation')
   15 continue
      do 16 i=kadd,top
         a(i)=0
   16 continue
      call intgrl(a(bfsym),a(orbtbf),a(kadd),a(ladd),a(ijadd),a(ijgrp)
     #,           a(inint),a(inext),a(jmnnxt),a(jmxnxt),a(ningrp)
     #,           a(ijxx),a(klxx),a(nklxx),a(ijww),a(klww),a(nklww),
     $            no34x)
      csav=ningrp+ngroup
      rowsv=csav+nlevs
      wtsav=rowsv+nlevs
      refwt=wtsav+nlevs
      top=refwt+nrefs
      mxref=cortop-refwt+1
      if (top.le.cortop) go to 20
      write (errout,19) top,cortop
   19 format (//,' not quite enough core -- need',i7,' and have'
     #,           i7,' words of core',//)
      call lnkerr('need a smidgeon more core')
   20 continue
cps      call ref(a(arc),a(wght),a(levnr),a(levpt),a(orbtbf),a(bfcode),
cps     #c(bfkey),a(csav),a(rowsv),a(wtsav),a(refwt),mxref)
c
c
      top=refwt+nvref
c
      call print(ops,a(ap),a(b),a(s),a(arc),a(nlwks),a(wght),a(wtab),
     $     a(wtw),a(wtx),a(wty),a(levnr),a(levpt),a(kadd),a(ladd)
     #,          a(inint),a(ijgrp),a(ijadd),a(orbtbf),a(bfcode)
     #,          c(bfkey),a(bfsym))
      call renumb(a(levpt),a(levnr),a(arc))
c
c
c
      call putdrt(a(bfsym),a(bfcode),a(orbtbf),a(iout),a(levpt),a(levnr)
     #,           a(ap),a(b),a(s),a(arc),a(nlwks),a(wght),a(wtab),a(wtw)
     #,           a(wtx),a(wty),a(kadd),a(ladd),a(ijadd),a(ijgrp)
     #,           a(inint),a(inext),a(jmnnxt),a(jmxnxt),a(ningrp)
     #,           a(orbsym),a(ijxx),a(klxx),a(nklxx),a(ijww),a(klww)
     #,           a(nklww),a(refwt))
c
c     write (out,17) top,cortop
   17 format(/,1x,20('*'),' used',i7,' of',i7,' words of core ',22('*'))
c
c     ----- end timing -----
c
c
c     ----- and exit gracefully -----
c
c________________________________________________________________
      call getmem(-cortop,p,idum,'m801',idum)
c________________________________________________________________
      call chainx(0)
      stop
c
c
      end
