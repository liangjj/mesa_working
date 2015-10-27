*deck @(#)m822.f	1.2  7/30/91
      program m822
c***begin prologue     pm822
c***date written       850601   (yymmdd)
c***revision date      910725   (yymmdd)
c
c  25 july    1991    rlm at lanl
c       the incoming transformed integrals may be picked up from
c       either 'zints' or 'tints' depending upon the circumstance.
c  18 march   1988    bhl at llnl
c       max. maxcor set to 2,000,000 for cos io
c
c  10 january 1988    bhl at brl
c       skipping read of orbital energies if mcscf-ci run
c
c  13 october 1987    pws at lanl
c       adding options to print out the vector before and after
c       reordering to the drt order.
c
c   3 december 1986   pws at lanl
c       changing 'namint' and iosys open to character
c
c***keywords           m822, link 822, integrals, transformation,
c***                   one-electron, two-electron
c***author             saxe, paul (lanl)
c***source             @(#)pm822.f	1.2   7/30/91
c***purpose            transforms the one- and two-electron integrals
c                      from the atomic-orbital to the molecular orbital
c                      basis.
c***description
c     the option substrings currently recognized are:
c         timing                collect and report timing statistics
c         xform=string          transform the orbitals given by string.
c                               at present this may be scf_vector or
c                               "no vector nnn", where nnn denotes the root
c                               number.
c
c
c***references
c
c***routines called
c***end prologue       pm822
c
      implicit integer (a-z)
c
c
      character*4096 ops
      character*8 unit
      character*128 tints,zints,gints,nmcnfg,namchk
      character*16 key, chrkey
      character*8 dsk
      character*3 ans
      logical toguga
      logical logkey
      integer a
      real*8 z
      data maxtri/1000000/
      pointer(p,z(1)),(p,a(1))
c
c
      common /io/inp,ioutpt
c
c
 1000 format(1x,'m822: sort to guga order',/)
 1010 format(5x,'orbital set:',a16,
     $     /,5x,'number of symmetries',14x,i4,
     $     /,5x,'orbitals transformed',14x,i4,
     $     /,5x,'orbitals frozen     ',14x,i4)
 1030 format(5x,'memory use        ',11x,i9)
c
c     ----- recover the options string -----
c
      call drum
      call iosys('read character options from rwf',-1,0,0,ops)
      call iosys('read character "checkpoint filename" from rwf',
     $            0,0,0,namchk)
      call iosys('open chk as old',0,0,0,namchk)
      key=chrkey(ops,'int=drt=key','drt',' ')
      call pakstr(key,lenkey)
      call iosys('does "drt file name '//key(1:lenkey)//
     1           '" exist on chk',0,0,0,ans)
      if(ans.ne.'yes') then
         write(ioutpt,*) '         drt file does not exist. quit'
         call lnkerr('quit m822')
      endif
      call iosys('read character "drt file name '//key(1:lenkey)
     1           //'" from chk',0,0,0,dsk)
      if(dsk.eq.'hconfig') then
         write(ioutpt,*) '            reading drt information '// 
     #                   'from hconfig'
         call iosys('read character "hamiltonian manipulation '//
     #               'filename" from rwf',0,0,0,nmcnfg)
         call iosys('open hconfig as unknown',0,0,0,nmcnfg)
      else
         call iosys('read character "drt file name '//key(1:lenkey)
     1               //'" from chk',0,0,0,dsk)
         call iosys('read character "drt unit file name '//dsk
     #               //'" from chk',0,0,0,nmcnfg)
         write(ioutpt,*) '            reading drt information '// 
     #                   'from '//dsk
         call iosys('open '//dsk//' as unknown',0,0,0,nmcnfg)
      call iosys('close chk',0,0,0,namchk)
      endif
c
c     ----- process the other options. -----
c
      toguga=.true.
c
c     ----- read dimensions etc from wherever they are -----
c
      call iosys('read integer "number of basis functions" from rwf',
     $            1,nbf,0,' ')
c
      write(ioutpt,1000)
      call iosys('read integer "symmetries in ci" from '//dsk,
     $            1,nsym,0,' ')
      call iosys('read integer norbs from '//dsk,1,norbs,0,' ')
      call iosys('read integer nrows from '//dsk,1,nrows,0,' ')
      call iosys('read integer nlevs from '//dsk,1,nlevs,0,' ')
      call iosys('read integer nrefs from '//dsk,1,nrefs,0,' ')
      call iosys('read integer orbfrm from '//dsk,1,levfrm,0,' ')
      levfrm=levfrm+1
      call iosys('read integer symorb from '//dsk,1,symorb,0,' ')
      call iosys('read integer ngroup from '//dsk,1,ngroup,0,' ')
      call iosys('read integer nmax from '//dsk,1,nmax,0,' ')
c
      call getmem(0,p,ngot,'first',0)
      call iosys('read integer mxcore from rwf',1,maxcor,0,' ')
c
c
      numij=(norbs+1)*norbs/2
      nnp=nbf*(nbf+1)/2
c
c     ----- divide core and get the drt arrays -----
c
      kadd=27
      ladd=kadd+symorb
      ijadd=ladd+symorb
      ijgrp=ijadd+numij
      iout=ijgrp+numij
      orbsym=iout+nbf
      hmo=iadtwp(orbsym+norbs)
      c=hmo+numij
      t1=c+nbf*norbs
      t2=t1+nbf**2
      val=t2+norbs*nbf
      lenbin=max(1920,nnp)
      lab=wpadti(val+lenbin)
      bin=lab+lenbin
      values=iadtwp(bin+lenbin)
      ntriang=min(nnp,maxtri/nnp+1)
      asort=values+nnp*ntriang
c
      lnsort=min(wptoin(max(nnp**2/2,nmax*ngroup)+5000),
     $           maxcor-wpadti(asort)-100)
      need=wpadti(asort)+lnsort + 100
      call getmem(need,p,ngot,'m822',0)
c
      write(ioutpt,1030) need
      if(ntriang.le.0) then
         call lnkerr('not enough memory')
      endif
c      
c
      call getdrt(a(kadd),a(ladd),a(ijadd),a(ijgrp),
     #            a(iout),a(orbsym),end,ngroup,nrefs,symorb,numij,
     $            nbf,norbs,levfrm,nlevs,nrows,nsym,dsk)
c
c     ----- now sort the integrals
c
      

c     ----- open the integral units -----
c
      if (logkey(ops,'int=zint',.false.,' ')) then
         unit='zints'
         call iosys('read character "zeroed integral filename"'
     $            //' from rwf',0,0,0,zints)
         call iosys('open zints as old',0,0,0,zints)
      else 
c        default to tint
         unit='tints'
         call iosys('read character "transformed integral filename"'
     $            //' from rwf',0,0,0,tints)
         call iosys('open tints as old',0,0,0,tints)
      end if
c
      len=3*(numij*numij+numij)/2
      call iosys('read character "guga integral filename"'
     $         //' from rwf',0,0,0,gints)
      call iosys('open gints as new',len,0,0,gints)
c
      call srt2e(z(values),nnp,ntriang,z(c),nbf,norbs,z(t1),z(t2),
     #           numij,z(asort),lnsort,ngroup,
     #           nmax,nsym,a(ijgrp),a(ijadd),a(kadd),a(ladd),
     #           a(iout),z(hmo),levfrm,
     #           z(val),a(lab),a(bin),lenbin,toguga,z(values),
     $           unit,ops)
      call trn2e(z(values),nnp,ntriang,z(c),nbf,norbs,z(t1),z(t2),
     #           numij,z(asort),lnsort,ngroup,
     #           nmax,nsym,a(ijgrp),a(ijadd),a(kadd),a(ladd),
     #           a(orbsym),z(hmo),levfrm,
     #           z(val),a(lab),a(bin),lenbin,toguga,z(values),ops)
      call getmem(-ngot,p,idum,'m822',idum)
c
c     ----- and exit with grace -----
c
      call iosys('close '//unit,0,0,0,' ')
      call iosys('close gints',0,0,0,' ')
      call chainx(0)
c
c
      stop
      end
