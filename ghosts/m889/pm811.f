*deck @(#)pm811.f	1.1  11/30/90
      subroutine pm811(z,a)
c***begin prologue     m811
c***date written       850601   (yymmdd)
c***revision date      900618   (yymmdd)
c
c   18 june   1990    rlm at lanl
c       reading mcscf orbital energies if mcscf-ci run.
c
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
c***keywords           m811, link 811, integrals, transformation,
c***                   one-electron, two-electron
c***author             saxe, paul (lanl)
c***source             @(#)pm811.f	1.1   11/30/90
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
c***end prologue       m811
c
      implicit integer (a-z)
c
      parameter (maxnbf=300)
c
      character*4096 ops
      character chrkey*32,xform*32
      character*32 occ
      character*128 namint
      character*16 bflabl(maxnbf)
      logical logkey
      logical toguga
      real*8 z,enuc,ecore,sdot
      integer a(*)
      dimension z(*)
c
c      common //       z(20000)
c
      common /io/     inp,ioutpt
c
c      equivalence (z,a)
c
      data maxcor /20000/
c
c
 1000 format(1x,'transformation:')
 1010 format(5x,'orbital set:',a16,
     $     /,5x,'number of symmetries',14x,i4,
     $     /,5x,'orbitals transformed',14x,i4,
     $     /,5x,'orbitals frozen     ',14x,i4)
 1020 format(5x,'frozen core energy',5x,f15.9)
 1030 format(5x,'memory use        ',11x,i9)
cdir$ fastmd
c
c     ----- open read-write and checkpoint files -----
c
c..bhl      call drum
c
c     ----- recover the options string -----
c
      ops=' '
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     ----- start timing routines ------
c
c
c     process the other options.
c
      toguga=.not.logkey(ops,'ci=full',.false.,' ')
      xform=chrkey(ops,'ci=tvector',' ',' ')
      if(xform.eq.' ') then
         call iosys('read character "transformation vector" from rwf',
     $        -1,0,0,xform)
      endif
      call iosys('write character "transformation vector" to rwf',
     $            -1,0,0,xform)
c
c     ----- work out the name of the corresponding orbital eigenvalues
c           or occupations
c
      call locase(xform,xform)
      if (xform(1:9).eq.'no vector') then
         occ='"no occ'//xform(10:31)//'"'
      else if (xform(1:12).eq.'mcscf vector') then
         occ='"mcscf orbital energies"'
      else
         occ='"orbital energies"'
      end if
c
c     ----- read dimensions etc from wherever they are -----
c
      call iosys('read integer "number of basis functions" from rwf',
     $     1,nbf,0,' ')
c
      write(ioutpt,1000)
      if (logkey(ops,'transformation=scf',.false.,' ').or.
     #    logkey(ops,'kohn',.false.,' ')) then
         toguga=.false.
         nsym=1
         norbs=nbf
         nrows=1
         nlevs=1
         nrefs=1
         levfrm=1
         symorb=1
         ngroup=1
         nmax=1
         numij=norbs*(norbs+1)/2
c
         kadd=27
         ladd=kadd+symorb
         ijadd=ladd+symorb
         ijgrp=ijadd+numij
         iout=ijgrp+numij
         orbsym=iout+nbf
         need=orbsym+norbs
c
         write (ioutpt,1050) xform,nbf
 1050    format (5x,'vector used:               ',a16,/,
     $           5x,'number of basis functions: ',i4)
      else
         call iosys('read integer "symmetries in ci" from rwf',
     $        1,nsym,0,' ')
         call iosys('read integer norbs from rwf',1,norbs,0,' ')
         call iosys('read integer nrows from rwf',1,nrows,0,' ')
         call iosys('read integer nlevs from rwf',1,nlevs,0,' ')
         call iosys('read integer nrefs from rwf',1,nrefs,0,' ')
         call iosys('read integer orbfrm from rwf',1,levfrm,0,' ')
         levfrm=levfrm+1
         call iosys('read integer symorb from rwf',1,symorb,0,' ')
         call iosys('read integer ngroup from rwf',1,ngroup,0,' ')
         call iosys('read integer nmax from rwf',1,nmax,0,' ')
c
         write(ioutpt,1010) xform,nsym,norbs,nbf-norbs
c
         numij=(norbs+1)*norbs/2
c
c        ----- divide core and get the drt arrays -----
c
         kadd=27
         ladd=kadd+symorb
         ijadd=ladd+symorb
         ijgrp=ijadd+numij
         iout=ijgrp+numij
         orbsym=iout+nbf
         need=orbsym+norbs
c
         call getdrt(a(kadd),a(ladd),a(ijadd),a(ijgrp),
     #        a(iout),a(orbsym),end,ngroup,nrefs,symorb,numij,nbf,
     #        norbs,levfrm,nlevs,nrows,nsym)
      end if
c
c     ----- now get the transformation vector -----
c
      hmo=iadtwp(need)
      c=hmo+numij
      tempc=c+nbf*norbs
      eigval=tempc+nbf**2
      e=eigval+nbf
      tempe=e+norbs
      need=wpadti(tempe+nbf)
      call getscm(need,z,ngot,'m811 xform vector',0)
c
      call iosys('read real '//xform//' from rwf',-1,z(tempc),0,' ')
      if(.not.logkey(ops,'mcscf',.false.,' ')) then
       call iosys('read real '//occ//' from rwf',nbf,z(tempe),0,' ')
      end if
c
c     ----- reorder the vector to the drt order and form a core
c           density matrix, if necessary
c
      nnp=(nbf+1)*nbf/2
      cdens=tempe+nbf
      need=wpadti(cdens+nnp)
      call getscm(need,z,ngot,'m811 core fock',0)
c
      if (nbf.gt.maxnbf) call lnkerr('length of bflabl in m811 '//
     $     'needs to be increased')
c
      call reordr(z(tempc),z(c),z(cdens),a(iout),nbf,norbs,nnp,
     $     ncore,ops,bflabl,z(tempe),z(e))
c
c     ----- open the integral file and check the orthonormality of
c                   the orbitals
c
      call iosys('read character "integral filename" from rwf',
     $     0,0,0,namint)
      call iosys('open ints as old',0,0,0,namint)
c
      s=iadtwp(need)
      h=s+nnp
      v=h+nnp
      need=wpadti(v+nnp)
      call getscm(need,z,ngot,'m811 one-electron',0)
c
      call iosys('read real "overlap integrals" from rwf',
     $            -1,z(s),0,' ')
      call iosys('read real "kinetic integrals" from rwf',
     $            -1,z(h),0,' ')
      call iosys('read real "potential integrals" from rwf',
     $     -1,z(v),0,' ')
c
      call vadd(z(h),z(h),z(v),nnp)
c
      t1=iadtwp(need)
      t2=t1+nbf**2
      t3=t2+nbf**2
      need=wpadti(t3+nbf**2)
      call getscm(need,z,ngot,'m811 chknrm',0)
c
      call chknrm(z(tempc),z(s),z(t1),z(t2),z(t3),nbf,nnp)
c
c     ----- form the frozen-core fock-matrix -----
c
      ecore=0.0d+00
      if (ncore.gt.0) then
         cf=iadtwp(need)
         values=cf+nnp
         ntriang=min(nnp,30000/nnp+1)
         need=wpadti(values+nnp*ntriang)
         call getscm(need,z,ngot,'m811 fock',0)
c
         call fock(z(values),z(cdens),z(cf),nnp,nbf,z(t1),z(t2),
     #             ntriang)
c
c     ----- form the frozen core energy -----
c
         call vadd(z(cf),z(cf),z(h),nnp)
         call trtosq(z(t1),z(cdens),nbf,nnp)
         call trtosq(z(t2),z(cf),nbf,nnp)
         call trtosq(z(t3),z(h),nbf,nnp)
         ecore=0.5d+00*(sdot(nbf**2,z(t1),1,z(t2),1)+
     $        sdot(nbf**2,z(t1),1,z(t3),1))
         write(ioutpt,1020) ecore
         call vmove(z(h),z(cf),nnp)
      end if
c
c     ----- store the frozen core energy -----
c
      call iosys('write real "frozen core energy" to rwf',
     $            1,ecore,0,' ')
      call iosys('write real "frozen core energy" to ints',
     $            1,ecore,0,' ')
c
c     ----- get nuclear repulsion energy -----
c
      call iosys('read real "nuclear repulsion energy" from rwf',
     $     1,enuc,0,' ')
c
c     ----- transform the one-electron integrals -----
c
      call trn1e(z(h),z(c),nbf,norbs,nnp,numij,z(t1),z(t2),z(hmo))
c
c     ----- reallocate core for the two-electron transformation -----
c
      t1=c+nbf*norbs
      t2=t1+nbf**2
      val=t2+norbs*nbf
      lenbin=max(1920,nnp)
      lab=wpadti(val+lenbin)
      bin=lab+lenbin
      values=iadtwp(bin+lenbin)
      ntriang=min(nnp,30000/nnp+1)
      asort=values+nnp*ntriang
      call getscm(0,z,maxcor,'possible',0)
c..bhl.cos.io
      maxcor=min(maxcor,2000000)
c..bhl.cos.io
c
      lnsort=min(wptoin(max(nnp**2/2,nmax*ngroup)+1000),
     $     maxcor-wpadti(asort)-100)
      need=wpadti(asort)+lnsort
      call getscm(need+100,z,ngot,'m811 two-electron',0)
c
      write(ioutpt,1030) need
      if(ntriang.le.0) call lnkerr('not enough memory')
c
      continue
c
      call trn2e(z(values),nnp,ntriang,z(c),nbf,norbs,z(t1),z(t2),
     #           numij,z(asort),lnsort,ngroup,
     #           nmax,nsym,a(ijgrp),a(ijadd),a(kadd),a(ladd),
     #           a(orbsym),z(hmo),levfrm,
     #           z(val),a(lab),a(bin),lenbin,toguga,z(values),ops)
c
c     ----- stop the timing routines -----
c
c
c     ----- and exit with grace -----
c
      call chainx(0)
c
c
      stop
      end
