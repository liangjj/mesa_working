*deck @(#)pm811.f	5.1 11/6/94
      subroutine pm811(z,a)
c***begin prologue     m811
c***date written       850601   (yymmdd)
c***revision date      910225   (yymmdd)
c
c  25 february 1991   rlm at lanl
c       deleting the limitation on maxcor implemented 18 march 1988.
c       this causes problems on workstations, and is machine dependent
c       anyway.
c   2 january 1991  mob at lanl
c       adding option to transform 1e spin-orbit integrals
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
c***source             @(#)pm811.f	5.1 11/6/94 
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
      parameter (maxnbf=2000)
c
      character*4096 ops
      character chrkey*32,xform*32
      character*32 occ
      character*128 namint,tints,gints
      character*16 bflabl(maxnbf)
      logical logkey
      logical toguga
      logical tocan
      logical d2e
      real*8 enuc,ecore,sdot,schmrr,fpkey
      integer a(*)
      real*8 z(*)
c
      common /io/     inp,ioutpt
c
      data tocan/.false./
      save tocan
c
 1000 format(1x,'m811:transformation')
 1010 format(5x,'orbital set:',a32,
     $     /,5x,'number of symmetries',14x,i4,
     $     /,5x,'orbitals transformed',14x,i4,
     $     /,5x,'orbitals frozen     ',14x,i4)
 1015 format(5x,'orbital set:',a32,
     $     /,5x,'number of basis functions:',8x,i4)
 1020 format(5x,'guga order')
 1025 format(5x,'canonical order')
 1030 format(5x,'frozen core energy',5x,f15.9)
 1040 format(5x,'memory use        ',11x,i9)
cdir$ fastmd
c
c     ----- recover the options string -----
c
      ops=' '
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     --- this link usually produces transformed integrals on 'gints' which are 
c     in guga order so they can be read by m901/m902.  in the case that
c     the ci to be used is m901, m811 should be followed by m821
c     to sort the guga integrals into supermatrix format.
c
c     if one is doing full ci(m903), the standard canonical transformed integral
c     set is needed, and so toguga is turned off and tocan turned on.
c     this is also the case if one is doing some form of ci optimization.
c     the canonical transformed integrals are need by the coupled perturbed
c     equation links.  this link produces canonical integrals, and is
c     followed by m820 which sorts to guga order, and then m821 which forms
c     supermatrices to be read by m901 or m902.
c
c
c     --- process other options
      xform=chrkey(ops,'ci=tvector',' ',' ')
      schmrr=fpkey(ops,'schmidt-tol',5.d-08,' ')
      if(xform.eq.' ') then
         call iosys('read character "transformation vector" from rwf',
     $        -1,0,0,xform)
      endif
      call iosys('write character "transformation vector" to rwf',
     $            0,0,0,xform)
c
c     --- work out the name of the corresponding orbital eigenvalues
c           or occupations
      call locase(xform,xform)
      if (xform(1:9).eq.'no vector') then
         occ='"no occ'//xform(10:31)//'"'
      else if (xform(1:12).eq.'mcscf vector') then
         occ='"mcscf orbital energies"'
      else
         occ='"orbital energies"'
      end if
c
c     --- read dimensions etc from wherever they are -----
      call iosys('read integer "number of basis functions" from rwf',
     $     1,nbf,0,' ')
c
c     --- determine whether to do canonical or guga transformation.
      write(ioutpt,1000)
      toguga=.true.
      if   ((logkey(ops,'ci=full',.false.,' '))
     $  .or.(logkey(ops,'ciopt',.false.,' '))
     $  .or.(logkey(ops,'opt=ci',.false.,' '))
     $  .or.(logkey(ops,'mrciopt',.false.,' '))
     $  .or.(logkey(ops,'opt=mrci',.false.,' '))) then
         toguga=.false.
      endif
      tocan=.not.toguga
c     --- may want to do canonical because we're doing analytical force 
c         constants
      d2e=.false.
      if(logkey(ops,'force-constants',.false.,' ')) then
         if(logkey(ops,'force-constants=numerical',.false.,' ')) then
c           numerical force constants, do nothing
         else
c           analytic force constants
            d2e=.true.
         endif 
      endif
      if(logkey(ops,'transformation=canonical',.false.,' ').or.
     $    d2e.or.
     #    logkey(ops,'kohn',.false.,' ').or.
     #    logkey(ops,'m819',.false.,' ')) then
         toguga=.false.
         tocan=.true.
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
         write (ioutpt,1025)
         write (ioutpt,1015) xform,nbf
c
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
         if(toguga) then
            write(ioutpt,1020)
         else
            write(ioutpt,1025)
         endif
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
c
      call reordr(z(tempc),z(c),z(cdens),a(iout),nbf,norbs,nnp,
     $     ncore,ops,bflabl,z(tempe),z(e),tocan)
c
c     ----- open the integral files and check the orthonormality of
c                   the orbitals
c
      call iosys('read character "integral filename" from rwf',
     $     0,0,0,namint)
      call iosys('open ints as old',0,0,0,namint)
c
      if(tocan) then
c        normal canonical transformation
         len=numij*numij+numij
         call iosys('read character "transformed integral filename"'
     $            //' from rwf',0,0,0,tints)
         if(logkey(ops,'unit=ssd=tint',.false.,' ')) then
            call iosys('open tints as new on ssd',len,0,0,tints)
         else
            call iosys('open tints as new',len,0,0,tints)
         endif
      endif
      if(toguga) then
         len=3*(numij*numij+numij)/2
         call iosys('read character "guga integral filename"'
     $            //' from rwf',0,0,0,gints)
         if(logkey(ops,'unit=ssd=gint',.false.,' ')) then
            call iosys('open gints as new on ssd',len,0,0,gints)
         else
            call iosys('open gints as new',len,0,0,gints)
         endif
      endif
c
      s=iadtwp(need)
      h=s+nnp
      v=h+nnp
      need=wpadti(v+nnp)
      call getscm(need,z,ngot,'m811 one-electron',0)
c
      call iosys('read real "overlap integrals" from rwf',
     $            nnp,z(s),0,' ')
      call iosys('read real "kinetic integrals" from rwf',
     $            nnp,z(h),0,' ')
      call iosys('read real "potential integrals" from rwf',
     $            nnp,z(v),0,' ')
c
      call vadd(z(h),z(h),z(v),nnp)
c
      t1=iadtwp(need)
      t2=t1+nbf**2
      t3=t2+nbf**2
      need=wpadti(t3+nbf**2)
      call getscm(need,z,ngot,'m811 chknrm',0)
c
      call chknrm(z(tempc),z(s),z(t1),z(t2),z(t3),schmrr,nbf,nnp)
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
         write(ioutpt,1030) ecore
         call scopy(nnp,z(cf),1,z(h),1)
      end if
c
c     ----- store the frozen core energy -----
c
      call iosys('write real "frozen core energy" to rwf',
     $            1,ecore,0,' ')
c
c     ----- get nuclear repulsion energy -----
c
      call iosys('read real "nuclear repulsion energy" from rwf',
     $     1,enuc,0,' ')
c
c     ----- transform the one-electron integrals -----
      call trn1e(z(h),z(c),nbf,norbs,nnp,numij,z(t1),z(t2),z(hmo))
c
c     ---- perhaps transform one-electron spin-orbit integrals ----
      if(logkey(ops,'spin-orbit',.false.,' ')) then
         call trn1so(z(c),nbf,norbs,numij,z(t1),z(t2),z(t3),bflabl)
      end if
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
c     maxcor=min(maxcor,2000000)
c..bhl.cos.io
c
      lnsort=min(wptoin(max(nnp**2/2,nmax*ngroup)+5000),
     $     maxcor-wpadti(asort)-100)
      need=wpadti(asort)+lnsort
      call getscm(need+100,z,ngot,'m811 two-electron',0)
c
      write(ioutpt,1040) need
      if(ntriang.le.0) call lnkerr('not enough memory')
c
      continue
c
      call trn2e(z(values),nnp,ntriang,z(c),nbf,norbs,z(t1),z(t2),
     #           numij,z(asort),lnsort,ngroup,
     #           nmax,nsym,a(ijgrp),a(ijadd),a(kadd),a(ladd),
     #           a(orbsym),z(hmo),levfrm,
     #           z(val),a(lab),a(bin),lenbin,toguga,z(values),
     $           tocan,ops)
c
c     ----- and exit with grace -----
c
      if(tocan) then
         call iosys('close tints',0,0,0,' ')
      endif
      if(toguga) then
         call iosys('close gints',0,0,0,' ')
      endif
      call chainx(0)
c
c
      stop
      end
