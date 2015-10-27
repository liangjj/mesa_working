*deck @(#)pm811.f	1.3  7/30/91
      program m811
c***begin prologue     m811
c***date written       850601   (yymmdd)
c***revision date      910225   (yymmdd)
c
c  25 february 1991   rlm at lanl
c       deleting the limitation on maxcor implemented 18 march 1988.
c       this causes problems on workstations, and is machine dependent anyway.
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
c***source             @(#)pm811.f	1.3   7/30/91
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
c***end prologue       pm811
c
      implicit integer (a-z)
c
      parameter (maxnbf=2000)
c
      character*4096 ops
      character chrkey*32,xform*32
      character*3 ans
      character*32 occ
      character*16 key
      character*5 dsk
      character*8 scatyp
      character*128 namint,tints,gints
      character*16 bflabl(maxnbf)
      character*128 nmpart
      logical scat
      logical logkey
      logical toguga
      logical tocan
      real*8 enuc,ecore,sdot,schmrr,fpkey
      integer need(6), ngot(6)
      integer a1, a2
      real*8 z1, z2, z3, z4
      pointer(pa1,a1(1))
      pointer(pz1,z1(1))
      pointer(pz2,z2(1)),(pz2,a2(1))
      pointer(pz3,z3(1))
      pointer(pz4,z4(1))
c
      common /io/     inp,ioutpt
c
      data tocan/.false./
      data maxtri/1000000/
c
 1000 format(1x,'m811:transformation',//)
 1010 format(5x,'orbital set:',a16,
     $     /,5x,'number of symmetries',14x,i4,
     $     /,5x,'orbitals transformed',14x,i4,
     $     /,5x,'orbitals frozen     ',14x,i4)
 1015 format(5x,'orbital set:',a16,
     $     /,5x,'number of basis functions:',8x,i4)
 1020 format(5x,'guga order')
 1025 format(5x,'canonical order')
 1030 format(5x,'frozen core energy',5x,f15.9)
 1040 format(5x,'memory use        ',11x,i9)
cdir$ fastmd
c
      call drum
c     ----- recover the options string -----
c
      write(ioutpt,1000)      
      ops=' '
      call iosys('read character options from rwf',-1,0,0,ops)
      scatyp=chrkey(ops,'scattering','none',' ')
      scat=.false.
      if(scatyp.eq.'kohn'.or.scatyp.eq.'r-matrix') then
         scat=.true.
      end if
      if(logkey(ops,'int=reuse=m811',.false.,' ') ) then
         write(ioutpt,*) 'skipping m811 transformation'
      else                                               
         key=chrkey(ops,'int=drt=key','drt',' ')
         call pakstr(key,lenkey)
         call iosys('does "drt file name '//key(1:lenkey)//
     1              '" exist on rwf',0,0,0,ans)
         if(ans.ne.'yes') then
            write(ioutpt,*) '         drt file does not exist. quit'
            call lnkerr('quit m811')
         endif
         call iosys('read character "drt file name '//key(1:lenkey)
     1               //'" from rwf',0,0,0,dsk)
         write(ioutpt,*) '            reading drt information '//
     $                   'from '//dsk
         if(dsk.eq.'hpart') then
            call iosys('read character "partitioning filename" '//
     $                 'from rwf',0,0,0,nmpart)
            call iosys('open hpart as unknown',0,0,0,nmpart)
         endif
c
c        process the other options.
c
c        this link usually produces transformed integrals on 'gints' 
c        which are in guga order so they can be read by m901/m902.  
c        in the case that the ci to be used is m901, m811 should be 
c        followed by m821 to sort the guga integrals into 
c        supermatrix format.
c
c        if one is doing full ci(m903), the standard canonical 
c        transformed integral set is needed, and so toguga is turned 
c        off and tocan turned on. this is also the case if one is doing 
c        some form of ci optimization.
c        the canonical transformed integrals are need by the coupled 
c        perturbed equation links.  this link produces canonical 
c        integrals, and is followed by m820 which sorts to guga order, 
c        and then m821 which forms supermatrices to be read by m901 or m902.
c
         toguga=.true.
         if   ((logkey(ops,'ci=full',.false.,' '))
     $     .or.(logkey(ops,'ciopt',.false.,' '))
     $     .or.(logkey(ops,'opt=ci',.false.,' '))
     $     .or.(logkey(ops,'mrciopt',.false.,' '))
     $     .or.(logkey(ops,'opt=mrci',.false.,' '))) then
            toguga=.false.
         endif
         tocan=.not.toguga
c
         xform=chrkey(ops,'ci=tvector',' ',' ')
         schmrr=fpkey(ops,'schmidt-tol',5.d-08,' ')
         lenbin=intkey(ops,'m811=sort-bin-size',2000000,' ')
         if(xform.eq.' ') then
            call iosys('read character "transformation vector" '//
     $                 'from rwf',-1,0,0,xform)
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
     $        1,nbf,0,' ')
c
         if (logkey(ops,'transformation=scf',.false.,' ')
     #                     .or.
     #                     scat
     $                     .or.
     #       logkey(ops,'m819',.false.,' ')) then
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
c                            a1(kadd)
            kadd=27
c                            a1(ladd)
            ladd=kadd+symorb
c                            a1(ijadd)
            ijadd=ladd+symorb
c                            a1(ijgrp)
            ijgrp=ijadd+numij
c                            a1(iout)
            iout=ijgrp+numij
c                            a1(orbsym)
            orbsym=iout+nbf
c
c           get the words
c
            need(1)=orbsym+norbs
            call getmem(need(1),pa1,ngot(1),'m811',0)
            write (ioutpt,1025)
            write (ioutpt,1015) xform,nbf
c
         else
            call iosys('read integer "symmetries in ci" from '//dsk,
     $           1,nsym,0,' ')
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
c                             a1(kadd)
            kadd=27
c                             a1(ladd)
            ladd=kadd+symorb
c                             a1(ijadd)
            ijadd=ladd+symorb
c                             a1(ijgrp)
            ijgrp=ijadd+numij
c                             a1(iout)
            iout=ijgrp+numij
c                             a1(orbsym)
            orbsym=iout+nbf
c   
c           get the words
c
            need(1)=orbsym+norbs
            call getmem(need(1),pa1,ngot(1),'m811',0)
c
            call getdrt(a1(kadd),a1(ladd),a1(ijadd),a1(ijgrp),
     #           a1(iout),a1(orbsym),end,ngroup,nrefs,symorb,numij,
     #           nbf,norbs,levfrm,nlevs,nrows,nsym,dsk)
         end if
c
c     ----- now get the transformation vector -----
c
c
c
c     allocate z1  this is done separately because they are needed in the
c                  two electron step
c
c                     z1(hmo)
         hmo=1
c                     z1(c)
         c=hmo+numij
c
c        get the words for z1
c
         need(2)=wpadti(c+nbf*norbs)
         call getmem(need(2),pz1,ngot(2),'m811',0)
c
c        allocate z2
c   
c                               z2(tempc) 
         tempc=1
c                               z2(eigval)
         eigval=tempc+nbf**2
c                               z2(e)
         e=eigval+nbf
c                               z2(tempe)
         tempe=e+norbs
c                 
         nnp=(nbf+1)*nbf/2
c                               z2(cdens)
         cdens=tempe+nbf
c
c        get the words for z2
c
         need(3)=wpadti(cdens+nnp)
c      call getscm(need,z,ngot,'m811 xform vector',0)
         call getmem(need(3),pz2,ngot(3),'m811',0)
c
         call iosys('read real '//xform//' from rwf',-1,z2(tempc),0,' ')
         if(.not.logkey(ops,'mcscf',.false.,' ')) then
            call iosys('read real '//occ//' from rwf',nbf,z2(tempe),
     #                  0,' ')
         end if
c
c     ----- reorder the vector to the drt order and form a core
c           density matrix, if necessary
c
c      call getscm(need,z,ngot,'m811 core fock',0)
c
c
         call reordr(z2(tempc),z1(c),z2(cdens),a1(iout),nbf,norbs,nnp,
     $        ncore,ops,bflabl,z2(tempe),z2(e),tocan)
c
c     ----- open the integral files and check the orthonormality of
c                   the orbitals
c
         call iosys('read character "integral filename" from rwf',
     $               0,0,0,namint)
         call iosys('open ints as old',0,0,0,namint)
c
         if(tocan) then
c           normal canonical transformation
            len=numij*numij+numij
            call iosys('read character "transformed integral filename"'
     $               //' from rwf',0,0,0,tints)
            if(logkey(ops,'unit=ssd=tint',.false.,' ')) then
               call iosys('open tints as new on ssd',len,0,0,tints)
            else
               call iosys('open tints as new',len,0,0,tints)
            endif
         endif
         if(toguga) then
            len=3*(numij*numij+numij)/2
            call iosys('read character "guga integral filename"'
     $                 //' from rwf',0,0,0,gints)
            if(logkey(ops,'unit=ssd=gint',.false.,' ')) then
               call iosys('open gints as new on ssd',len,0,0,gints)
            else
               call iosys('open gints as new',len,0,0,gints)
            endif
         endif
c
c         s=iadtwp(need)
c
c        allocate z3 for the one-electron transformation
c
c                      z3(s)
         s=1
c                      z3(h)
         h=s+nnp
c                      z3(v)      
         v=h+nnp
c                      z3(t1)      
         t1=v+nnp
c                      z3(t2)      
         t2=t1+nbf**2
c                      z3(t3)         
         t3=t2+nbf**2
         need(4)=wpadti(t3+nbf**2)
         call getmem(need(4),pz3,ngot(4),'m811',0)
c      call getscm(need,z,ngot,'m811 one-electron',0)
c
         call iosys('read real "overlap integrals" from rwf',
     $               nnp,z3(s),0,' ')
         call iosys('read real "kinetic integrals" from rwf',
     $               nnp,z3(h),0,' ')
         call iosys('read real "potential integrals" from rwf',
     $               nnp,z3(v),0,' ')
c
         call vadd(z3(h),z3(h),z3(v),nnp)
c
c     call getscm(need,z,ngot,'m811 chknrm',0)
c
         call chknrm(z2(tempc),z3(s),z3(t1),z3(t2),z3(t3),schmrr,
     $               nbf,nnp)
c
c     ----- form the frozen-core fock-matrix -----
c
         ecore=0.0d+00
         if (ncore.gt.0) then
c            cf=iadtwp(need)
c
c           allocate z4
c
c                       z4(cf)
            cf=1
c                       z4(values)         
            values=cf+nnp
            ntriang=min(nnp,maxtri/nnp+1)
            need(5)=wpadti(values+nnp*ntriang)
c         call getscm(need,z,ngot,'m811 fock',0)
            call getmem(need(5),pz4,ngot(5),'m811',0)
c
            call fock(z4(values),z2(cdens),z4(cf),nnp,nbf,
     #                z3(t1),z3(t2),ntriang)
c
c     ----- form the frozen core energy -----
c
            call vadd(z4(cf),z4(cf),z3(h),nnp)
            call trtosq(z3(t1),z2(cdens),nbf,nnp)
            call trtosq(z3(t2),z4(cf),nbf,nnp)
            call trtosq(z3(t3),z3(h),nbf,nnp)
            ecore=0.5d+00*(sdot(nbf**2,z3(t1),1,z3(t2),1)+
     $           sdot(nbf**2,z3(t1),1,z3(t3),1))
            write(ioutpt,1030) ecore
            call scopy(nnp,z4(cf),1,z3(h),1)
c
c           free z4
c
            call getmem(-ngot(5),pz4,idum,'m811',idum)
         end if
c
c     ----- store the frozen core energy -----
c
         call iosys('write real "frozen core energy" to rwf',
     $               1,ecore,0,' ')
c
c     ----- get nuclear repulsion energy -----
c
         call iosys('read real "nuclear repulsion energy" from rwf',
     $               1,enuc,0,' ')
c
c     ----- transform the one-electron integrals -----
c
         call trn1e(z3(h),z1(c),nbf,norbs,nnp,numij,
     #              z3(t1),z3(t2),z1(hmo))
c
c     ----- reallocate core for the two-electron transformation -----
c
c
c        free z3 and z2
c
         call getmem(-ngot(4),pz3,idum,'m811',idum)
      call getmem(-ngot(3),pz2,idum,'m811',idum)
c         t1=c+nbf*norbs
c                            z2(t1)
         t1=1
c                            z2(t2)      
         t2=t1+nbf**2
c                            z2(val)      
         val=t2+norbs*nbf
c     ----- lenbin now input parameter -----
c                            a2(lab)
         lab=wpadti(val+lenbin)
c                            a2(bin)        
         bin=lab+lenbin
c                            z2(values)      
         values=iadtwp(bin+lenbin)
c     ----- ntriang now based on value of lenbin -----
         ntriang=min(nnp,lenbin/nnp+1)
c                            z2(asort)       
         asort=values+nnp*ntriang
c
c     find out how much memory we can use for the sort
c      
         call iosys('read integer mxcore from rwf',1,maxcor,0,' ')
c      
c      call getscm(0,z,maxcor,'possible',0)
c..bhl.cos.io
c     maxcor=min(maxcor,2000000)
c..bhl.cos.io
c
         lnsort=min(wptoin(max(nnp**2/2,nmax*ngroup)+5000),
     $        maxcor-wpadti(asort)-100)
         write(ioutpt,*) '    no. triangles/read = ',ntriang
         write(ioutpt,*) '    length of sort bin = ',lenbin
         write(ioutpt,*) '    memory used for sort = ',lnsort
         need(6)=wpadti(asort)+lnsort + 100
c      call getscm(need+100,z,ngot,'m811 two-electron',0)
         call getmem(need(6),pz2,ngot(6),'m811',0)
c
         write(ioutpt,1040) need(6)
         if(ntriang.le.0) call lnkerr('not enough memory')
c
         continue
c
         call trn2e(z2(values),nnp,ntriang,z1(c),nbf,norbs,
     #              z2(t1),z2(t2),numij,z2(asort),lnsort,ngroup,
     #              nmax,nsym,a1(ijgrp),a1(ijadd),a1(kadd),a1(ladd),
     #              a1(orbsym),z1(hmo),levfrm,
     #              z2(val),a2(lab),a2(bin),lenbin,toguga,z2(values),
     $              tocan,ops)
         call getmem(-ngot(6),pz2,idum,'m811',idum)
         call getmem(-ngot(1),pa1,idum,'m811',idum)
c
c     ----- and exit with grace -----
c
         if(tocan) then
            call iosys('close tints',0,0,0,' ')
         endif
         if(toguga) then
            call iosys('close gints',0,0,0,' ')
         endif
      endif
      call chainx(0)
c
c
      stop
      end
