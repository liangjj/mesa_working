*deck m6002
c***begin prologue     m6002
c***date written       890407   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           m6002, link 6003, direct potential
c***author             schneider, barry (lanl)
c***source             m6002
c***purpose            driver for numerical direct potential
c***description        calculates static potential on a numerical grid
c***                   using technology of nuclear attraction integrals
c*** 
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       m6002
      program vstat
      implicit integer (a-z)
      parameter (dimpr=300 , dimcen=10 , dimst=25 , dimst2=325)
      parameter (mmax=2, lmax=2, maxmom=20)
      logical logkey, logky, grdtyp, vmom, posinp
      character *4096 ops
      character *1600 card
      character *13 grdnam
      character *8 chrkey, cpass, whrden
      character *128 filkne, filgrd, filpot
      character *3 itoc, ans, aosym, group, psym, sym
      character *3 roexst
      character *12 diag
      real *8 z, alf, cont, rloc, pi, piterm, pitern
      real *8 charge, fpkey, acrcy, scale, sumpot, reschg
c----------------------------------------------------------------------c
c                unicos memory management                              c
      common a(1)
      dimension z(1)
      common /memory / ioff 
      equivalence (z,a)
c----------------------------------------------------------------------c
      common /io/ inp,iout
      common/chrpss/ aosym(dimpr)
      common /aosi/ npr, ncon, nxyzc(dimpr,4), nprc(dimpr), 
     1              sncon, smllst(dimpr,2)
      common /aosr/ alf(dimpr), cont(dimpr)
      common /rloc/ charge(dimcen), rloc(3,dimcen)
      common/logprt/logky(10)
      common /nmbrs /  pi, piterm, pitern, acrcy, scale, icanon
      dimension  psym(dimst2), sym(dimst)
      dimension diag(dimst2), roexst(dimst2), sumpot(1200)
      data pi /3.14159265358979323846d+00/
      piterm=2.d+00/pi**0.5d+00
      pitern=pi**1.5d+00
      call drum
      call iosys ('read character options from rwf',-1,0,0,ops)
      if( posinp ('$vstat',cpass) ) then
          call cardin (card)
      endif
      call iosys ('read character "kohn data filename" from rwf',-1,0,0,
     1             filkne)
      call iosys ('open kohndt as old',0,0,0,filkne)
c----------------------------------------------------------------------c
c                    set options                                       c
c----------------------------------------------------------------------c 
      logky(1)=logkey(ops,'print=m6002=vectors',.false.,' ')
      logky(2)=logkey(ops,'print=m6002=initial-density-matrix',
     1               .false.,' ')
      logky(3)=logkey(ops,'print=m6002=final-density-matrix',.false.,
     1                ' ')
      logky(4)=logkey(ops,'no-ssd',.false.,' ')
      logky(5)=logkey(ops,'print=m6002=potential',.false.,' ')
      logky(6)=logkey(ops,'print=m6002=basis',.false.,' ')
      write (iout,1020)
c---------------------------------------------------------------------c
c            initialize x y and z values for cartesian aos            c
c---------------------------------------------------------------------c
      call nlmxyz
c---------------------------------------------------------------------c
c                read in basic data                                   c
c---------------------------------------------------------------------c
      vmom=logkey(card,'long-range-moments',.false.,' ')
      grdtyp=logkey(card,'untransformed-grid',.false.,' ')
      grdnam='"trns grid"'
      if (grdtyp) then
          grdnam='"untrns grid"'
      endif
      call iosys ('read integer "no. channels" from kohndt',1,nsts,
     1             0,' ')
      nstri=nsts*(nsts+1)/2
      call iosys ('read character "grid filename" from rwf',-1,0,0,
     1             filgrd)
c----------------------------------------------------------------------c
c               open grid file and get no. of points                   c
c               need four arrays of dimension pntbuf                   c
c----------------------------------------------------------------------c
      call iosys ('open grid as old',0,0,0,filgrd)
      call iosys ('read integer "no. grid pts" from grid',1,npnts,0,
     1            ' ')
      call iosys ('read integer "point buffer" from grid',1,
     1            pntbuf,0,' ')
      if (vmom) then
          call iosys ('read integer ilong from grid',1,ilong,0,' ')
      else
          ilong=npnts+1
      endif
      call iosys ('read character "potential filename" from rwf',-1,0,0,
     1             filpot)
      call iosys ('read real charge from kohndt',1,reschg,0,' ')
      if (reschg.ne.0.d+00) then
          write (iout,1040) reschg
      endif
      whrden=chrkey(card,'density-matrix-input','kohndt',' ')
c----------------------------------------------------------------------c
c                read in geometry                                      c
c----------------------------------------------------------------------c
      call iosys ('read "number of atoms" from kohndt',1,ncen,0,' ')
      call iosys ('read real "nuclear charges" from kohndt',ncen,charge,
     1             0,' ')
      call iosys ('read integer "no. primitives" from kohndt',1,
     1             npr,0,' ')
      call iosys ('read integer "no. contracted" from kohndt',1,
     1             ncon,0,' ')
      call iosys ('read integer "total no. mos" from kohndt',1,nmo,0,
     1             ' ')
      sumpot(1)=0.0d0      
      do 4 i=1,ncen
         sumpot(1)=sumpot(1)+charge(i)
         call iosys ('read real "x-y-z atom-'//itoc(i)//'" from kohndt',
     1               3,rloc(1,i),0,' ')
    4 continue
      sumpot(1)=sumpot(1)-reschg
      nel=sumpot(1)
      size=npnts*nstri
      call iosys ('open vstat as new on ssd',size,0,0,filpot)
      if ( posinp('$geom',cpass) ) then
           call cardin(card)
           ncen=intkey(card,'no-centers',1,' ')
           do 10 i=1,ncen
              call fparr(card,'x-y-z-atom-'//itoc(i),rloc(1,i),3,' ')
              charge(i)=fpkey(card,'charge-atom-'//itoc(i),1.d+00,' ')
   10      continue
      endif
c----------------------------------------------------------------------c
c                 read in basis set and density matrix                 c
c                            information                               c
c----------------------------------------------------------------------c
      if ( posinp('$basis',cpass) ) then
           call cardin(card)
           group=chrkey(card,'group','c2v',' ')
           call locase(group,group)
           npr=intkey(card,'no-prim',1,' ')
           ncon=intkey(card,'no-contracted',1,' ')
           nmo=intkey(card,'no-mos',ncon,' ')
           nel=intkey(card,'no-target-electrons',1,' ')
           call intarr (card,'no-prim/cont',nprc,ncon,' ')
c----------------------------------------------------------------------c
c             old code will use someday when symmetry in               c
c                                                                      c
c          do 20 i=1,nsts
c             sym(i)=chrkey(card,'symmetry-state-'//itoc(i),'a1',' ')
c             call locase(sym(i),sym(i))
c   20     continue
c          ij=0
c          do 30 i=1,nsts
c             do 40 j=1,i
c                ij=ij+1
c                if (group.eq.'c2v') then
c                    call c2v(sym(i),sym(j),psym(ij))
c                elseif (group.eq.'d2h') then
c                    call d2h(sym(i),sym(j),psym(ij))
c                else
c                    call lnkerr('group type screwed up')
c                endif
c                diag(ij)='off diagonal'
c   40        continue
c             diag(ij)='diagonal'
c   30     continue
c     write (iout,1030) (psym(ij),ij=1,nstri)
      call iosys ('write character group to vstat',0,0,0,group)
      call iosys ('write character symmetries to vstat',0,0,0,psym)
      call iosys ('write character "diagonal ?" to vstat',0,0,0,diag)
c----------------------------------------------------------------------c
      endif
      call iosys ('write character "grid type" to vstat',0,0,0,grdnam)
      call iosys ('write integer "no. states" to vstat',1,nsts,0,' ')
c----------------------------------------------------------------------c
c                 read in the basis set parameters                     c
c----------------------------------------------------------------------c
      call rdbsis(ntot,nprmx,aosym)
      write (iout,2000)
      do 90 i=1,ncen
         write (iout,2010) (rloc(j,i),j=1,3), charge(i)
   90 continue
      nmotri=nmo*(nmo+1)/2
      ncntri=ncon*(ncon+1)/2
      write (iout,2020) npr,ncon
      write (iout,2030) (nprc(i),i=1,ncon)
      if (logkey(ops,'print=m6002=basis',.false.,' ')) then
          call wrtbas (ntot)
      endif
      if (ntot.ne.npr) call lnkerr('no. primitives error')
c----------------------------------------------------------------------c
c                 calculate memory requirements                        c
c----------------------------------------------------------------------c
      call iosys ('read integer maxsiz from rwf',1,canget,0,' ')
c     call getscm(0,z,canget,'how much core',0)
      if (vmom) then
          words0=ncntri*nstri+iadtwp(3*nstri*maxmom+nstri)+nstri*maxmom
          nadd=pntbuf
      else
          words0=0
          nadd=ncntri*nstri
      endif
      words1=words0+nadd
      words0=words0+nadd+ncon*nmo+nmotri+nmo*nmo+nmo*ncon+ncon*ncon
      words1=words1+(37+nstri)*pntbuf+nstri+iadtwp(2*pntbuf)
      words=max(words0,words1)
      if (words.gt.canget) call lnkerr('cannot get required memory:'//
     1                                 'will quit')
      call iosys ('write integer maxsiz to rwf',1,words,0,' ')
      call getscm(z,words,ngot,'m6002',0)
      rhocn=ioff
      if (vmom) then
          momind=rhocn+ncntri*nstri
          nomom=wpadti(momind+3*nstri*maxmom)
          mom=nomom+nstri
          ylm=iadtwp(mom+nstri*maxmom)
          nadd=pntbuf
      else
          ylm=rhocn
          nadd=ncntri*nstri
      endif
      vec=ylm+nadd
      rhomo=vec+ncon*nmo
      temp=rhomo+nmotri
      temp1=temp+nmo*nmo
      temp2=temp1+nmo*ncon
      grid=vec
      valint=grid+4*pntbuf
      pc=valint+nstri*pntbuf
      g=pc+3*pntbuf
      fvec=g+21*pntbuf
      arg=fvec+7*pntbuf
      facr=arg+pntbuf
      pointr=wpadti(facr+pntbuf)
c----------------------------------------------------------------------c
c           calculate number of passes of grid file needed             c
c----------------------------------------------------------------------c
      npass=npnts/pntbuf
      nolst=npnts-npass*pntbuf
      if (nolst.ne.0) then
          npass=npass+1
      else
          nolst=pntbuf
      endif
c---------------------------------------------------------------------c
c          calculate which buffer load begins long range region       c
c---------------------------------------------------------------------c
      iright=0
      nwd=pntbuf
      do 95 region=1,npass
         if (region.eq.npass) nwd=nolst
         iright=iright+nwd
         if (iright.ge.ilong) go to 96
   95 continue
      passlr=npass+1
      go to 97
   96 passlr=region
   97 write(iout,3000) pntbuf, npass, words
      if (passlr.le.npass) then
          write (iout,3001) passlr
      endif     
c----------------------------------------------------------------------c
c        read in transformation vectors and mo density matrix          c
c                 transform to contracted ao basis                     c
c        when done only need rhocn for rest of calculation             c
c----------------------------------------------------------------------c
      call vecin(z(vec),nmo,ncon,logky(1))
      call dmatin (z(rhocn),z(rhomo),z(vec),z(temp),z(temp1),z(temp2),
     1             roexst,smllst,nel,ncon,nmo,nsts,nstri,nmotri,dimpr,
     2             ncntri,whrden,logky(2))
c----------------------------------------------------------------------c
c           read in long range moments if required                     c
c----------------------------------------------------------------------c
      if (vmom) then
          call rdmom(z(mom),a(nomom),a(momind),nsts,nstri,maxmom)
      endif    
      call iosys ('create real "static potential" on vstat',size,0,
     1            0,' ')
c---------------------------------------------------------------------c
c                   generate error function table                     c
c---------------------------------------------------------------------c
      call generf
      write (iout,4000)
      nwds=0
      nwrite=0
      call rzero(sumpot,1200)
      cntpnt=0
      callsr=0
      calllr=0
      do 200 ipass=1,npass
         noptrg=pntbuf
         if (ipass.eq.npass) then
             noptrg=nolst
         endif
c----------------------------------------------------------------------c
c           read in this block of grid points                          c
c----------------------------------------------------------------------c
         call iosys ('read real '//grdnam//' from grid without '//
     1               'rewinding',4*noptrg,z(grid),0,' ')   
         call rzero(z(valint),noptrg*nstri)
         if (ipass.lt.passlr) then
             callsr=callsr+1
             call vints (z(valint),z(rhocn),roexst,z(grid),z(g),z(pc),
     1                   z(fvec),z(facr),z(arg),a(pointr),reschg,diag,
     2                   sumpot,ncen,noptrg,nsts,nstri,ncntri,noptrg)
         elseif (ipass.eq.passlr) then
             nop0=ilong-cntpnt-1
             callsr=callsr+1
             call vints (z(valint),z(rhocn),roexst,z(grid),z(g),z(pc),
     1                   z(fvec),z(facr),z(arg),a(pointr),reschg,diag,
     2                   sumpot,ncen,nop0,nsts,nstri,ncntri,noptrg)
             nop1=cntpnt+noptrg-ilong+1
             calllr=calllr+1
             call vlr (z(valint+nop0),z(ylm),z(mom),a(nomom),a(momind),
     1                 z(grid+4*nop0),nop1,nstri,noptrg,maxmom)
         elseif (ipass.gt.passlr) then
             callsr=callsr+1
             call vints (z(valint),z(rhocn),roexst,z(grid),z(g),z(pc),
     1                   z(fvec),z(facr),z(arg),a(pointr),reschg,diag,
     2                   sumpot,ncen,noptrg,nsts,nstri,ncntri,noptrg)
         endif
         call wrtpot(z(valint),z(grid),noptrg,nstri,noptrg,nwrite,nwds,
     1               ipass,logky(5))      
         cntpnt=cntpnt+noptrg
  200 continue
      if (nwds.ne.size) call lnkerr('error in words to vstat')
      write (iout,5000) nwrite
      write (iout,6000) nwds
      write (iout,7000) (sumpot(i),i=1,nstri)
      write (iout,8050) callsr, calllr
      call iosys ('rewind all on grid read-and-write',0,0,0,' ')
      call iosys ('rewind all on vstat read-and-write',0,0,0,' ')
      call iosys ('close grid',0,0,0,' ')
      call iosys ('close vstat',0,0,0,' ')
      call iosys ('write integer maxsiz to rwf',1,canget,0,' ')
      call chainx(0)
      stop
 1000 format(a80)
 1020 format(//,20x,'*****  m6002:numerical static potential program ***
     1**')
 1030 format (/,5x,'symmetry of density matrices',(/,10(a3,1x)))
 1040 format (/,5x,'residual charge is',1x,f10.5)
 2000 format (//,10x,'position and charge of atomic centers',
     1        //,7x,'x',17x,'y',17x,'z',17x,'charge')
 2010 format (/,1x,e15.8,3x,e15.8,3x,e15.8,3x,e15.8)
 2020 format (//,5x,'data on atomic orbitals',//,5x,'no. prim.',
     1 2x,i4,5x,'no. of cont.',2x,i4)
 2030 format (//,5x,'no. primitives / contracted function',(/,10x,20(1x
     1 ,i2)))
 3000 format (//,5x,'point buffer',2x,i6,2x,'no. passes of grid file',
     1        1x,i3,1x,'no. words memory required',1x,i8)
 3001 format (/,5x,'long range region begins at pass',1x,i4)
 3010 format (//,5x,'no. writes to pot',2x,i4)
 3020 format (//,5x,'no. words written to pot',2x,i8)
 4000 format (/,5x,'pass',5x,'density matrix',5x,'words written')
 5000 format(/,5x,'no. disk writes',1x,i8)
 6000 format(/,5x,'no. words written',1x,i8)
 7000 format (/,5x,'potential sum',(/,5x,5e15.8))
 8050 format (/,5x,'calls to vints',1x,i3,2x,'calls to vlr',2x,i3)
      end
