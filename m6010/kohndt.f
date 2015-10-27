*deck @(#)kohndt.f	1.1 9/7/91
c***begin prologue     m6010
c***date written       890629   (yymmdd)
c***revision date      921113   (yymmdd)
c                      revised to read data from mesa which is then
c                      made compatible with bis kohn codes.
c                      revised to delete channels.                   
c***keywords           m6010, link 6010, kohn data
c***author             schneider, barry (lanl)
c***source             m6010
c***purpose            enter data for kohn calculation
c***
c
c***references
c
c***routines called    iosys, util and mdutil
c***end prologue       m6010
      program kohndt
      implicit integer (a-z)
      character *1600 card
      character *8  cpass, presnt, chrkey, drctv
      character *16 fptoc
      character *16 bflabl(2000)
      character *24 ftitl
      character *4096 ops
      character *128 namkne, filkne
      character *32 xform
      character *4 itoc
      character *80 title
      logical logkey, prbaso, prbasn, propt, prden, provlp, prkn
      logical prpot, prdir1, prdir2, prdir, prvec, diagdr
      logical noopt, rearr, rearrc, cmpopt, vdrctv
      real *8 z, fpkey, rmin, rmax, ksqmax, echan, energy, reschg
      real *8 ene
      dimension echan(100), energy(100), list(300), listc(100)
      common /io/ inp, iout
      pointer(p,z(1)),(p,a(1))
      call drum
      call iosys ('read character options from rwf',-1,0,0,ops)
      prbaso=logkey(ops,'print=m6010=old-basis',.false.,' ')
      prbasn=logkey(ops,'print=m6010=new-basis',.false.,' ')
      prvec=logkey(ops,'print=m6010=transformation-matrix',.false.,' ')
      propt=logkey(ops,'print=m6010=optical-potential',.false.,' ')
      prden=logkey(ops,'print=m6010=density-matrix',.false.,' ')
      provlp=logkey(ops,'print=m6010=overlap',.false.,' ')
      prdir1=logkey(ops,'print=m6010=1-electron-direct-hamiltonian',
     1              .false.,' ')
      prdir2=logkey(ops,'print=m6010=2-electron-direct-hamiltonian',
     1              .false.,' ')
      prdir=logkey(ops,'print=m6010=direct-hamiltonian',.false.,' ')
      diagdr=logkey(ops,'m6010=diagonalize',.false.,' ')
      noopt=logkey(ops,'kohn=noopt',.false.,' ')
      vdrctv=logkey(ops,'only-one-electron',.false.,' ')
      prkn=logkey(ops,'print=m6010=kinetic-energy',.false.,' ')
      prpot=logkey(ops,'print=m6010=potential-energy',.false.,' ')      
      cmpopt=logkey(ops,'m6010=complex-optical-potential',.false.,' ')
      call iosys ('read character "kohn filename" from rwf',-1,
     1             0,0,namkne)
      call iosys ('open kohn as old',0,0,0,namkne)
      write (iout,500)
      if (cmpopt) then
          write(iout,501)
      endif
      drctv='all'
      if (vdrctv) then
          drctv='one'
          write(iout,502)   
      endif
      call posinp('$kohndt',cpass)
      call cardin(card)
c     filkne=chrkey(card,'kohn-data-file-name','kohndt',' ')
c     call locase(filkne,filkne)
c     call iosys ('write character "kohn data filename" to rwf',0,0,0,
c    1             filkne)
      call iosys ('read character "kohn data filename" from rwf',-1,
     1             0,0,filkne)
      pntbuf=intkey(card,'point-buffer',20000,' ')
c----------------------------------------------------------------------c
c      get channel and scattering energies. put on kohndt              c
c----------------------------------------------------------------------c
      call iosys ('open kohndt as new on ssd',262144,0,0,filkne)
      nchano=intkey(card,'no-open-channels',1,' ')
      nchanc=intkey(card,'no-closed-channels',0,' ')
      nchan=nchano+nchanc
      ncokpt=intkey(card,'no-open-channels-kept',nchano,' ')
      ncckpt=intkey(card,'no-closed-channels-kept',nchanc,' ')
      nchann=ncokpt+ncckpt
      call iosys ('write character "integral directive" to kohndt',0,0,
     1             0,drctv)
      call iosys ('write integer "no. channels" to kohndt',1,
     1             nchann,0,' ')
      call iosys ('read integer "no. energies" from kohn',1,nener,0,' ')
      do 10 i=1,nchan
         call iosys ('read real "ci energy '//itoc(i)//'" from kohn',
     1               1,echan(i),0,' ')
   10 continue
      call iosys ('write integer "no. energies" to kohndt',1,
     1             nener,0,' ')
      call izero(listc,nchan)
      if (nchann.ne.nchan) then
          if (ncokpt.ne.nchano) then
              write(iout,1501)
              call intarr(card,'open-channels-kept',list,ncokpt,' ')
              write(iout,1601) (list(i),i=1,ncokpt)
              do 1800 i=1,ncokpt
                 listc(list(i))=1
 1800         continue
              rearrc=.true.
          endif
          if (ncckpt.ne.nchanc) then
              write(iout,1502)
              call intarr(card,'closed-channels-kept',list,ncckpt,' ')
              write(iout,1602) (list(i),i=1,ncckpt)
              do 1801 i=1,ncckpt
                 listc(list(i))=1
 1801         continue
              rearrc=.true.
          endif
      else
          rearrc=.false.
      endif
      if (rearrc) then
          count=0
          do 1605 i=1,nchan
             if(listc(i).eq.1) then
                count=count+1
                listc(i)=count 
                echan(count)=echan(i)
             endif
 1605     continue
          if (count.ne.nchann) then
              call lnkerr('error in channel count')
          endif
      endif
      call iosys ('write real "chan energies" to kohndt',nchann,
     1             echan,0,' ')
      call iosys ('read real "scattering energies" from kohn',nener,
     1            energy,0,' ')
c----------------------------------------------------------------------c
c               convert from hartrees to rydbergs                      c
c----------------------------------------------------------------------c
      do 1 iene=1,nener
         energy(iene)=2.d+00*energy(iene)
    1 continue 
      call iosys ('write real "scatt energies" to kohndt',nener,
     1             energy,0,' ')
c----------------------------------------------------------------------c
c             get basis function information                           c
c----------------------------------------------------------------------c
      call iosys('read integer "number of basis functions"'//
     $ ' from rwf',1,nbf,0,' ')
      call iosys('read integer npvec from kohn',1,npvec,0,' ')
      nsmall=intkey(ops,'kohn=nsmall',nsmall,' ')
      nl2=intkey(ops,'kohn=nl2',npvec,' ')
      norb=npvec/nchan
      if (norb.ne.nl2) then
          call lnkerr('error in orbitals when deleting channels')
      endif
      npvecn=nl2*nchann
      call iosys ('write integer "size of p space" to kohndt',1,
     1             npvecn,0,' ')
      write(iout,700) nsmall,npvec,npvecn
      call iosys ('write integer nsmall to kohndt',1,nsmall,0,' ')
      nmotot=nl2+nsmall
      write(iout,800) nl2,nmotot
      if (nmotot.ne.nbf) then
          write (iout,1500)
          call intarr(card,'orbitals-kept',list,nmotot,' ')
          write (iout,1600) (list(i),i=1,nmotot) 
          rearr=.true.
      else
          rearr=.false.
      endif
      call iosys ('write integer "no. l2 orbitals" to kohndt',1,nl2,
     1             0,' ')
      call iosys ('write integer "total no. mos" to kohndt',1,nmotot,
     1             0,' ')
      if(logkey(ops,'drop',.false.,' ')) then
         idrop=1
         call iosys('read integer "old nbf" from kohn',
     $              1,oldnbf,0,' ')
         write (iout,900)
      else
         oldnbf=nbf
         idrop=0
         write (iout,1000)
      end if
c----------------------------------------------------------------------c
c                 get the basis function labels                        c
c----------------------------------------------------------------------c
      call iosys('read character "basis function labels" from rwf',
     $          -1,0,0,bflabl)
      call iosys('read integer "number of atoms" from rwf',1,nat,0,' ')
      call iosys('write integer "number of atoms" to kohndt',1,
     1            nat,0,' ')
      call iosys('read integer "number of basis functions" from rwf',
     $     1,nbasis,0,' ')
      call iosys('length of exponents on rwf',nprim,0,0,' ')
      call iosys('length of "contraction coefficients" on rwf',
     $     ncont,0,0,' ')
      call iosys('length of "number of pure functions" on rwf',
     $     ntypes,0,0,' ')
      call iosys('read integer "number basis types" from rwf',
     $     1,nbtype,0,' ')
      call iosys('length of "power of x" on rwf',lenxyz,0,0,' ')
c----------------------------------------------------------------------c
c              divide core for basis set information                   c
c----------------------------------------------------------------------c
      nnpbs=(nbasis+1)*nbasis/2
      nwreal=oldnbf+4*nat+ncont+nprim
      nwint=5*ntypes*nat+5*ntypes+3*lenxyz
      need0=nwreal+iadtwp(nwint)
c----------------------------------------------------------------------c
c              divide core for matrix output                           c
c----------------------------------------------------------------------c
      nnp=nbf*(nbf+1)/2
      call iosys('read integer nwks from rwf',1,mdim,0,' ')
      fac=1
      if (cmpopt) then
          fac=2
      endif 
      need1=4*nbf*nbf+fac*npvec*npvec+npvec*npvec+2*nnp
      if (rearrc) then
          need1=need1+npvec*npvec
      endif
c----------------------------------------------------------------------c
c               get larger of two words                                c
c----------------------------------------------------------------------c
      words=max(need0,need1)
      words=wptoin(words)
      call iosys ('read integer mxcore from rwf',1,maxcor,0,' ')
      if (words.gt.maxcor) then
          call lnkerr('requested too much core')
      else
          call memory(words,p,ngot,'m6010',0)
      endif
c----------------------------------------------------------------------c
c            put out transformation matrix                             c
c                        and other                                     c
c                    required matrices                                 c
c----------------------------------------------------------------------c
      c=1
      scr=c+nbf*nbf
      s=scr+nbf*nbf
      vec=s+nbf*nbf
      hopt=vec+nbf*nbf
      hpp=hopt+fac*npvec*npvec
      dir1=hpp+npvec*npvec
      dir2=dir1+nnp
      if (rearrc) then
          remat=dir2+nnp
      endif
c----------------------------------------------------------------------c
c             all matrices are written in ao and mo form               c
c                          to kohndt                                   c
c----------------------------------------------------------------------c
      call iosys('read character "transformation vector" from kohn',
     $           -1,0,0,xform)
      call iosys('write character "transformation vector" to kohndt',
     $            0,0,0,xform)
      call iosys('read real '//xform//' from kohn',nbf*nbf,z(c),
     1            0,' ')
      call copy(z(c),z(s),nbf*nbf)
      if (rearr) then
          call newmat(z(c),z(scr),z(scr),list,nbf,nmotot,'ao-mo')
      endif
      if (prvec) then
          title='transformation vectors'
          call prntrm(title,z(c),nbf,nmotot,nbf,nmotot,iout)
      endif
      call iosys ('write real '//xform//' to kohndt',nbf*nmotot,
     1             z(c),0,' ')
      call copy(z(s),z(c),nbf*nbf)
      call iosys ('read real "overlap integrals" from kohn',nnp,z(scr),
     1             0,' ')
      call trtosq(z(s),z(scr),nbf,nnp)
      if (provlp) then
          title='overlap integrals'
          call prntrm(title,z(s),nbf,nbf,nbf,nbf,iout)
      endif
      call iosys ('write real "overlap integrals" to kohndt',nbf*nbf,
     1             z(s),0,' ')
      call iosys ('read real "kinetic integrals" from kohn',
     1             nnp,z(scr),0,' ')
      call trtosq(z(s),z(scr),nbf,nnp)
      if (prkn) then
          title='kinetic energy integrals'
          call prntrm(title,z(s),nbf,nbf,nbf,nbf,iout)
      endif
      call iosys ('write real "kinetic integrals" to kohndt',
     1             nbf*nbf,z(s),0,' ')
c----------------------------------------------------------------------c
c         put out kinetic energy integrals in mo form they may         c
c         be needed later                                              c
c----------------------------------------------------------------------c
      call ebc(z(vec),z(s),z(c),nbf,nbf,nbf)
      call ebtc(z(scr),z(c),z(vec),nbf,nbf,nbf)
      if (rearr) then
          call newmat(z(scr),z(vec),z(vec),list,nbf,nmotot,'mo-mo')
      endif
      call iosys ('write real "mo kinetic energy integrals" to kohndt',
     1             nmotot*nmotot,z(scr),0,' ')
      call iosys ('read real "potential integrals" from kohn',
     1             nnp,z(scr),0,' ')
      call trtosq(z(s),z(scr),nbf,nnp)
      if (prpot) then
          title='potential energy integrals'
          call prntrm(title,z(s),nbf,nbf,nbf,nbf,iout)
      endif
      call iosys ('write real "potential integrals" to kohndt',
     1             nbf*nbf,z(s),0,' ')
c----------------------------------------------------------------------c
c            get ao one-electron matrix and transform                  c
c                       to mo basis                                    c
c----------------------------------------------------------------------c
      call iosys('read real "direct ao 1e kohn hamiltonian"'//
     $' from kohn',nnp,z(dir1),0,' ')
      call trtosq(z(scr),z(dir1),nbf,nnp)
      call iosys ('write real "ao 1e kohn hamiltonian" to kohndt',
     1             nbf*nbf,z(scr),0,' ')
      if (prdir1) then
          title='ao 1e kohn hamiltonian'
          call prntrm(title,z(scr),nbf,nbf,nbf,nbf,iout)
      endif
      call ebc(z(vec),z(scr),z(c),nbf,nbf,nbf)
      call ebtc(z(scr),z(c),z(vec),nbf,nbf,nbf)
      if (rearr) then
          call newmat(z(scr),z(vec),z(vec),list,nbf,nmotot,'mo-mo')
      endif
      call iosys ('write real "mo 1e kohn hamiltonian" to kohndt',
     1             nmotot*nmotot,z(scr),0,' ')
      if (prdir1) then
          title='mo 1e kohn hamiltonian'
          call prntrm(title,z(scr),nmotot,nmotot,nmotot,nmotot,iout)
      endif
c----------------------------------------------------------------------c
c         now do same with part depending on density matrix            c
c----------------------------------------------------------------------c
      call iosys('rewind "direct ao 2e kohn hamiltonian" on kohn',
     $ 0,0,0,' ')
      do 20 i=1,nchan
         do 30 j=1,i
            if (rearrc) then
                ii=listc(i)
                jj=listc(j)
                prd=ii*jj
            else
                ii=i
                jj=j
                prd=1
            endif
                call iosys('read real "ao t1pdm:'//itoc(i)//
     1                      itoc(j)//'" from kohn',nbf*nbf,z(scr),0,' ')
                if (prd.ne.0) then
                    call iosys('write real "ao t1pdm:'//itoc(ii)//
     1                          itoc(jj)//'" to kohndt',nbf*nbf,
     2                          z(scr),0,' ')
                    if (prden) then
                        title='ao density matrix for '//
     1                        'state 1 ='//itoc(ii)//' state 2'//
     2                        ' ='//itoc(jj)
                        call prntrm(title,z(scr),nbf,nbf,nbf,nbf,iout)
                    endif
                endif
                call iosys('read real "mo t1pdm:'//itoc(i)//
     1                      itoc(j)//'" from kohn',nbf*nbf,z(scr),0,' ')
                if (prd.ne.0) then
                    if (rearr) then
                        call newmat(z(scr),z(vec),z(vec),list,nbf,
     1                              nmotot,'mo-mo')
                    endif
                    call iosys('write real "mo t1pdm:'//itoc(ii)//
     1                          itoc(jj)//'" to kohndt',nmotot*nmotot,
     2                          z(scr),0,' ')
                    if (prden) then
                        title='mo density matrix for state 1 ='//
     1                         itoc(ii)//' state 2 ='//itoc(jj)
                        call prntrm(title,z(scr),nmotot,nmotot,nmotot,
     1                              nmotot,iout)
                    endif
                endif
                call iosys('read real "direct ao 2e kohn hamiltonian"'//
     1                     ' from kohn without rewinding',nnp,z(dir2),
     2                       0,' ')
                if (prd.ne.0) then
                    call trtosq(z(scr),z(dir2),nbf,nnp)
                    call iosys('write real "ao 2e kohn hamiltonian'//
     1                          itoc(ii)//itoc(jj)//'" to kohndt',
     2                          nbf*nbf,z(scr),0,' ')
                    if (prdir2) then
                        title='ao 2 electron kohn hamiltonian for '//
     1                        'state 1 ='//itoc(ii)//' state 2 ='//
     2                         itoc(jj)
                        call prntrm(title,z(scr),nbf,nbf,nbf,nbf,iout)
                    endif
                    call ebc(z(vec),z(scr),z(c),nbf,nbf,nbf)
                    call ebtc(z(scr),z(c),z(vec),nbf,nbf,nbf)
                    if (rearr) then
                        call newmat(z(scr),z(vec),z(vec),list,nbf,
     1                              nmotot,'mo-mo')
                    endif
                    call iosys('write real "mo 2e kohn hamiltonian'//
     1                          itoc(ii)//itoc(jj)//'" to kohndt',
     2                          nmotot*nmotot,z(scr),0,' ')
                    if (prdir2) then
                        title='mo 2 electron kohn hamiltonian for '//
     1                        'state 1 ='//itoc(ii)//' state 2 ='//
     2                         itoc(jj)
                        call prntrm(title,z(scr),nmotot,nmotot,nmotot,
     1                              nmotot,iout)
                    endif
c----------------------------------------------------------------------c
c      if channel i = channel j add one-electron contribution to       c
c      direct hamiltonian                                              c
c----------------------------------------------------------------------c
                    if(ii.eq.jj) then
                       do 40 iii=1,nnp
   40                     z(dir2+iii-1) = z(dir2+iii-1) + z(dir1+iii-1)
                    endif
                    call trtosq(z(scr),z(dir2),nbf,nnp)
                    call iosys('write real "ao kohn hamiltonian'//
     1                          itoc(ii)//itoc(jj)//'" to kohndt',
     2                          nbf*nbf,z(scr),0,' ')
                    if (prdir) then
                        title='ao kohn hamiltonian for state 1 ='//
     1                         itoc(ii)//' state 2 ='//itoc(jj)
                        call prntrm(title,z(scr),nbf,nbf,nbf,nbf,iout)
                    endif
                    call ebc(z(vec),z(scr),z(c),nbf,nbf,nbf)
                    call ebtc(z(scr),z(c),z(vec),nbf,nbf,nbf)
                    if (rearr) then
                        call newmat(z(scr),z(vec),z(vec),list,nbf,
     1                              nmotot,'mo-mo')
                    endif
                    call iosys('write real "mo kohn hamiltonian'//
     1                          itoc(ii)//itoc(jj)//'" to kohndt',
     2                          nmotot*nmotot,z(scr),0,' ')
                    if (prdir) then
                        title='mo kohn hamiltonian for state 1 ='//
     1                         itoc(ii)//' state 2 ='//itoc(jj)
                        call prntrm(title,z(scr),nmotot,nmotot,nmotot,
     1                              nmotot,iout)
                    endif
                    if (diagdr) then
                        call tred2(nmotot,nmotot,z(scr),z(vec),
     1                             z(s),z(scr))
                        call tql2(nmotot,nmotot,z(vec),z(s),z(scr),ierr)
                        title='eigenvalues of direct hamiltonian i='
     1                        //itoc(ii)//' j='//itoc(jj)
                        write (iout,2) title
                        write (iout,3) (z(ii),ii=vec,vec+nmotot-1)
                    endif
                endif
   30    continue
   20 continue
c----------------------------------------------------------------------c
c         now output the optical potential and hpp                     c
c----------------------------------------------------------------------c
      call iosys('rewind hopt on kohn',0,0,0,' ')
      call iosys('rewind hpp on kohn',0,0,0,' ')
      ntot = npvec*npvec
      mtot=npvecn*npvecn
      do 50 iene=1,nener
         call iosys('read real hpp from kohn without rewinding',
     x               ntot,z(hpp),0,' ')
         if (rearrc) then
             call trimc(z(hpp),z(hpp),z(remat),listc,npvec,npvecn,
     x                  nchan,nl2)
         endif  
         if (.not.noopt) then
              call iosys('read real hopt from kohn without rewinding',
     x                    ntot,z(hopt),0,' ')
              if (rearrc) then
                  call trimc(z(hopt),z(hopt),z(remat),listc,npvec,
     x                       npvecn,nchan,nl2)
              endif  
         else
              call rzero(z(hopt),mtot)
         endif
c----------------------------------------------------------------------c
c          mesa computes hopt = hpq (hqq - e)**-1 hqp                  c
c          therefore change its sign here to agree with standard       c
c                         convention.                                  c
c----------------------------------------------------------------------c
         do 60 iip = 1, mtot
            z(hopt+iip-1) = -z(hopt+iip-1)
   60    continue
         ene=echan(1)+.5d+00*energy(iene)
         if (propt) then
             title='hpp energy ='//fptoc(ene)
             call prntrm(title,z(hpp),npvecn,npvecn,npvecn,npvecn,iout)
             title='vpp energy ='//fptoc(ene)
             call prntrm(title,z(hopt),npvecn,npvecn,npvecn,npvecn,iout)
         endif
         ftitl='h(pp)-'//fptoc(ene)
         call iosys ('write real '//ftitl//' to kohndt',mtot,
     1               z(hpp),0,' ')
c----------------------------------------------------------------------c
c              add hpp and hopt together and store                     c
c              if a complex optical potential is needed later          c
c              leave the space on the file to store it                 c
c----------------------------------------------------------------------c
         do 65 iip=1,mtot
            z(hopt+iip-1)=z(hopt+iip-1)+z(hpp+iip-1)
   65    continue
         ftitl='v(pp)-'//fptoc(ene)
         call iosys ('write real '//ftitl//' to kohndt',fac*mtot,
     1               z(hopt),0,' ')
         if (propt) then
             title='hpp + vpp  energy ='//fptoc(ene)
             call prntrm(title,z(hopt),npvecn,npvecn,npvecn,npvecn,iout)
         endif
   50 continue
c----------------------------------------------------------------------c
c           get basis function information                             c
c           from rwf file and put in new format                        c
c----------------------------------------------------------------------c
      nnp=nnpbs
      index=1
      zan=index+oldnbf
      c=zan+nat
      cont=c+3*nat
      ex=cont+ncont
      ptprim=wpadti(ex+nprim)
      noprim=ptprim+ntypes*nat
      nocont=noprim+ntypes*nat
      ptcont=nocont+ntypes*nat
      start=ptcont+ntypes*nat
      nocart=start+ntypes*nat
      nobf=nocart+ntypes
      minmom=nobf+ntypes
      maxmom=minmom+ntypes
      mintyp=maxmom+ntypes
      nx=mintyp+ntypes
      ny=nx+lenxyz
      nz=ny+lenxyz
      call basout(z(c),z(ex),z(cont),a(ptprim),a(noprim),a(nocont),
     1            a(ptcont),nat,nprim,ntypes,nbtype,nnp,ncont,a(start),
     2            nbasis,z(zan),a(nocart),a(nobf),a(maxmom),a(mintyp),
     3            a(nx),a(ny),a(nz),a(minmom),a(index),idrop,oldnbf,
     4            ncon,prbaso,prbasn)
      lmax=intkey(card,'maximum-l-value',30,' ')
      mumax=intkey(card,'maximum-m-value',5,' ')
      rmin=fpkey(card,'minimum-r-value',1.d-04,' ')
      rmax=fpkey(card,'maximum-r-value',10.0d+00,' ')
      ksqmax=fpkey(card,'maximum-k**2',1.0d+00,' ')
      ksqmax=sqrt(ksqmax)
      nper=intkey(card,'points-per-interval-for-spline-coefficients',
     1            4,' ')
      maxlm=intkey(card,'maximum-lm-chan',1,' ')
      reschg=fpkey(card,'net-charge',0.0d+00,' ')
      nstrio=nchano*(nchano+1)/2
      nstric=nchanc*(nchanc+1)/2
      nstri=nchan*(nchan+1)/2
      rmin=rmin*ksqmax
      rmax=rmax*ksqmax
      nr=(rmax-rmin)*nper
      presnt='present'
      call iosys ('write character kohndata to kohndt',0,0,0,presnt)
      call iosys ('write integer "no. open chan" to kohndt',1,nchano,
     1             0,' ')
      call iosys ('write integer "no. closed chan" to kohndt',1,
     1            nchanc,0,' ')
      call iosys ('write integer "no. states" to kohndt',1,nchan,0,' ')
      call iosys ('write integer "max l in ylm" to kohndt',1,lmax,
     1             0,' ')
      call iosys ('write integer "max m in ylm" to kohndt',1,mumax,
     1             0,' ')
      call iosys ('write real rmin to kohndt',1,rmin,0,' ')
      call iosys ('write real rmax to kohndt',1,rmax,0,' ')
      call iosys ('write integer "total points" to kohndt',1,nr,0,' ')
      call iosys ('write integer maximum-l to kohndt',1,lmax,0,' ')
      call iosys ('write real charge to kohndt',1,reschg,0,' ')
c----------------------------------------------------------------------c
c               count up for words needed in link m6005                c
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c                there are four types of functions                     c
c            complex open, real open, closed and bound                 c
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c                 the matrices are:                                    c
c                 +                           +                        c
c        1. hpbo=h open-bound open   2. hpbc=h open-bound closed       c
c        3. hpbom=same as 1. for mo  4. same as 2. for mo              c
c                 +      +                   +      -                  c
c        5. hphp=h open-h open       6.hphm=h open-h open              c
c                 +                                                    c
c        7. hphc=h open-closed       8. hchc=closed-closed             c
c        9. hcbo=closed-bound open  10. hcbc=closed-bound closed       c
c        11.ofbo=overlap open-bound open 12. ofbom=same as 11 for mo   c
c        13.ocbc=overlap closed-bound closed                           c
c        14.ocbcm=same as 13 for mo                                    c
c----------------------------------------------------------------------c
 
      hpbow=2*maxlm*nchano*ncon*nchano
      hpbcw=2*maxlm*nchano*ncon*nchanc
      hpbowm=2*maxlm*nchano*nbf*nchano
      hpbcwm=2*maxlm*nchano*nbf*nchanc
      hphpw=2*maxlm*maxlm*nstrio
      hphmw=2*maxlm*maxlm*nchano*nchano
      hphcw=2*maxlm*maxlm*nchano*nchanc
      hchcw=maxlm*maxlm*nstric
      hcbow=maxlm*nchanc*ncon*nchano
      hcbcw=maxlm*nchanc*ncon*nchanc
      ofbow=2*maxlm*ncon*nchano
      ofbowm=2*maxlm*nbf*nchano
      ocbcw=maxlm*ncon*nchanc
      ocbcwm=maxlm*nbf*nchanc
      vec=1
      hpbo=vec+ncon*nbf
      hpbc=hpbo+hpbow
      hphp=hpbc+hpbcw
      hphm=hphp+hphpw
      hphc=hphm+hphmw
      hchc=hphc+hphcw
      hcbo=hchc+hchcw
      hcbc=hcbo+hcbow
      ofbo=hcbc+hcbcw
      ocbc=ofbo+ofbow
      hpbom=ocbc+ocbcw
      hpbcm=hpbom+hpbowm
      ofbom=hpbcm+hpbcwm
      ocbcm=ofbom+ofbowm
      scrc=ocbcm+ocbcwm
      scrr=scrc+2*maxlm*max(maxlm,ncon)
      ipow=scrr+maxlm*max(maxlm,ncon)
      x=ipow+lmax+1
      beswdc=2*nr*(lmax+1)
      if (reschg.ne.0.0) then
          beswdc=beswdc*nchano
      endif
      hs=x+nr
      hsder=hs+beswdc
      cj=hsder+beswdc
      cy=cj+beswdc
      musthv=cy+beswdc
      if(nchanc.ne.0) then
         beswdr=nr*(lmax+1)
         ms=musthv
         msder=ms+beswdr
         cms=msder+beswdr
         cmsd=cms+beswdr
         musthv=cmsd+beswdr
      endif
      noarr=4*maxlm*nchano+2*maxlm*nchanc+6+nstri
     1                    +(lmax+1)*(mumax+mumax+1)+ncon
     2                    +4*maxlm+max(4*maxlm,2*ncon)
      avail=maxcor-musthv
      bigbuf=avail/noarr
      pntbuf=min(pntbuf,bigbuf)
      write (iout,600) ncon, nbf, ncokpt, ncckpt, lmax, mumax,
     1                  maxlm, reschg
      write (iout,650) rmin, rmax
      write (iout,100) avail, bigbuf, pntbuf
      words=musthv+noarr*pntbuf
      write (iout,200) words
      if (words.gt.maxcor) then
          call lnkerr('not enough memory for calculation:quit')
      endif
      call iosys ('write integer "point buffer" to kohndt',1,pntbuf,
     1             0,' ')
    2 format (//,a80)
    3 format ( (/,5x,5(e15.8,1x) ) )
  100 format(/,5x,i8,1x,'words available',5x,'largest buffer possible',
     1       1x,i8,//,5x,'point buffer selected for now',1x,i8)
  200 format(/,5x,'need',1x,i8,1x,'words in m6005')
  500 format (/,15x,'m6010: process mesa data and memory allocation rout
     1ine')
  501 format(/,25x,'using a complex optical potential') 
  502 format (/,15x,'only one-electron integrals will be computed')
  600 format(/,5x,'no. contracted',1x,i4,2x,'no. mos',1x,i4,2x,
     1            'no. open channels',1x,i3,/,5x,
     2            'no. closed channels',1x,i3,2x,'l max',1x,i3,2x,
     3            'm max',1x,i3,2x,'maxlm/chan',1x,i3,/5x,
     4            'residual charge',1x,f8.3)
  650 format(/,5x,'minimum r value',1x,e15.8,1x,'maximum r value',
     1             1x,e15.8)        
  700 format (/,5x,'no. valence orbitals',1x,i4,2x,'no. original p-space
     1 vectors',1x,i6,/5x,'no. new p-space vectors',1x,i6)
  800 format (/,5x,'no. l**2 expansion orbitals',1x,i4,2x,'total no. l**
     12 functions',1x,i4)
  900 format (/,20x,'***** dropping basis functions *****')
 1000 format (/,20x,'***** no basis functions dropped *****')
 1500 format (/,5x,'total of nsmall and nl2 less than nbf: functions wil
     1l be deleted')
 1501 format(/,20x,'open channels will be dropped')
 1502 format(/,20x,'closed channels will be dropped')
 1600 format (/,5x,'functions kept:',(/,5x,10(i3,1x)))
 1601 format (/,5x,'open channels kept:',(/,5x,10(i3,1x)))
 1602 format (/,5x,'closed channels kept:',(/,5x,10(i3,1x)))
      call memory(-ngot,p,idum,'m6010',idum)
      call chainx(0)
      stop
      end
