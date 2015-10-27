*deck @(#)onemat.f	1.1 9/8/91
c***begin prologue     m6018
c***date written       930321   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6018, link 6018
c***author             schneider, barry (nsf)
c***source             m6018
c***purpose            kinetic energy free-free and bound-free matrix
c***                   elements.
c***description        free-free and bound-free kinetic energy
c***                   assembled into channel form.
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       m6018
      program onemat
      implicit integer (a-z)
      parameter (dimc=30 , dimmo=200, dime=100)
      real*8 z, echan, energy
      common a(1)
      dimension z(1)
      equivalence (z,a)
      common /io / inp, iout
      common / memory/ ioff
c----------------------------------------------------------------------c
c        the small arrays below are the only explicitly dimensioned    c
c             arrays in the code (see parameter statement)             c
c----------------------------------------------------------------------c
      dimension nlm(dimc), ntrgt(dimc), nbscat(dimc), echan(dimc)
      dimension orblst(dimmo,dimc), energy(dime)
      dimension nbtot(dimc)
      logical logkey, ptest, bsym, posinp
      logical prntbf, prntff
      character *4096 ops
      character *1600 card
      character *24 filnm
      character *3 ans
      character *128 filint, filkne
      character *10 cpass
      character *16 fptoc
      character *3 itoc
      call drum
      call iosys ('read character options from rwf',-1,0,0,ops)
      ptest=logkey(ops,'print=m6018=all',.false.,' ')
      if (ptest) then
          prntbf=.true.
          prntff=.true.
      else
          prntbf=logkey(ops,'print=m6018=bound-free',.false.,' ')
          prntff=logkey(ops,'print=m6018=free-free',.false.,' ')
      endif
      write (iout,1)
      filint='knints'
      if( posinp('$onemat',cpass) ) then
          call cardin(card)
      endif
      call iosys ('read character "kohn data filename" from rwf',-1,0,0,
     1             filkne)
      call iosys ('open kohndt as old',0,0,0,filkne)
      call iosys('read character "kohn integral filename" from rwf',
     1            -1,0,0,filint)
      call iosys ('open kohnint as old',0,0,0,filint)
      call iosys ('read character symmetry from kohnint',-1,0,0,ans)
      if (ans.eq.'on') then
          bsym=.true.
      else
          bsym=.false.
      endif
      call iosys ('read integer "no. energies" from kohnint',1,
     1            nener,0,' ')
      call iosys ('read real "scatt energies" from kohnint',nener,
     1             energy,0,' ')
      call iosys ('read integer "no. channels" from kohnint',1,
     1            nchan,0,' ')
      call iosys ('read integer "lm-vals-chan" from kohnint',nchan,
     1            nlm,0,' ')
      call iosys ('read real "chan energies" from kohnint',nchan,echan,
     1            0,' ')
      call iosys ('read integer "no. contracted" from kohndt',1,
     1             naot,0,' ')
      call iosys ('read integer "total no. mos" from kohndt',1,
     1             nmot,0,' ')
c----------------------------------------------------------------------c
c                  channel information                                 c
c----------------------------------------------------------------------c
      bbsiz=0
      bfsiz=0
      maxlm=0                   
      maxb=0
      ntchn=0
      do 10 ch1=1,nchan
         if( posinp ('$chan-'//itoc(ch1),cpass) ) then
             call cardin(card)
             maxlm=max(maxlm,nlm(ch1))
             ntchn=ntchn+nlm(ch1)
             if (.not.bsym) then
                 nbscat(ch1)=intkey(card,'no-expansion-mos',1,' ')
                 call intarr(card,'expansion-mos',orblst(1,ch1),
     1                       nbscat(ch1),' ')
                 ntrgt(ch1)=intkey(card,'no-target-mos',1,' ')
                 ii=nbscat(ch1)+1
                 call intarr(card,'target-mos',orblst(ii,ch1),
     1                       ntrgt(ch1),' ')
             else
                 nbscat(ch1)=nmot
                 ntrgt(ch1)=0
                 do 20 orb=1,nmot
                    orblst(orb,ch1)=orb
   20            continue
             endif
             nbtot(ch1)=nbscat(ch1)+ntrgt(ch1)
             maxb=max(maxb,nbtot(ch1))
             bbsiz=bbsiz+nbtot(ch1)*nbtot(ch1)
             bfsiz=bfsiz+nbtot(ch1)*nlm(ch1)
             write(iout,30) ch1
             write (iout,40) (orblst(i,ch1),i=1,nbscat(ch1))
             if (ntrgt(ch1).ne.0) then
                 write (iout,50) (orblst(i,ch1),i=ii,nbtot(ch1))
             endif
         endif
   10 continue
      call iosys ('write integer "total channels" to kohnint',1,ntchn,
     1            0,' ')
c----------------------------------------------------------------------c
c               get memory and begin calculation one energy            c
c                           at a time                                  c
c----------------------------------------------------------------------c
      call iosys ('read integer maxsiz from rwf',1,canget,0,' ')
      ovbfwd=2*maxlm*nmot*nchan
      hpvbwd=ovbfwd
      hpvhpw=2*maxlm*nchan
      ovppw=hpvhpw
      ovpmw=ovppw
      words=nmot*nmot+bbsiz+maxb+maxb+ovbfwd+hpvbwd+6*bfsiz+hpvhpw+
     1      ovppw+ovpmw+7*ntchn
      if (words.gt.canget) then
          call lnkerr('not enough memory for calculation:quit')
      endif
      call iosys ('write integer maxsiz to rwf',1,words,0,' ')
      call getscm(words,z,ngot,'m6018',0)
      write (iout,100) words
      kein=ioff
      keout=kein+nmot*nmot
      eig=keout+bbsiz
      dum=eig+maxb
      ovbf=dum+maxb 
      hpvb=ovbf+ovbfwd
      opb=hpvb+hpvbwd
      hpb=opb+bfsiz+bfsiz
      omb=hpb+bfsiz+bfsiz
      hmb=omb+bfsiz
      hpvhp=hmb+bfsiz
      ovpp=hpvhp+hpvhpw
      ovpm=ovpp+ovppw
      opp=ovpm+ovpmw
      opm=opp+ntchn+ntchn
      omm=opm+ntchn+ntchn
      hpp=omm+ntchn
      words=hpp+ntchn+ntchn
c----------------------------------------------------------------------c
c           read in the kinetic energy bound-bound matrices            c
c           they are energy independent and can be calculated once.    c
c           also diagonalize the kinetic energy matrix and put out     c
c           the eigenvalues and transformation matrix for each         c
c                               channel                                c 
c----------------------------------------------------------------------c
      call bbkin(z(kein),z(keout),z(eig),z(dum),echan,nchan,nmot,
     1           nbtot,orblst,dimmo,dimc,prntbb)
c----------------------------------------------------------------------c
c                      loop over energies                              c
c----------------------------------------------------------------------c
      do 500 ene=1,nener
c----------------------------------------------------------------------c
c            read in bound-free integrals                              c
c----------------------------------------------------------------------c
         filnm='ovlp-'//fptoc(energy(ene))
         call iosys ('read real '//filnm//' from kohnint',ovbfwd,
     1               z(ovbf),0,' ')
         filnm='bf-'//fptoc(energy(ene))
         call iosys ('read real '//filnm//' from kohnint',hpvbwd,
     1               z(hpvb),0,' ')
c----------------------------------------------------------------------c
c          transform to the bound basis diagonalizing the kinetic      c
c                       energy operator                                c
c----------------------------------------------------------------------c
         call tonwbf (z(ovbf),z(hpvb),z(ovbf),z(hpvb),z(keout),
     1                z(opb),z(hpb),nbtot,orblst,nlm,nchan,
     2                maxlm,nmot,dimmo,dimc)
         call finbf (z(ovbf),z(hpvb),z(opb),z(omb),z(hpb),z(hmb),
     1               nlm,nbtot,nchan,dimc,words,prntbf)
c----------------------------------------------------------------------c
c                   write out the bound-free integrals                 c
c----------------------------------------------------------------------c      
         filnm='opb-'//fptoc(energy(ene))
         call iosys ('write real '//filnm//' to kohnint',2*words,
     1                z(opb),0,' ')
         filnm='hpb-'//fptoc(energy(ene))
         call iosys ('write real '//filnm//' to kohnint',2*words,
     1                z(hpb),0,' ')
         filnm='omb-'//fptoc(energy(ene))
         call iosys ('write real '//filnm//' to kohnint',words,
     1                z(omb),0,' ')
         filnm='hmb-'//fptoc(energy(ene))
         call iosys ('write real '//filnm//' to kohnint',words,
     1                z(hmb),0,' ')
c----------------------------------------------------------------------c
c                  read in free-free integrals                         c
c----------------------------------------------------------------------c
         filnm='ffp-'//fptoc(energy(ene))
         call iosys ('read real '//filnm//' from kohnint',hpvhpw,
     1               z(hpvhp),0,' ')
         filnm='ppo-'//fptoc(energy(ene))
         call iosys ('read real '//filnm//' from kohnint',ovppw,
     1               z(ovpp),0,' ')
         filnm='pmo-'//fptoc(energy(ene))
         call iosys ('read real '//filnm//' from kohnint',ovpmw,
     1               z(ovpm),0,' ')
c----------------------------------------------------------------------c
c                  reformat free-free integrals                        c
c----------------------------------------------------------------------c
         call finff(z(ovpp),z(ovpm),z(hpvhp),z(opp),z(opm),z(omm),
     1              z(hpp),nlm,nchan,maxlm,dimc,words,prntff)
c----------------------------------------------------------------------c
c            write out the matrices                                    c
c----------------------------------------------------------------------c
         filnm='opp-'//fptoc(energy(ene))
         call iosys ('write real '//filnm//' to kohnint',2*words,
     1                z(opp),0,' ')
         filnm='opm-'//fptoc(energy(ene))
         call iosys ('write real '//filnm//' to kohnint',2*words,
     1                z(opm),0,' ')
         filnm='omm-'//fptoc(energy(ene))
         call iosys ('write real '//filnm//' to kohnint',words,
     1                z(omm),0,' ')
         filnm='hpp-'//fptoc(energy(ene))
         call iosys ('write real '//filnm//' to kohnint',2*words,
     1                z(hpp),0,' ')
  500 continue
      call iosys ('write integer maxsiz to rwf',1,canget,0,' ')
      call iosys ('rewind all on kohndt read-and-write',0,0,0,' ')
      call iosys ('close kohndt',0,0,0,' ')
      call iosys ('rewind all on kohnint read-and-write',0,0,0,' ')
      call iosys ('close kohnint',0,0,0,' ')
    1 format(//,20x,'***** m6018:one electron matrix elements *****'
     1      )
   30 format(/,5x,'information for channel',1x,i3)
   40 format(/,5x,'expansion molecular orbitals',(/,17x,10(i3,1x)))
   50 format(/,5x,'target molecular orbitals',(/,17x,10(i3,1x)))
   60 format(//,20x,'using a complex optical potential')
   70 format(//,20x,'using a real optical potential')  
   11 format(a80)
  100 format(/,5x,'need',1x,i8,1x,'words for calculation')
      call chainx(0)
      stop
      end
