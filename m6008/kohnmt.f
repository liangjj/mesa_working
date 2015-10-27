*deck @(#)kohnmt.f	1.1 9/8/91
c***begin prologue     m6008
c***date written       890524   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6008, link 6008, kohn
c***author             schneider, barry (lanl)
c***source             m6008
c***purpose            hamiltonian matrix elements for kohn
c***                   variational calculations.
c***description        bound-bound, free-free and bound-free kohn
c***                   one-body matrix elements assembled into form
c***                   needed in kohn method. input matrices in mo
c***                   representation.
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       m6008
      program kohnmt
      implicit integer (a-z)
      parameter (dimc=30 , dimmo=200, dime=100)
      real*8 z, e0, dele, eshft, echan, energy
      common a(1)
      dimension z(1)
      equivalence (z,a)
      common /io / inp, iout
      common / memory/ ioff
c----------------------------------------------------------------------c
c        the small arrays below are the only explicitly dimensioned    c
c             arrays in the code (see parameter statement)             c
c----------------------------------------------------------------------c
      dimension  nlm(dimc), ntrgt(dimc), nbscat(dimc), echan(dimc)
      dimension orblst(dimmo,dimc), finlst(dimmo,dimc), energy(dime)
      dimension nbtot(dimc), zeroc(dimc)
      logical logkey, ptest, bsym, zeroc, posinp, statex
      logical prntbb, prntbf, prntff, prntop, prntov, prntfn
      character *4096 ops
      character *1600 card
      character *24 filnm
      character *3 ans, copt
      character *8 chrkey
      character *128 filint, filkne
      character *10 cpass
      character *16 fptoc
      character *3 itoc
      call drum
      call iosys ('read character options from rwf',-1,0,0,ops)
      statex=logkey(ops,'static-exchange',.false.,' ')
      ptest=logkey(ops,'print=m6008=all',.false.,' ')
      if (ptest) then
          prntbf=.true.
          prntff=.true.
          prntbb=.true.
          prntop=.true.
          prntov=.true.
          prntfn=.true.
      else
          prntbf=logkey(ops,'print=m6008=bound-free',.false.,' ')
          prntff=logkey(ops,'print=m6008=free-free',.false.,' ')
          prntbb=logkey(ops,'print=m6008=bound-bound',.false.,' ')
          prntov=logkey(ops,'print=m6008=overlaps',.false.,' ')
          prntop=logkey(ops,'print=m6008=optical-potential',.false.,' ')
          prntfn=logkey(ops,'print=m6008=final-ints',.false.,' ')
      endif
      write (iout,1)
      filint='knints'
      if( posinp('$kohnmt',cpass) ) then
          call cardin(card)
      endif
c----------------------------------------------------------------------c
c                     some comments:                                   c
c           1. full hamiltonian matrix elements are gotten from        c
c              the optical potential file even if energy independent   c
c----------------------------------------------------------------------c
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
      call iosys ('read integer "no. l2 orbitals" from kohndt',1,
     1             nl2,0,' ')
      call iosys ('read integer "size of p space" from kohndt',1,npvec,
     1             0,' ')
      call iosys ('read integer nsmall from kohndt',1,nsmall,0,' ')
      chkit=nl2*nchan
      if (chkit.ne.npvec) then
          call lnkerr('error: npvec and nl2*nchan not equal')
      endif
      call iosys ('does "complex optical potential" exist on kohndt',
     1             -1,0,0,ans)
      fac=1
      copt='no'
      if (ans.eq.'yes') then
          copt='yes'
          fac=2
          call iosys ('write character "complex optical potential" '//
     1               'to kohnint',0,0,0,ans)
          write(iout,60)
      else
          write(iout,70)
      endif
c----------------------------------------------------------------------c
c                  channel information                                 c
c----------------------------------------------------------------------c
      maxlm=0                   
      matbf=0
      matbv=0
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
             matbv=matbv+nbscat(ch1)
             matbf=matbf+nbtot(ch1)
             write(iout,30) ch1
             write (iout,40) (orblst(i,ch1),i=1,nbscat(ch1))
             if (ntrgt(ch1).ne.0) then
                 write (iout,50) (orblst(i,ch1),i=ii,nbscat(ch1)+
     1                            ntrgt(ch1))
             endif
             zeroc(ch1)=logkey(card,'zero-channel',.false.,' ')
         endif
   10 continue
      call iosys ('write integer "total channels" to kohnint',1,ntchn,
     1            0,' ')
      call iosys ('write integer "total bound" to kohnint',1,matbv,
     1            0,' ')
c----------------------------------------------------------------------c
c             newlst makes finlst which gives finder array             c
c                for bound portion of matrix                           c
c----------------------------------------------------------------------c 
      call newlst (nbtot,nbscat,ntrgt,finlst,nchan,dimc,dimmo)
c----------------------------------------------------------------------c
c               get memory and begin calculation one energy            c
c                           at a time                                  c
c----------------------------------------------------------------------c
      call iosys ('read integer maxsiz from rwf',1,canget,0,' ')
      nstri=nchan*(nchan+1)/2
      ovbfwd=2*maxlm*nmot*nchan
      hpvbwd=ovbfwd*nchan
      hpvhpw=2*maxlm*maxlm*nstri
      hpvhmw=2*maxlm*maxlm*nchan*nchan 
      lrg=max(matbf,nmot*nmot,fac*npvec*npvec)
      words=matbf*matbf+fac*matbv*matbv+hpvhpw+hpvhmw+hpvbwd+ovbfwd+
     1      2*ntchn*ntchn+2*ntchn*ntchn+2*ntchn*matbf+ntchn*matbf+
     2      2*matbf*ntchn+matbf*ntchn+lrg+3*ntchn*matbf+ntchn*ntchn
      if (words.gt.canget) then
          call lnkerr('not enough memory for calculation:quit')
      endif
      call iosys ('write integer maxsiz to rwf',1,words,0,' ')
      call getscm(words,z,ngot,'m6008',0)
      write (iout,100) words
      hamnum=ioff
      hambb=hamnum+matbf*matbf
      hpvhp=hambb+fac*matbv*matbv 
      hpvhm=hpvhp+hpvhpw
      hpvb=hpvhm+hpvhmw
      ovbf=hpvb+hpvbwd
      hpp=ovbf+ovbfwd
      hpm=hpp+2*ntchn*ntchn
      hmm=hpm+2*ntchn*ntchn
      hpb=hmm+ntchn*ntchn
      hmb=hpb+2*ntchn*matbf
      ovpb=hmb+ntchn*matbf
      ovmb=ovpb+2*matbf*ntchn
      scrb=ovmb+ntchn*matbf
      tonow=scrb+matbf*ntchn
      hamin=tonow
      diag=tonow
      scrc=tonow+lrg
c----------------------------------------------------------------------c
c           read in the direct (non-exchange) bound-bound matrices     c
c           they are energy independent and can be calculated once.    c
c----------------------------------------------------------------------c
      call bbnum(z(hamin),z(hamnum),echan,zeroc,nchan,nmot,matbf,nbtot,
     1           orblst,finlst,dimmo,dimc,prntbb)
c----------------------------------------------------------------------c
c           write out numerator bound-bound matrix to iosys file       c
c----------------------------------------------------------------------c
      call wrnum(z(hamnum),matbf)             
c----------------------------------------------------------------------c
c                      loop over energies                              c
c----------------------------------------------------------------------c
      e0=0.d+00
      do 500 ene=1,nener
            eshft=echan(1)+.5d+00*energy(ene)
            dele=eshft-e0
            e0=eshft
c----------------------------------------------------------------------c
c           read in bound-bound hamiltonian or modify energy           c
c                         denominator                                  c
c----------------------------------------------------------------------c
            call  rdham(z(hamin),z(hambb),z(hamin),z(hambb),orblst,
     1                  nchan,nbscat,npvec,matbv,nl2,nsmall,eshft,
     2                  dimmo,dimc,prntop,statex,copt)
c----------------------------------------------------------------------c
c              write denominator matrix to iosys file                  c
c----------------------------------------------------------------------c
            call wrdeno(z(hambb),z(hambb),energy(ene),matbv,copt)        
c----------------------------------------------------------------------c
c             modify bound-bound numerator matrix                      c
c----------------------------------------------------------------------c
            call modnum(z(hamnum),z(diag),eshft,matbf)
c----------------------------------------------------------------------c
c            read in bound-free integrals                              c
c----------------------------------------------------------------------c
            filnm='ovlp-'//fptoc(energy(ene))
            call iosys ('read real '//filnm//' from kohnint',ovbfwd,
     1                  z(ovbf),0,' ')
            filnm='bf-'//fptoc(energy(ene))
            call iosys ('read real '//filnm//' from kohnint',hpvbwd,
     1                  z(hpvb),0,' ')
c----------------------------------------------------------------------c
c            read in free-free integrals                               c
c----------------------------------------------------------------------c
            filnm='ffp-'//fptoc(energy(ene))
            call iosys ('read real '//filnm//' from kohnint',hpvhpw,
     1                  z(hpvhp),0,' ')
            filnm='ffm-'//fptoc(energy(ene))
            call iosys ('read real '//filnm//' from kohnint',hpvhmw,
     1                  z(hpvhm),0,' ')
c----------------------------------------------------------------------c
c                  transform free-free integrals                       c
c----------------------------------------------------------------------c
         call finff(z(hpvhp),z(hpvhm),z(hpvb),z(ovbf),z(ovpb),
     1              z(ovmb),z(hpb),z(hmb),z(hpp),z(hpm),z(hmm),
     2              z(hamnum),z(scrb),z(scrc),nbtot,orblst,finlst,
     3              zeroc,nlm,nchan,ntchn,matbf,maxlm,nstri,nmot,dimmo,
     4              dimc,prntff,prntbf,prntov,prntfn)
         call finbf (z(hpb),z(hmb),z(ovpb),z(ovmb),z(hamnum),ntchn,
     1               matbf,prntfn)
c----------------------------------------------------------------------c
c             re-arrange bound-free integrals keeping only             c
c             the variational orbitals                                 c
c----------------------------------------------------------------------c
         call rearrg(z(hpb),z(hmb),z(ovpb),z(ovpb),matbf,matbv,ntchn)
c----------------------------------------------------------------------c
c            restore bound-bound numerator diagonal elements           c
c----------------------------------------------------------------------c
         call resnum(z(hamnum),z(diag),matbf)
c----------------------------------------------------------------------c
c            write out the matrices                                    c
c----------------------------------------------------------------------c
         filnm='hpp-'//fptoc(energy(ene))
         call iosys ('write real '//filnm//' to kohnint',2*ntchn*ntchn,
     1               z(hpp),0,' ')
         filnm='hpm-'//fptoc(energy(ene))
         call iosys ('write real '//filnm//' to kohnint',2*ntchn*ntchn,
     1               z(hpm),0,' ')
         filnm='hmm-'//fptoc(energy(ene))
         call iosys ('write real '//filnm//' to kohnint',ntchn*ntchn,
     1               z(hmm),0,' ')
         filnm='hpb-'//fptoc(energy(ene))
         call iosys ('write real '//filnm//' to kohnint',2*ntchn*matbv,
     1               z(hpb),0,' ')
         filnm='hmb-'//fptoc(energy(ene))
         call iosys ('write real '//filnm//' to kohnint',ntchn*matbv,
     1               z(hmb),0,' ')
  500 continue
      call iosys ('write integer maxsiz to rwf',1,canget,0,' ')
      call iosys ('rewind all on kohndt read-and-write',0,0,0,' ')
      call iosys ('close kohndt',0,0,0,' ')
      call iosys ('rewind all on kohnint read-and-write',0,0,0,' ')
      call iosys ('close kohnint',0,0,0,' ')
    1 format(//,20x,'***** m6008:kohn variational matrix elements *****'
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
