*deck @(#)mrgmt.f
c***begin prologue     mrgmt
c***date written       920409   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6070, link 6060
c***author             schneider, barry (nsf)
c***source             m6070
c***purpose            merging matrix elements from ordinary Kohn
c***                   calculation with those og complex optical potential
c***references
c
c***routines called    iosys, util and mdutil
c***end prologue       m6070
      program mrgmt
      implicit integer (a-z)
      parameter (dime=100)
      common a(1)
      dimension z(1), energy(dime)
      equivalence (z,a)
      common /io / inp, iout
      common / memory / ioff
      real *8 z, energy
      logical logkey, ptest, ignore, realop
      character *4096 ops
      character *128  filkne, filint, filopt
      character *16 fptoc
      character *24 filnm
      character *80 title
      call drum
      call iosys ('read character options from rwf',-1,0,0,ops)
      ptest=logkey(ops,'print=m6070',.false.,' ')
      ignore=logkey(ops,'m6070=ignore',.false.,' ')
      realop=logkey(ops,'m6070=real-optical-potential',.false.,' ')
      call iosys ('read character "optical potential filename" '//
     1            'from rwf',-1,0,0,filopt)         
      call iosys ('read character "kohn integral filename" '//
     1            'from rwf',-1,0,0,filint)         
      call iosys ('read character "kohn data filename" '//
     1             'from rwf',0,0,0,filkne)
      write (iout,1)
      call iosys ('open kohndt as old on ssd',0,0,0,filkne)
      call iosys ('open optint as old on ssd',0,0,0,filopt)
      call iosys ('open kohnint as old on ssd',0,0,0,filint)
c----------------------------------------------------------------------c
c         read in information about energies and channels              c
c         the scattering energies are incident electron energies       c
c----------------------------------------------------------------------c
      call iosys ('read integer "no. channels" from optint',1,
     1             nchan,0,' ')
      call iosys ('read integer "no. energies" from optint',1,
     1             nener,0,' ')
      call iosys ('read real "chan energies" from optint',nchan,
     1             echan,0,' ')
      call iosys ('read real "scatt energies" from optint',nener,
     1             energy,0,' ')
      call iosys ('read integer "total no. channels" from optint',
     1             1,ntchn,0,' ')
      call iosys ('read integer "no. variational orbitals" from optint',
     1             1,matbv,0,' ')
      ntri=nchan*(nchan+1)/2
c----------------------------------------------------------------------c
c               get memory and begin calculation one energy            c
c                           at a time                                  c
c----------------------------------------------------------------------c
      call iosys ('read maxsiz from rwf',1,maxcor,0,' ')
      words=max(ntchn*ntchn,ntchn*matbv,matbv*matbv)
      words=4*words
      if (words.gt.maxcor) then
          call lnkerr('not enough memory for calculation:quit')
      endif
      call iosys ('write integer maxsiz to rwf',1,words,0,' ')
      call getscm(words,z,ngot,'m6070',0)
      vppo=ioff
      hpp=vppo+2*ntchn*ntchn
      vpmo=vppo
      hpm=hpp
      vmmo=vppo
      hmm=hpp
      vpbo=vppo
      hpb=vpbo+2*ntchn*matbv
      vmbo=vppo
      hmb=vmbo+2*ntchn*matbv
      vbbo=vppo
      hbb=vbbo+2*matbv*matbv      
      do 400 ene=1,nener
c----------------------------------------------------------------------c
         call rzero(z(vppo),2*ntchn*ntchn)
         if (.not.ignore) then
         filnm='vppo-'//fptoc(energy(ene))
         call iosys ('read real '//filnm//' from optint',2*ntchn*ntchn,
     1                z(vppo),0,' ')
         if (ptest) then
             title='vppo-'//fptoc(energy(ene))
             call prntcmn(title,z(vppo),ntchn,ntchn,ntchn,ntchn,iout,
     1                    'e')
         endif
         endif
         filnm='hpp-'//fptoc(energy(ene))
         call iosys ('read real '//filnm//' from kohnint',2*ntchn*ntchn,
     1               z(hpp),0,' ')
         if (ptest) then
             title='hpp-'//fptoc(energy(ene))
             call prntcmn(title,z(hpp),ntchn,ntchn,ntchn,ntchn,iout,'e')
         endif
         call cvadd(z(hpp),ntchn,z(vppo),ntchn,ntchn,ntchn)
         call iosys ('write real '//filnm//' to kohnint',2*ntchn*ntchn,
     1               z(hpp),0,' ')
         if (ptest) then
             title='hppc-'//fptoc(energy(ene))
             call prntcmn(title,z(hpp),ntchn,ntchn,ntchn,ntchn,iout,'e')
         endif
c----------------------------------------------------------------------c
         call rzero(z(vpmo),2*ntchn*ntchn)
         if(.not.ignore) then
         filnm='vpmo-'//fptoc(energy(ene))
         call iosys ('read real '//filnm//' from optint',2*ntchn*ntchn,
     1                z(vpmo),0,' ')
         if (ptest) then
             title='vpmo-'//fptoc(energy(ene))
             call prntcmn(title,z(vpmo),ntchn,ntchn,ntchn,ntchn,iout,
     1                    'e')
         endif
         endif
         filnm='hpm-'//fptoc(energy(ene))
         call iosys ('read real '//filnm//' from kohnint',2*ntchn*ntchn,
     1               z(hpm),0,' ')
         if (ptest) then
             title='hpm-'//fptoc(energy(ene))
             call prntcmn(title,z(hpm),ntchn,ntchn,ntchn,ntchn,iout,'e')
         endif
         call cvadd(z(hpm),ntchn,z(vpmo),ntchn,ntchn,ntchn)
         call iosys ('write real '//filnm//' to kohnint',2*ntchn*ntchn,
     1               z(hpm),0,' ')
         if (ptest) then
             title='hpmc-'//fptoc(energy(ene))
             call prntcmn(title,z(hpm),ntchn,ntchn,ntchn,ntchn,iout,'e')
         endif
c----------------------------------------------------------------------c
         call rzero(z(vmmo),2*ntchn*ntch)
         if(.not.ignore) then
         filnm='vmmo-'//fptoc(energy(ene))
         call iosys ('read real '//filnm//' from optint',2*ntchn*ntchn,
     1                z(vmmo),0,' ')
         if (ptest) then
             title='vmmo-'//fptoc(energy(ene))
             call prntcmn(title,z(vmmo),ntchn,ntchn,ntchn,ntchn,iout,
     1                    'e')
         endif
         endif
         filnm='hmm-'//fptoc(energy(ene))
         call iosys ('read real '//filnm//' from kohnint',ntchn*ntchn,
     1                z(hmm),0,' ')
         if (ptest) then
             title='hmm-'//fptoc(energy(ene))
             call prntrm(title,z(hmm),ntchn,ntchn,ntchn,ntchn,iout,
     1                   'e')
         endif
         call crvadd(z(vmmo),ntchn,z(hmm),ntchn,ntchn,ntchn) 
         call iosys ('write real '//filnm//' to kohnint',2*ntchn*ntchn,
     1                z(vmmo),0,' ')
         if (ptest) then
             title='hmmc-'//fptoc(energy(ene))
             call prntcmn(title,z(vmmo),ntchn,ntchn,ntchn,ntchn,iout,
     1                    'e')
         endif
c----------------------------------------------------------------------c
         call rzero(z(vpbo),2*ntchn*matbv)
         if (.not.ignore) then
         filnm='vpbo-'//fptoc(energy(ene))
         call iosys ('read real '//filnm//' from optint',2*ntchn*matbv,
     1                z(vpbo),0,' ')
         if (ptest) then
             title='vpbo-'//fptoc(energy(ene))
             call prntcmn(title,z(vpbo),ntchn,matbv,ntchn,matbv,iout,
     1                    'e')
         endif
         endif
         filnm='hpb-'//fptoc(energy(ene))
         call iosys ('read real '//filnm//' from kohnint',2*ntchn*matbv,
     1               z(hpb),0,' ')
         if (ptest) then
             title='hpb-'//fptoc(energy(ene))
             call prntcmn(title,z(hpb),ntchn,matbv,ntchn,matbv,iout,'e')
         endif
         call cvadd(z(hpb),ntchn,z(vpbo),ntchn,ntchn,matbv)
         call iosys ('write real '//filnm//' to kohnint',2*ntchn*matbv,
     1                z(hpb),0,' ')
         if (ptest) then
             title='hpbc-'//fptoc(energy(ene))
             call prntcmn(title,z(hpb),ntchn,matbv,ntchn,matbv,iout,'e')
         endif
c----------------------------------------------------------------------c
         call rzero(z(vmbo),2*ntchn*matbv)
         if (.not.ignore) then          
         filnm='vmbo-'//fptoc(energy(ene))
         call iosys ('read real '//filnm//' from optint',2*ntchn*matbv,
     1                z(vmbo),0,' ')
         if (ptest) then
             title='vmbo-'//fptoc(energy(ene))
             call prntcmn(title,z(vmbo),ntchn,matbv,ntchn,matbv,iout,
     1                    'e')
         endif
         endif
         filnm='hmb-'//fptoc(energy(ene))
         call iosys ('read real '//filnm//' from kohnint',
     1                ntchn*matbv,z(hmb),0,' ')
         if (ptest) then
             title='hmb-'//fptoc(energy(ene))
             call prntrm(title,z(hmb),ntchn,matbv,ntchn,matbv,iout,
     1                   'e')
         endif
         call crvadd(z(vmbo),ntchn,z(hmb),ntchn,ntchn,matbv) 
         call iosys ('write real '//filnm//' to kohnint',2*ntchn*matbv,
     1                z(vmbo),0,' ')
         if (ptest) then
             title='hmbc-'//fptoc(energy(ene))
             call prntcmn(title,z(vmbo),ntchn,matbv,ntchn,matbv,iout,
     1                    'e')
         endif
c----------------------------------------------------------------------c
         filnm='vbbo-'//fptoc(energy(ene))
         call iosys ('read real '//filnm//' from optint',2*matbv*matbv,
     1                z(vbbo),0,' ')
         if (ptest) then
             title='vbbo-'//fptoc(energy(ene))         
             call prntcmn(title,z(vbbo),matbv,matbv,matbv,matbv,iout,
     1                    'e')
         endif
         filnm='bbdn-'//fptoc(energy(ene))
         call iosys ('read real '//filnm//' from kohnint',matbv*matbv,
     1                z(hbb),0,' ')
         if (ptest) then
             title='bbdd-'//fptoc(energy(ene))         
             call prntrm(title,z(hbb),matbv,matbv,matbv,matbv,iout,
     1                   'e')
         endif
         call crvadd(z(vbbo),matbv,z(hbb),matbv,matbv,matbv)
         filnm='bbdc-'//fptoc(energy(ene))         
         call iosys ('write real '//filnm//' to kohnint',2*matbv*matbv,
     1                z(vbbo),0,' ')
         if (realop) then
             ic=vbbo
             jc=hbb
             do 10 i=1,2*matbv*matbv,2
                z(jc)=z(i+ic-1)
                jc=jc+1
   10        continue
             filnm='bbdn-'//fptoc(energy(ene))
             call iosys ('write real '//filnm//' to kohnint',
     1                    matbv*matbv,z(hbb),0,' ')
         endif     
         if (ptest) then
             title='bbdd-'//fptoc(energy(ene))         
             call prntrm(title,z(hbb),matbv,matbv,matbv,matbv,iout,
     1                   'e')
         endif
         if (ptest) then
             title='bbdc-'//fptoc(energy(ene))         
             call prntcmn(title,z(vbbo),matbv,matbv,matbv,matbv,iout,
     1                    'e')
         endif
c----------------------------------------------------------------------c
  400 continue
c----------------------------------------------------------------------c
c                   rewind and close all files                         c
c----------------------------------------------------------------------c
      call iosys ('rewind all on kohndt read-and-write',0,0,0,' ')
      call iosys ('close kohndt',0,0,0,' ')
      call iosys ('rewind all on optint read-and-write',0,0,0,' ')
      call iosys ('close optint',0,0,0,' ')
      call iosys ('rewind all on kohnint read-and-write',0,0,0,' ')
      call iosys ('close kohnint',0,0,0,' ')
      call iosys ('write integer maxsiz to rwf',1,maxcor,0,' ')
      call chainx(0)
      stop
    1 format(//,20x,'*****  m6070:merge matrix elements *****')
      end




