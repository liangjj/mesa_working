*deck @(#)meopt.f
c***begin prologue     meopt
c***date written       920409   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6060, link 6060
c***author             schneider, barry (nsf)
c***source             m6060
c***purpose            complex optical potential
c***references
c
c***routines called    iosys, util and mdutil
c***end prologue       m6060
      program meopt
      implicit integer (a-z)
      parameter (dimc=30 , dime=100 , dimbf=200, dimmo=200)
      real*8 z, echan, energy, etot
      common a(1)
      dimension z(1)
      equivalence (z,a)
      common /io / inp, iout
      common / memory / ioff
c----------------------------------------------------------------------c
c        the small arrays below are the only explicitly dimensioned    c
c             arrays in the code (see parameter statement)             c
c----------------------------------------------------------------------c
      dimension echan(dimc), energy(dime), list(dimbf)
      logical logkey, prme, propt, typ, posinp, north, rtest, zroset
      logical test
      character *4096 ops
      character *1600 card
      character *24 filnm
      character *32 xform
      character *13 grdtyp
      character *128 filgrd, filorb, filopt, filkne
      character *13 chrkey, srtfl
      character *10 cpass
      character *16 fptoc
      character *80 title
      data pi /3.14159265358979323846d+00/
c----------------------------------------------------------------------c
      call drum
      call iosys ('read character options from rwf',-1,0,0,ops)
      prme=logkey(ops,'print=m6060=matrix-elements',.false.,' ')
      propt=logkey(ops,'print=m6060=optical-potential',.false.,' ')
      rtest=logkey(ops,'real-check',.false.,' ')
      zroset=logkey(ops,'zero-potential',.false.,' ')
      north=logkey(ops,'no-orthogonality',.false.,' ')
      test=logkey(ops,'test-integral',.false.,' ')
      srtfl='"sorted orbs"'
      typ=.true.
      if( posinp('$optint',cpass) ) then
          call cardin(card)
	  nchan=intkey(card,'no-channels',1,' ')
          srtfl=chrkey(card,'orbital-file-type','con array',' ')
          if (srtfl.eq.'con array') then
              srtfl='"con array"'
              typ=.false.
          else
              srtfl='"sorted orbs"'
              typ=.true.
          endif
c----------------------------------------------------------------------c
c        read boundary condition specification , filenames etc         c
c----------------------------------------------------------------------c
      endif
      call iosys ('read character "grid filename" from rwf',-1,
     1             0,0,filgrd)
      call iosys ('read character "orbital filename" from rwf',-1,
     1             0,0,filorb)
      call iosys ('read character "kohn data filename" from rwf',-1,0,0,
     1             filkne)
      call iosys ('read character "optical potential filename" from '//
     1            'rwf',-1,0,0,filopt)
      write (iout,1)
      call iosys ('open kohndt as old',0,0,0,filkne)
      call iosys ('open grid as old',0,0,0,filgrd)
      call iosys ('open orbs as old',0,0,0,filorb)
      call iosys ('open optint as old',0,0,0,filopt)
c----------------------------------------------------------------------c
c             read grid parameters from basis set file                 c
c----------------------------------------------------------------------c
      call iosys ('read character "grid type" from orbs',0,0,0,grdtyp)
      call iosys ('read integer "no. grid pts" from orbs',1,ngrid,0,
     1            ' ')
      call iosys ('read integer "point buffer" from orbs',1,pntbuf,
     1            0,' ')
      call iosys ('read integer "no. cont" from orbs',1,ncon,0,' ')
      call iosys ('read integer "no. regions" from orbs',1,nreg,0,' ')
      call iosys ('read integer "final pts" from orbs',1,nolst,0,' ')
      nkept=ncon
      do 100 i=1,nkept
         list(i)=i
  100 continue
      if (typ) then
          write (iout,1500)
          call iosys ('read integer "no. kept" from orbs',1,nkept,
     1                0,' ')
          call iosys ('read integer "function list" from orbs',dimbf,
     1                list,0,' ')
      endif
      call iosys ('read integer "no. vlamdas" from optint',1,nolam,
     1             0,' ')
c----------------------------------------------------------------------c
c         read in information about energies and channels              c
c         the scattering energies are incident electron energies       c
c----------------------------------------------------------------------c
      call iosys ('read integer "no. open chan" from kohndt',1,nchan,
     1             0,' ')
      ntri=nchan*(nchan+1)/2
      call iosys ('read integer "no. energies" from kohndt',1,
     1             nener,0,' ')
      call iosys ('read real "chan energies" from kohndt',nchan,
     1             echan,0,' ')
      call iosys ('read real "scatt energies" from kohndt',nener,
     1             energy,0,' ')
      call iosys ('read integer "total no. mos" from kohndt',1,nmotot,
     1             0,' ')
      call iosys ('read integer "no. l2 orbitals" from kohndt',1,nl2,
     1             0,' ')
      call iosys ('read integer nsmall from kohndt',1,nsmall,0,' ')
      call iosys ('read integer "size of p space" from kohndt',1,
     1             npvec,0,' ')
      call iosys ('write character "complex optical potential" to '//
     1            'kohndt',0,0,0,'yes')
      ntot = npvec*npvec
      ntest=nchan*nl2
      if (ntest.ne.npvec) then
          call lnkerr('size of optical potential not equal to '//
     1                'nchan*nl2')
      endif
      ntest=nsmall+nl2
      if(ntest.ne.nmotot) then
         call lnkerr('nmotot and nl2+nsmall are not equal')
      endif
      if( posinp('$energy',cpass) ) then
          call cardin(card)
          call fparr(card,'channel-energies',echan,nchan,' ')
          nener=intkey(card,'no-energies',1,' ')
          call fparr(card,'scattering-energies',energy,nener,' ')
      endif
      write (iout,3) nchan,(echan(i),i=1,nchan)
c----------------------------------------------------------------------c
c               get memory and perform calculation one energy          c
c                           at a time                                  c
c----------------------------------------------------------------------c
      call iosys ('read maxsiz from rwf',1,maxcor,0,' ')
      wdlmb=2*nolam*nkept
      wdlmbm=2*nolam*nl2
      baswds=nkept*pntbuf
      wrdlst=2*nolam+2*ncon*nmotot+wdlmb+wdlmbm
      wrds1=wrdlst+2*nolam*pntbuf+4*pntbuf+nkept*pntbuf
      if (test) then
          wrds1=wrds1+2*nolam
      endif
      wrds2=wrdlst+4*ntot
      words=max(wrds1,wrds2)
      if (words.gt.maxcor) then
          call lnkerr('not enough memory for calculation:quit')
      endif
      call iosys ('write integer maxsiz to rwf',1,words,0,' ')
      call getscm(words,z,ngot,'m6060',0)
      write (iout,300) words
      eigval=ioff
      vec=eigval+2*nolam
      vectmp=vec+ncon*nmotot
      ovlmb=vectmp+ncon*nmotot
      ovlmbm=ovlmb+wdlmb
      if (test) then
         sinnt=ovlmbm+wdlmbm
         vlamda=sinnt+2*nolam
      else
         vlamda=ovlmbm+wdlmbm
      endif
      grid=vlamda+2*nolam*pntbuf      
      basis=grid+4*pntbuf
c----------------------------------------------------------------------c
c                   core reallocated here                              c
c----------------------------------------------------------------------c
      hpp=vlamda      
      vbb=hpp+2*ntot
c----------------------------------------------------------------------c
c               read then re-arrange transformation vector             c
c               after this call there are nkept ao's and nl2 mo's      c
c                     in the transformation matrix                     c
c----------------------------------------------------------------------c
      call iosys('read character "transformation vector" from kohndt',
     $           -1,0,0,xform)
      call iosys ('read real '//xform//' from kohndt',ncon*nmotot,
     1             z(vec),0,' ')
      call newvec(z(vec),z(vectmp),list,ncon,nmotot,nkept,nl2)
      call iosys ('read real "complex eigenvalues" from optint',2*nolam,
     1             z(eigval),0,' ')
      if (rtest) then
          call zimag(z(eigval),nolam,1)
      endif
      if (zroset) then
          write(iout,*) '   zeroing potential'
      endif
c----------------------------------------------------------------------c
c                      loop over energies                              c
c----------------------------------------------------------------------c
      do 400 ene=1,nener
         etot=.5d+00*energy(ene) + echan(1)
c----------------------------------------------------------------------c
c              rewind all energy-independent files                     c
c----------------------------------------------------------------------c
         call iosys ('rewind all on grid read-and-write',0,0,0,' ')
         call iosys ('rewind all on orbs read-and-write',0,0,0,' ')
         call iosys ('rewind "complex vlamdas" on optint '//
     1               'read-and-write',0,0,0,' ')
c----------------------------------------------------------------------c
c                   zero matrix elements                               c
c----------------------------------------------------------------------c
         call rzero(z(ovlmb),wdlmb)
         call rzero(z(ovlmbm),wdlmbm)
         if (test) then
             tstint=dcmplx(0.d0,0.d0)
             call rzero(z(sinnt),2*nolam)
         endif         
c----------------------------------------------------------------------c
c                   loop over grid                                     c
c----------------------------------------------------------------------c
         npnts=pntbuf
         do 430 ireg=1,nreg
            if (ireg.eq.nreg) then
                npnts=nolst
            endif
            call iosys ('read real '//grdtyp//' from grid without '//
     1                  'rewinding',4*npnts,z(grid),0,' ')
c----------------------------------------------------------------------c
c                read in a block of gaussians                          c
c----------------------------------------------------------------------c
            call iosys ('read real '//srtfl//' from orbs without '//
     1                  'rewinding',npnts*nkept,z(basis),0,' ')
c----------------------------------------------------------------------c
c              read vlamdas for this block                             c
c----------------------------------------------------------------------c
            call iosys ('read real "complex vlamdas" from optint '//
     1                  'without rewinding',2*npnts*nolam,z(vlamda),
     2                   0,' ')
            if (rtest) then
                call zimag(z(vlamda),npnts,nolam)
            endif
c----------------------------------------------------------------------c
c             multiply integration weight into vlamdas                 c
c                       and sqrt(2.) into definition of                c
c                               vlamda                                 c
c----------------------------------------------------------------------c
            call sclvlm(z(vlamda),z(grid),npnts,nolam)
c----------------------------------------------------------------------c
c               calculate overlaps with vlamdas                        c
c----------------------------------------------------------------------c
            call ovrlp(z(ovlmb),z(basis),z(vlamda),nkept,npnts,nolam)
            if (test) then
                call intgl(z(sinnt),z(vlamda),z(grid),energy(ene),
     1                     npnts,nolam)
            endif
  430    continue
c----------------------------------------------------------------------c
c              transform to molecular orbital basis                    c
c----------------------------------------------------------------------c
         if(prme) then
            call aoout(z(ovlmb),energy(ene),nolam,nkept)
         endif
         call tomo(z(vec),z(ovlmb),z(ovlmbm),nolam,nkept,nl2)
         filnm='blam-'//fptoc(energy(ene))
         call iosys ('write real '//filnm//' to optint',wdlmbm,
     1                z(ovlmbm),0,' ')
c----------------------------------------------------------------------c
c                    output section                                    c
c----------------------------------------------------------------------c
         if(prme) then
            call moout(z(ovlmbm),energy(ene),nolam,nl2)
         endif
c----------------------------------------------------------------------c
c                construct optical potential matrix elements           c
c----------------------------------------------------------------------c
         if (test) then
             call fnlint(z(sinnt),etot,z(eigval),nolam)
         endif
         call optpot(z(ovlmbm),z(vbb),etot,z(eigval),nolam,
     1               nl2,propt)
c----------------------------------------------------------------------c
c            read in hpp matrix elements from mesa                     c
c----------------------------------------------------------------------c
         filnm='h(pp)-'//fptoc(etot)
         call iosys ('read real '//filnm//' from kohndt',ntot,
     1                z(hpp),0,' ')
         if (propt) then
             title='h(pp)'
             call prntrm(title,z(hpp),nl2,nl2,nl2,nl2,iout)
         endif
c----------------------------------------------------------------------c
c               add two to form complex optical potential              c
c----------------------------------------------------------------------c
         call crvadd(z(vbb),nl2,z(hpp),nl2,nl2,nl2)
         filnm='v(pp)-'//fptoc(etot)
         call iosys ('write real '//filnm//' to kohndt',2*ntot,
     1                z(vbb),0,' ')
         call iosys ('write real '//filnm//' to optint',2*ntot,
     1                z(vbb),0,' ')
  400 continue
c----------------------------------------------------------------------c
c                   rewind and close all files                         c
c----------------------------------------------------------------------c
      call iosys ('rewind all on kohndt read-and-write',0,0,0,' ')
      call iosys ('close kohndt',0,0,0,' ')
      call iosys ('rewind all on grid read-and-write',0,0,0,' ')
      call iosys ('close grid',0,0,0,' ')
      call iosys ('rewind all on optint read-and-write',0,0,0,' ')
      call iosys ('close optint',0,0,0,' ')
      call iosys ('rewind all on orbs read-and-write',0,0,0,' ')
      call iosys ('close orbs',0,0,0,' ')
      call iosys ('write integer maxsiz to rwf',1,maxcor,0,' ')
    1 format(//,20x,'*****  m6060:optical potential matrix elements ***
     1**')
    3 format(/,5x,'target energies for ',i3,' channels:',/,(2x,5e15.8))
    4 format(/,5x,'incident energies:',/,(2x,5e15.8))
    5 format(/,5x,'spline data',/,5x,'rmin',1x,f10.5,1x,'rmax',f10.5,
     11x,'step size',1x,f10.5,/,5x,'max. l',1x,i3,1x,'inner cutoff',1x,
     2 f6.3,1x,'outer exponential cutoff',1x,f6.3,1x,'outer n cutoff',
     3 1x,i3)
    6 format(/,5x,'boundary conditions are',1x,a10)
    7 format(/,5x,'information for channel',1x,i3)
    8 format(/,5x,'l-values',(/,15x,10(i2,1x)))
    9 format(/,5x,'m-values',(/,15x,10(i2,1x)))
   10 format(/,5x,'bound-atomic-orbitals',(/,17x,10(i3,1x)))
   25 format(/,5x,'bound-molecular-orbitals',(/,17x,10(i3,1x)))
   40 format(/,5x,'expansion molecular orbitals',(/,17x,10(i3,1x)))
   50 format(/,5x,'target molecular orbitals',(/,17x,10(i3,1x)))
  300 format(/,5x,'need',1x,i8,1x,'words for calculation')
  420 format(/,5x,'incident energy(ry)',1x,f10.6,/,5x,'channel momenta',
     1       (/,10x,5(e12.5,1x)))
 1000 format(/,5x,'atomic orbital representation')
 1500 format(/,5x, 'will use short orbital list')
 2000 format(/,5x,'molecular orbital representation')
      call chainx(0)
      stop
      end














