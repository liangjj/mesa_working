*deck @(#)makefn.f
c***begin prologue     makefn
c***date written       920409   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m7000, link 7000
c***author             schneider, barry (nsf)
c***source             m7000
c***purpose            driver to produce bound and free functions on
c***                   numerical grid.
c***references
c
c***routines called    iosys, util and mdutil
c***end prologue       m7000
      program makefn
      implicit integer (a-z)
      parameter (dimpt=1000,dime=50,dimc=20,dimlf=20,dimmf=20,dimsym=49)
      real*8 z, acc, alfa, rmin, rmax, rstr, rfnr, rfs, rgswtg
      real *8 echan, energy, rdel, smldel, rd26, rl, kchan, ec
      real *8 dummy(2*dimpt), fpkey, charge, maxksq
      character *80 title
      common a(1)
      dimension z(1), echan(dimc), energy(dime), rgswtg(dimpt,2)
      dimension nlmfre(dimc), lfree(dimlf,dimc), mfree(dimmf,dimc)
      dimension nsym(dimsym*dimc), types(dimsym), ltyp(dimsym)
      dimension mtyp(dimsym)
      equivalence (z,a), (rgswtg,dummy)
      character *4096 ops
      character *1600 card
      character *128  filci
      character *16 cpass, chrkey, mshtyp
      character *3 itoc
      character *8 typefn, types
      logical posinp, nfree, prrb, prcb, prfr, prgr, logkey, scale
      logical prbs
      common /io / inp, iout
      common / memory / ioff
c
      data pi /3.14159265358979323846d+00/
      call drum
      write(iout,1)
      call iosys ('read character options from rwf',-1,0,0,ops)
      prgr=logkey(ops,'print=m7000=grid',.false.,' ')
      prbs=logkey(ops,'print=m7000=basis-information',.false.,' ')      
      prrb=logkey(ops,'print=m7000=real-basis-functions',.false.,' ')
      prcb=logkey(ops,'print=m7000=complex-basis-functions',.false.,' ')
      prfr=logkey(ops,'print=m7000=free-basis-functions',.false.,' ')
      scale=logkey(ops,'m7000=scale-all-functions',.false.,' ')
      call iosys ('read character "atomic ci filename" from rwf',-1,
     1             0,0,filci)
      call iosys ('open atomci as new',0,0,0,filci)
c----------------------------------------------------------------------c
c                   basic information on basis                         c
c                          and integration points                      c
c----------------------------------------------------------------------c
      call typylm(types,ltyp,mtyp)
      if ( posinp('$functions',cpass) ) then
           call cardin(card)
      endif
      typefn=chrkey(card,'type-functions','gaussian',' ')
      nbndr=intkey(card,'number-real-functions',1,' ')
      nbndc=intkey(card,'number-complex-functions',1,' ')
      nfree=logkey(card,'free-functions',.false.,' ')
      nrg=intkey(card,'number-integration-regions',1,' ')
      charge=fpkey(card,'atomic-charge',1.d0,' ')
      call iosys ('write real charge to atomci',1,charge,0,' ')
      write(iout,2) nbndr, nbndc
      if (nfree) then
          write(iout,9)
          nn=1
          call iosys ('write integer "free on" to atomci',1,
     1                 nn,0,' ')
c----------------------------------------------------------------------c
c        variables needed for bessel function routines                 c
c----------------------------------------------------------------------c
          acc=intkey(card,'accuracy-parameter',30,' ')
          lmax=intkey(card,'maximum-l-for-bessel',10,' ')
          nper=intkey(card,'points-per-interval-for-spline-'//
     1                     'coefficients',4,' ')      
          ltop=intkey(card,'top-recursion-l',100,' ')       
          alfa=fpkey(card,'exponential-parameter',1.d+00,' ')
          subint=intkey(card,'integration-subinterval',3,' ')
          rmin=fpkey(card,'rmin',1.d-04,' ')          
          rmax=fpkey(card,'rmax',10.d+00,' ')          
          maxksq=fpkey(card,'maximum-k**2',1.d0,' ')
          rmin=rmin*sqrt(maxksq)
          rmax=rmax*sqrt(maxksq)
          rl=lmax*acc
          strtl=lmax+sqrt(rl)
          strtl=max(strtl,lmax+1)
          ltop=min(ltop,strtl)
      endif
      npts=0
      do 10 iz=1,nrg
         if ( posinp('$region-'//itoc(iz),cpass) ) then
              call cardin(card)
              mshtyp=chrkey(card,'type-quad','legendre',' ')
              npnr=intkey (card,'number-points',3,' ')
              rstr=fpkey (card,'starting-r',0.d0,' ')
              rfnr=fpkey (card,'ending-r',10.d0,' ')
              rfs=rfnr-rstr
              nps=npnr
              do 20 jz=1,nps
                 npts=npts+1
                 call lgndrx (nps,jz,rgswtg(npts,2),rgswtg(npts,1))
                 rgswtg(npts,2)=rgswtg(npts,2)*rfs
                 rgswtg(npts,1)=rgswtg(npts,1)*rfs+rstr
   20         continue
         endif
   10 continue
      if (rmin.gt.rgswtg(1,1)) then
          rmin=rgswtg(1,1)
      endif
      if (rmax.lt.rgswtg(npts,1)) then
          rmax=rgswtg(npts,1)
      endif
      call iosys ('write integer "number quadrature points" to atomci',
     1             1,npts,0,' ')
      call iosys('write real "radial points" to atomci',npts,
     1            rgswtg(1,1),0,' ')
      call iosys('write real "radial weights" to atomci',npts,
     1            rgswtg(1,2),0,' ')
      call copy(rgswtg(1,2),dummy(npts+1),npts)
      if (prgr) then
          title='radial points and weights'
          call prntrm(title,rgswtg,npts,2,npts,2,iout)
      endif
      if (scale) then
          count=npts
          do 25 i=1,npts
             count=count+1
             dummy(count)=dummy(i)*sqrt(dummy(count))
   25     continue
      endif
      wdbndr=0
      if (nbndr.ne.0) then
          wdbndr=2*nbndr*npts+4*nbndr
      endif
      wdbndc=0
      if(nbndc.ne.0) then
         wdbndc=4*nbndc*npts+5*nbndc
      endif
      wdfre=0
c----------------------------------------------------------------------c
c                 channel information and energies needed for          c
c                 generation of free functions                         c
c----------------------------------------------------------------------c
      if (nfree) then
         call izero(nsym,nchan*dimsym)
         if ( posinp('$scatt',cpass) ) then
              call cardin(card)
              nchan=intkey(card,'number-channels',1,' ')
              call fparr(card,'channel-energies',echan,nchan,' ')
              nener=intkey(card,'number-energies',1,' ')
              call fparr(card,'scattering-energies',energy,nener,' ')
              call iosys('write integer "number of channels" to atomci'
     1                    ,1,nchan,0,' ')
              call iosys('write real "channel energies" to atomci',
     1                    nchan,echan,0,' ')
              call iosys('write integer "number of scattering '//
     1                   'energies" to atomci',1,nener,0,' ')
              call iosys('write real "scattering energies" to atomci',
     1                    nener,energy,0,' ')
              write(iout,3) nchan, nener
              write(iout,4) (echan(i),i=1,nchan)
              write(iout,5) (energy(i),i=1,nener)
              sze=0
              is=1
              do 30 i=1,nchan
                 if (posinp('$chan-'//itoc(i),cpass)) then
                     call cardin(card)
                     nlmfre(i)=intkey(card,'number-lm-values',1,' ')
                     call frelm(card,nsym(is),nlmfre(i),lfree(1,i),
     1                          mfree(1,i),types,ltyp,mtyp,dimsym)
                     is=is+dimsym
                     sze=sze+nlmfre(i)
                 endif
   30         continue
              call iosys ('write integer "size of channel space" to '//
     1                    'atomci',1,sze,0,' ')           
              call iosys('write integer "number of lm pairs" to atomci',
     1                    nchan,nlmfre,0,' ')
              call iosys('create integer "channel l values" on rwf',
     1                    sze,0,0,' ')
              call iosys('create integer "channel m values" on rwf',
     1                    sze,0,0,' ')
              do 15 i=1,nchan
                 call iosys('write integer "channel l values" to '//
     1                      'atomci without rewinding',nlmfre(i),
     2                      lfree(1,i),0,' ')
                 call iosys('write integer "channel m values" to '//
     1                      'atomci without rewinding',nlmfre(i),
     2                       mfree(1,i),0,' ')
   15         continue     
         endif
         call iosys ('write integer "free symmetry list" to atomci',
     1                dimsym*nchan,nsym,0,' ')
         nr=nper*(rmax-rmin)
         nint=nr*subint
         rdel=(rmax-rmin)/nr
         smldel=(rmax-rmin)/nint
         nr=nr+1
         nint=nint+1
         rd26=rdel*rdel/6.d+00
         rmin=rmin+1.d-10
         wdfre=9*nint+8*(ltop+1)+2+4*nr+6*nr*(lmax+1)+3*npts*sze
      endif
      words=max(wdbndr,wdbndc,wdfre)
      call iosys ('read maxsiz from rwf',1,maxcor,0,' ')
      if (words.gt.maxcor) then
          call lnkerr('not enough memory for calculation:quit')
      endif
      call iosys ('write integer maxsiz to rwf',1,words,0,' ')
      call getscm(words,z,ngot,'m7000',0)
c----------------------------------------------------------------------c
c                tabulate the real bound orbitals                      c
c----------------------------------------------------------------------c
      if (nbndr.ne.0) then
          fns=ioff
          ddfns=fns+npts*nbndr
          lval=wpadti(ddfns+npts*nbndr)
          mval=lval+nbndr
          power=mval+nbndr
          alpha=iadtwp(power+nbndr)
          call rbasis(z(fns),z(ddfns),a(lval),a(mval),z(alpha),
     1                a(power),rgswtg,nsym,npts,nbndr,dimsym,
     2                prbs,typefn)
          if (prrb) then
             title='real bound functions' 
             call prntrm(title,z(fns),npts,nbndr,npts,nbndr,iout)
             title='kinetic energy operator acting on real bound '//
     1             'functions' 
             call prntrm(title,z(ddfns),npts,nbndr,npts,nbndr,iout)
          endif
          if (scale) then
              call scalfn(z(fns),z(fns),dummy(npts+1),'real',
     1                    nbndr,npts)
              call scalfn(z(ddfns),z(ddfns),dummy(npts+1),'real',
     1                    nbndr,npts)
              call nrmlz(z(fns),z(fns),z(ddfns),z(ddfns),'real',
     1                   nbndr,npts)
          endif
          call iosys('write real "bound basis functions" to atomci',
     1                npts*nbndr,z(fns),0,' ') 
          call iosys('write real "kinetic energy of bound basis '//
     1               'functions" to atomci',npts*nbndr,z(ddfns),0,' ') 
      endif       
c----------------------------------------------------------------------c
c                tabulate the complex bound orbitals                   c
c----------------------------------------------------------------------c
      if (nbndc.ne.0) then
          fns=ioff
          ddfns=fns+2*npts*nbndc
          lval=wpadti(ddfns+2*npts*nbndc)
          mval=lval+nbndc
          power=mval+nbndc
          alpha=iadtwp(power+nbndc)
          call cbasis(z(fns),z(ddfns),a(lval),a(mval),z(alpha),
     1                a(power),rgswtg,nsym,npts,nbndc,dimsym,
     2                prbs,typefn)
          if (prcb) then
              title='complex bound functions'
              call prntcmn(title,z(fns),npts,nbndc,npts,nbndc,iout,'e')
              title='kinetic energy operator acting on complex '//
     1              'bound functions'
              call prntcmn(title,z(ddfns),npts,nbndc,npts,nbndc,
     1                     iout,'e')
          endif
          if (scale) then
              call scalfn(z(fns),z(fns),dummy(npts+1),'complex',
     1                    nbndc,npts)
              call scalfn(z(ddfns),z(ddfns),dummy(npts+1),'complex',
     1                    nbndc,npts)
              call nrmlz(z(fns),z(fns),z(ddfns),z(ddfns),'complex',
     1                   nbndc,npts)
          endif
          call iosys('write real "complex basis functions" to atomci',
     1                2*npts*nbndc,z(fns),0,' ') 
          call iosys('write real "kinetic energy of complex basis '//
     1               'functions" to atomci',2*npts*nbndc,z(ddfns),0,' ') 
      endif       
c----------------------------------------------------------------------c
c             tabulate the free functions                              c
c----------------------------------------------------------------------c
      if (nfree) then
          x=ioff
          xinv=x+nint
          expfn=xinv+nint
          v=expfn+nint
          j=v+nint
          jp=j+2*(ltop+1)
          y=jp+2*(ltop+1)
          yp=y+2*(ltop+1)
          norm=yp+2*(ltop+1)
          psi0=norm+2
          psi1=psi0+nint
          driver=psi1+nint
          fun=driver+nint
          fnl=fun+2*nint
          ddfnl=fnl+2*(lmax+1)*nr
          coefc=ddfnl+(lmax+1)*nr
          coefr=coefc+2*(lmax+1)*nr
          scr=coefr+(lmax+1)*nr
c**********************************************************************c
c                  get the spline coefficients needed for              c
c                  interpolating the universal function                c
c                  at the l values and energies needed                 c
c                           for each channel                           c
c**********************************************************************c 
          write(iout,6) rmin, rmax, rdel, rd26, smldel, alfa, nint, nr,
     1                  lmax, ltop
          call genbes(z(x),z(xinv),z(expfn),z(v),z(psi0),z(psi1),
     1                z(driver),z(j),z(jp),z(y),z(yp),z(norm),z(fun),
     2                z(fnl),z(ddfnl),z(coefr),z(coefc),z(scr),rmin,
     3                rmax,rdel,smldel,alfa,nint,subint,nr,lmax,ltop)
c**********************************************************************c
c                 add arrays needed for getting universal free         c
c                 functions interpolated on real grid                  c
c**********************************************************************c
          krvec=scr+2*nr
          frefn0=krvec+npts
          ddfref=frefn0+2*npts*sze
c**********************************************************************c
c                   make required functions                            c
c**********************************************************************c
          call iosys('create real "free functions" on atomci',nener*sze*
     1                2*npts,0,0,' ')                        
          call iosys('create real "kinetic energy of free '//
     1               'functions" on atomci',nener*sze*npts,0,0,' ')
          if (prbs) then
              write(iout,7)
          endif
          do 40 ene=1,nener
             frst0=frefn0
             ddfrst=ddfref
             do 50 ch=1,nchan
                ec=energy(ene)-echan(ch)+echan(1)
                kchan=sqrt(2.d0*ec)
                do 60 lchn=1,nlmfre(ch)
                   if (prbs) then
                       write(iout,8) ch, ec, lfree(lchn,ch), 
     1                                mfree(lchn,ch)
                   endif
                   fnllc=fnl+2*nr*lfree(lchn,ch)
                   ddfnlc=ddfnl+nr*lfree(lchn,ch) 
                   coclc=coefc+2*nr*lfree(lchn,ch)
                   corlc=coefr+nr*lfree(lchn,ch)
                   call mkbes(z(fnllc),z(ddfnlc),z(coclc),z(corlc),
     1                        z(frst0),z(ddfrst),rgswtg,z(x),z(krvec),
     2                        kchan,rmin,rdel,rd26,lfree(lchn,ch),
     3                        npts,nr)
                   if (prfr) then
                       title='complex free function chan-'//itoc(ch)//
     1                       ' lval-'//itoc(lfree(lchn,ch))
                       call prntcmn(title,z(frst0),npts,1,npts,1,
     1                              iout,'e')
                       title='( H0 - E ) on complex free function chan-'
     1                       //itoc(ch)//' lval-'//itoc(lfree(lchn,ch))
                       call prntrm(title,z(ddfrst),npts,1,npts,1,
     1                              iout,'e')
                   endif
                   frst0=frst0+2*npts
                   ddfrst=ddfrst+npts
   60           continue     
   50        continue        
             if (scale) then
                 call scalfn(z(frefn0),z(frefn0),dummy(npts+1),
     1                       'complex',sze,npts)
                 call scalfn(z(ddfref),z(ddfref),dummy(npts+1),
     1                       'real',sze,npts)
             endif
             call iosys('write real "free functions" on atomci '//
     1                  'without rewinding',2*npts*sze,z(frefn0),0,' ')
             call iosys('write real "kinetic energy of free '//
     1                  'functions" on atomci without rewinding',
     2                   npts*sze,z(ddfref),0,' ')

   40     continue
      endif
      call iosys ('write integer maxsiz to rwf',1,maxcor,0,' ')
      call iosys ('rewind all on atomci read-and-write',0,0,0,' ')
      call iosys ('close atomci',0,0,0,' ')
    1 format(//,10x,'***** m7000:calculation of bound and free radial fu
     1nctions *****')   
    2 format(/,15x,'no. real bound orbitals',1x,i3,1x,'no. complex bound 
     1 orbitals',1x,i3)    
    3 format(/,15x,'no. channels',1x,i3,1x,'no. energies',1x,i3)
    4 format(/,15x,'channel energies(Hartrees)',/,(/,5x,5(e15.8,1x)))
    5 format(/,15x,'scattering energies(Hartrees)',/,(/,5x,5(e15.8,1x)))
    6 format(//,15x,'parameters for bessel function generation',//,5x,
     1           'rmin = ',e15.8,1x,'rmax = ',e15.8,1x,'rdel = ',e15.8,
     2      /,5x,'rd26 = ',e15.8,1x,'smldel = ',e15.8,1x,'alpha = ',
     3        e15.8,/,5x,'nint = ',i5,1x,'nr = ',i5,1x,'lmax = ',i3,
     4        1x,'ltop = ',i3)                   
    7 format(//,10x,'free functions',//,5x,' channel ',10x,'  energy  ',
     1           5x,' l ',5x,' m ',/)        
    8 format(7x,i3,8x,f15.8,5x,i3,5x,i3)  
    9 format(//,15x,'the basis set contains free functions')  
      call chainx(0)
      stop
      end





