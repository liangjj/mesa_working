*deck @(#)mkonel.f
c***begin prologue     mkonel
c***date written       920409   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m7001, link 7001
c***author             schneider, barry (nsf)
c***source             m7001
c***purpose            driver for one electron atomic radial integrals
c***                   and density matrices.
c***                   
c***references
c
c***routines called    iosys, util and mdutil
c***end prologue       m7001
      program mkonel
      implicit integer (a-z)
      parameter (dimpt=1000,dime=50,dimc=20,dimlf=20,dimmf=20,dimb=500)
      parameter (dimsym=49) 
      real*8 z, echan, energy, charge, pt(dimpt), tol, fpkey
      real *8 ecdiag
      character *80 title
      common a(1)
      dimension z(1), echan(dimc), energy(dime), nlmfre(dimc)
      dimension lfree(dimlf*dimc), mfree(dimmf*dimc), nsymb(dimsym)
      dimension lb(dimb), mb(dimb), lc(dimb), mc(dimb), nsymc(dimsym)
      dimension nsymf(dimsym,dimc), ecdiag(dimc)
      equivalence (z,a)
      character *4096 ops
      character *128  filci
      character *3 ans
      logical logkey, nfree, prnt, den, orth, tstrho
      common /io / inp, iout
      common / memory / ioff
      call drum
      title='     ***** m7001:calculation of bound and free one '//
     1                       'electron integrals *****'   
      write(iout,100) title
      call iosys ('read character options from rwf',-1,0,0,ops)
      den=logkey(ops,'m7001=no-densities',.false.,' ')
      orth=logkey(ops,'m7001=schmidt',.false.,' ')
      tol=fpkey(ops,'overlap-tolerance',1.d-06,' ')
      tstrho=logkey(ops,'m7001=test-densities',.false.,' ')
      call iosys ('read character "atomic ci filename" from rwf',-1,
     1             0,0,filci)
c**********************************************************************c
c              this section reads in the basic info on the             c
c              orbitals, quadrature points, energies etc.              c
c**********************************************************************c
      call iosys ('open atomci as old',0,0,0,filci)
      call iosys ('read real charge from atomci',1,charge,0,' ')
      call iosys ('read integer "number quadrature points" from atomci',
     1             1,npts,0,' ')
      call iosys('read real "radial points" from atomci',npts,pt,0,' ')
      nbndr=0
      call iosys ('does "number bound functions" exist on atomci',
     1             0,0,0,ans)   
      if (ans.ne.'no') then
          call iosys ('read integer "number bound functions" '//
     1                'from atomci',1,nbndr,0,' ')   
      endif
      nbndc=0
      call iosys ('does "number complex functions" exist on '//
     1            'atomci',0,0,0,ans)   
      if (ans.ne.'no') then
          call iosys ('read integer "number complex functions" '//
     1                'from atomci',1,nbndc,0,' ')   
      endif
      sze=0
      call iosys ('does "free on" exist on atomci',0,0,0,ans)
      if (ans.ne.'no') then
          nfree=.true.
          call iosys('read integer "number of channels" from atomci',
     1                1,nchan,0,' ')
          call iosys('read integer "number of scattering '//
     1               'energies" from atomci',1,nener,0,' ')
          call iosys('read integer "number of lm pairs" from atomci',
     1                nchan,nlmfre,0,' ')
          call iosys ('read integer "size of channel space" from '//
     1                'atomci',1,sze,0,' ')           
          cntlm=1
          maxfrl=0   
          do 10 i=1,nchan
              call iosys('read integer "channel l values" from '//
     1                   'atomci without rewinding',nlmfre(i),
     2                    lfree(cntlm),0,' ')
              call iosys('read integer "channel m values" from '//
     1                   'atomci without rewinding',nlmfre(i),
     2                    mfree(cntlm),0,' ')
              pnt=cntlm
              do 15 j=1,nlmfre(i)
                 maxfrl=max(maxfrl,lfree(pnt))
                 pnt=pnt+1
   15         continue
              cntlm=cntlm+nlmfre(i)
   10     continue     
      endif
      prnt=logkey(ops,'print=m7001',.false.,' ')
c**********************************************************************c
c                 this section handles the memory allocation           c
c**********************************************************************c
      maxfn=max(nbndr,nbndc,sze)
      maxsiz=max(maxfn,npts)
      maxarr=max(maxfn*maxfn,maxfn*maxsiz)
      words=0
      wcalc=0
      if (nbndr.ne.0) then
          wcalc=wcalc+nbndr*npts
          words=words+2*nbndr*npts
      endif
      if(nbndc.ne.0) then
         wcalc=wcalc+2*nbndc*npts
         words=words+4*nbndc*npts
      endif
      if(nfree) then
         wcalc=wcalc+4*sze*npts
         words=words+8*sze*npts
      endif
      wrds1=words+9*maxarr+8*maxfn
      call iosys ('read maxsiz from rwf',1,maxcor,0,' ')
      if (.not.den) then
          left=maxcor-wcalc-1000-maxarr*2
          npair=left/npts
          npair=min(npair,maxfn*maxfn)
          wrds2=wcalc+npair*npts
          write(iout,7) npair
      endif            
      words=max(wrds1,wrds2)
      if (words.gt.maxcor) then
          call lnkerr('not enough memory for calculation:quit')
      endif
      call iosys ('write integer maxsiz to rwf',1,words,0,' ')
      call getscm(words,z,ngot,'m7001',0)
      fns=ioff
      fnsc=fns+nbndr*npts
      frefn0=fnsc+2*nbndc*npts
      frefn1=frefn0+2*npts*sze
      ddfns=frefn1+2*npts*sze
      ddfnsc=ddfns+nbndr*npts
      ddfre0=ddfnsc+2*npts*nbndc
      ddfre1=ddfre0+2*npts*sze
      msk=ddfre1+2*npts*sze
      s=msk+maxarr
      t=s+2*maxarr
      v=t+2*maxarr
      scr=v+2*maxarr
      eig=scr+2*maxarr
      work=eig+2*maxfn
      val=ddfns
      rho=ddfns+maxarr*2
c**********************************************************************c
c              read in bound type functions and do all bound           c
c                            type integrals                            c
c                  these are energy independent                        c
c**********************************************************************c
      if(nbndr.ne.0) then
         call iosys('read integer "l values of bound orbitals" '//
     1               'from atomci',nbndr,lb,0,' ')
         call iosys('read integer "m values of bound orbitals" '//
     1              'from atomci',nbndr,mb,0,' ')
         call iosys('read real "bound basis functions" from atomci',
     1               npts*nbndr,z(fns),0,' ') 
         call iosys('read real "kinetic energy of bound basis '//
     1              'functions" from atomci',npts*nbndr,z(ddfns),0,' ') 
         call mskone(z(msk),lb,mb,nbndr,lb,mb,nbndr)
c**********************************************************************c
c              transform bound functions to an orthonormal basis.      c
c              then compute all one electron integrals in this         c
c                               basis.                                 c
c**********************************************************************c
         call onebb(z(fns),z(ddfns),lb,mb,z(s),z(t),z(v),z(scr),z(msk),
     1              z(eig),z(work),pt,charge,tol,nsymb,npts,
     2              nbndr,noutb,dimsym,prnt,orth)
         write(iout,1)
      endif
c**********************************************************************c
c              read in complex type functions and do all complex       c
c              complex type integrals. these are energy independent    c
c**********************************************************************c
      if(nbndc.ne.0) then
         call iosys('read integer "l values of complex orbitals" '//
     1               'from atomci',nbndc,lc,0,' ')
         call iosys('read integer "m values of complex orbitals" '//
     1              'from atomci',nbndc,mc,0,' ')
         call iosys('read real "complex basis functions" from atomci',
     1               2*npts*nbndc,z(fnsc),0,' ') 
         call iosys('read real "kinetic energy of complex basis '//
     1              'functions" from atomci',2*npts*nbndc,z(ddfnsc),
     1               0,' ')
c**********************************************************************c
c               transform complex functions to basis orthogonal        c
c                           to bound functions.                        c
c**********************************************************************c 
         call tomocc(z(fns),z(ddfns),z(fnsc),z(ddfnsc),z(scr),lb,mb,lc,
     1               mc,z(msk),z(s),z(s),z(eig),z(v),z(work),z(work),
     2               tol,nsymc,npts,noutb,nbndc,noutc,dimsym,prnt,orth)
c**********************************************************************c
c               compute bound-complex integrals in new basis           c
c**********************************************************************c
         if (noutc.ne.0) then
             if (noutb.ne.0) then
                 call mskone(z(msk),lb,mb,noutb,lc,mc,noutc)
                 call onebc(z(fns),z(fnsc),z(ddfnsc),z(s),z(t),z(v),
     1                      z(scr),z(msk),pt,charge,npts,noutb,noutc,
     2                      prnt)
             endif
             if (prnt) then
                 if(noutb.ne.0) then
                    write(iout,2)
                 endif
             endif
c**********************************************************************c
c            compute complex integrals in new basis                    c
c**********************************************************************c
             call mskone(z(msk),lc,mc,noutc,lc,mc,noutc)
             call onecc(z(fnsc),z(ddfnsc),z(s),z(t),z(v),z(scr),z(msk),
     1                  pt,charge,npts,noutc,prnt)
             if (prnt) then      
                 write(iout,3)
             endif
         endif
      endif
c**********************************************************************c
c             integrals involving the free functions                   c
c             these are energy dependent                               c
c**********************************************************************c
      if (nfree) then
          call iosys ('create real "free molecular basis functions" '//
     1                'on atomci',4*npts*sze*nener,0,0,' ')
          call iosys ('create real "kinetic energy of free '//
     1                'molecular basis functions" on atomci',
     1                4*npts*sze*nener,0,0,' ')
          call iosys ('read integer "free symmetry list" from atomci',
     1                 dimsym*nchan,nsymf,0,' ')
          call iosys ('write integer "free0 orthogonal symmetry '//
     1                'list" to atomci',dimsym*nchan,nsymf,0,' ')
          call iosys ('write integer "free1 orthogonal symmetry '//
     1                'list" to atomci',dimsym*nchan,nsymf,0,' ')
          do 30 ene=1,nener
             call iosys ('read real "free functions" from atomci '//
     1                   'without rewinding',2*npts*sze,z(frefn0),0,' ')
             call iosys ('read real "kinetic energy of free '//
     1                   'functions" from atomci without rewinding',
     2                    npts*sze,z(ddfre0),0,' ')
             do 40 ch=1,nchan
                ecdiag(ch)=energy(ene)-echan(ch)+echan(1)
   40        continue     
c**********************************************************************c
c               transform free functions and ( H0 - E ) to basis       c
c               orthogonal to bound and complex functions.             c
c**********************************************************************c 
             call mkfrs1(z(frefn0),z(frefn1),npts,sze)
             call rzero(z(ddfre1),2*npts*sze)
             call tomoff(z(frefn0),z(frefn1),z(ddfre0),z(ddfre1),
     1                   z(fns),z(ddfns),z(fnsc),z(ddfnsc),z(s),z(t),
     2                   z(v),z(scr),z(t),z(v),z(scr),z(msk),ecdiag,
     3                   lfree,mfree,lb,mb,lc,mc,nlmfre,nchan,npts,
     4                   sze,noutb,noutc,prnt)
             call iosys ('write real "free molecular basis '//
     1                   'functions" to atomci without rewinding',
     2                   2*npts*sze,z(frefn0),0,' ')
             call iosys ('write real "free molecular basis '//
     1                   'functions" to atomci without rewinding',
     2                   2*npts*sze,z(frefn1),0,' ')
             call iosys ('write real "kinetic energy of free '//
     1                   'molecular basis functions" to atomci '//
     2                   'without rewinding',2*npts*sze,z(ddfre0),
     3                    0,' ')
             call iosys ('write real "kinetic energy of free '//
     1                   'molecular basis functions" to atomci '//
     2                   'without rewinding',2*npts*sze,z(ddfre1),
     3                    0,' ')
             call mskone(z(msk),lfree,mfree,sze,lfree,mfree,sze)
             call oneff(z(frefn0),z(frefn1),z(ddfre0),z(ddfre1),z(s),
     1                  z(t),z(v),z(scr),z(msk),pt,charge,npts,sze,prnt)
             write(iout,4)
             if (noutb.ne.0) then
                 call mskone(z(msk),lb,mb,noutb,lfree,mfree,sze)
                 call onebf(z(fns),z(frefn0),z(frefn1),z(ddfre0),
     1                      z(ddfre1),z(s),z(t),z(v),z(scr),z(msk),pt,
     2                      charge,npts,noutb,sze,prnt)
                 write(iout,5)
             endif
             if (noutc.ne.0) then
                 call mskone(z(msk),lc,mc,noutc,lfree,mfree,sze)
                 call onecf(z(fnsc),z(frefn0),z(frefn1),z(ddfre0),
     1                      z(ddfre1),z(s),z(t),z(v),z(scr),z(msk),
     2                      pt,charge,npts,noutc,sze,prnt)
                 write(iout,6)
             endif
   30     continue
      endif
      write(iout,*)
      write(iout,*)
      write (iout,*) '                    all one electron '//
     1               'integrals finished'
      write(iout,*)
      write(iout,*)
c**********************************************************************c
c                all one-electron integrals are now calculated         c
c                in an orthogonal basis and out on disk. now           c
c                           compute the densities                      c
c**********************************************************************c
      if (.not.den) then
          write(iout,*)
          write(iout,*) '                    density matrix '//
     1                  'calculation requested'
          write(iout,*)
          write(iout,*)
c**********************************************************************c
c                       bound-bound density                            c
c**********************************************************************c
          if (noutb.ne.0) then
              call bbden(z(fns),z(rho),z(val),nsymb,npts,noutb,
     1                   npair,dimsym,tstrho)
c**********************************************************************c
c                       bound-complex density                          c
c**********************************************************************c
              if (noutc.ne.0) then
                  call bcden(z(fns),z(fnsc),z(rho),z(val),nsymb,nsymc,
     1                       npts,noutb,noutc,npair,dimsym,tstrho)
              endif
c**********************************************************************c
c                       bound-free density                             c
c**********************************************************************c
              if (nfree) then
                  call bfden(z(fns),z(frefn0),z(rho),z(val),nsymb,
     1                       nsymf,npts,noutb,nchan,sze,npair,
     2                       dimsym,'free0',tstrho)
                  call bfden(z(fns),z(frefn1),z(rho),z(val),nsymb,
     1                       nsymf,npts,noutb,nchan,sze,npair,
     2                       dimsym,'free1',tstrho)
              endif
          endif
c**********************************************************************c
c                       complex-complex density                        c
c**********************************************************************c
          if (noutc.ne.0) then
              call ccden(z(fnsc),z(rho),z(val),nsymc,npts,noutc,
     1                   npair,dimsym,tstrho)
c**********************************************************************c
c                       complex-free density                           c
c**********************************************************************c
              if (nfree) then
                  call cfden(z(fnsc),z(frefn0),z(rho),z(val),nsymc,
     1                       nsymf,npts,noutc,nchan,sze,npair,
     2                       dimsym,'free0',tstrho)
                  call cfden(z(fnsc),z(frefn1),z(rho),z(val),nsymc,
     1                       nsymf,npts,noutc,nchan,sze,npair,
     2                       dimsym,'free1',tstrho)
              endif
          endif
          if (nfree) then
              call ffden(z(frefn0),z(frefn0),z(rho),z(val),nsymf,
     1                   nlmfre,npts,nchan,sze,npair,dimsym,'free0',
     2                   'free0',tstrho)
              call ffden(z(frefn1),z(frefn0),z(rho),z(val),nsymf,
     1                   nlmfre,npts,nchan,sze,npair,dimsym,'free1',
     2                   'free0',tstrho)
              call ffden(z(frefn1),z(frefn1),z(rho),z(val),nsymf,
     1                   nlmfre,npts,nchan,sze,npair,dimsym,'free1',
     2                   'free1',tstrho)
          endif
      endif
    1 format(//,'*bound-bound one electron integrals finished*')    
    2 format(//,'*bound-complex one electron integrals finished*')       
    3 format(//,'*complex-complex one electron integrals finished*')       
    4 format(//,'*free-free one electron integrals finished*')       
    5 format(//,'*bound-free one electron integrals finished*')       
    6 format(//,'*complex-free one electron integrals finished*')       
    7 format(//,'*calculating densities*',/,' can store',1x,i5,1x,
     1            'in main storage',/)
      call iosys ('rewind all on atomci read-and-write',0,0,0,' ')
      call iosys ('close atomci',0,0,0,' ')
      call iosys ('write integer maxsiz to rwf',1,maxcor,0,' ')
      call chainx(0)
  100 format(//,a80)    
      stop
      end






