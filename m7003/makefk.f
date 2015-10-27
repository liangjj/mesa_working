*deck @(#)makefk.f
c***begin prologue     makefk
c***date written       920409   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m7003, link 7003
c***author             schneider, barry (nsf)
c***source             m7003
c***purpose            compute the radial integrals of the density
c***                   needed to form the two electron integrals.
c***                   
c***description        the radial functions f(k) are prepared by numerical
c***                   integration of the density with the appropriate
c***                   multipole interaction.
c***references
c
c***routines called    iosys, util and mdutil
c***end prologue       m7003
      program makefk
      implicit integer (a-z)
      parameter (dimpt=1000,dime=50,dimc=20,dimlf=20,dimmf=20,dimb=500)
      parameter (dimsym=25) 
      real*8 z, echan, energy, charge, pt(dimpt), tol, fpkey
      real *8 ecdiag
      character *80 title
      common a(1)
      dimension z(1), echan(dimc), energy(dime), nlmfre(dimc)
      dimension lfree(dimlf*dimc), mfree(dimmf*dimc), nsymb(dimsym)
      dimension nsymc(dimsym),lb(dimb), mb(dimb), lc(dimb), mc(dimb)
      dimension ecdiag(dimc*dimlf), nsymf(dimsym)
      equivalence (z,a)
      character *4096 ops
      character *128  filci
      character *3 ans
      logical logkey, nfree, prnt, den, orth
      common /io / inp, iout
      common / memory / ioff
      call drum
      title='     *****m7003:calculation of f(k) functions*****'   
      write(iout,100) title
      call iosys ('read character options from rwf',-1,0,0,ops)
      call iosys ('read character "atomic ci filename" from rwf',-1,
     1             0,0,filci)
      call iosys ('open atomci as old',0,0,0,filci)
      call iosys ('read integer "number quadrature points" from atomci',
     1             1,npts,0,' ')
      call iosys('read real "radial points" from atomci',npts,pt,0,' ')
      maxl=0
      maxblk=0
      nbndr=0
      call iosys ('does "number bound functions" exist on atomci',
     1             0,0,0,ans)   
      nprbb=0
      if (ans.ne.'no') then
          call iosys ('read integer "number bound functions" '//
     1                'from atomci',1,nbndr,0,' ')   
          call iosys ('read integer "bound symmetry list" from atomci',
     1                 dimsym,nsymb,0,' ')
          call iosys('read integer "l values of bound orbitals" from '//
     1               'atomci',nbndr,lb,0,' ')
          call iosys ('read integer "largest bound-bound pair block" '//
     1                'from atomci',1,mxprbb,0,' ')
          maxblk=max(maxblk,mxprbb*npts)
          nprbb=nbndr*(nbndr+1)/2                   
          do 1 i=1,nbndr
             maxl=max(maxl,lb(i))
    1     continue
      endif
      call iosys ('does "number complex functions" exist on '//
     1            'atomci',0,0,0,ans)   
      nprcc=0
      nprbc=0
      if (ans.ne.'no') then
          call iosys ('read integer "number complex functions" '//
     1                'from atomci',1,nbndc,0,' ')   
          call iosys ('read integer "complex symmetry list" from '//
     1                'atomci',dimsym,nsymc,0,' ')
          call iosys('read integer "l values of complex orbitals" '//
     1               'from atomci',nbndc,lc,0,' ')
          call iosys ('read integer "largest complex-complex pair '//
     1                'block" from atomci',1,mxprcc,0,' ')
          maxblk=max(maxblk,2*mxprcc*npts)          
          nprcc=nbndc*(nbndc+1)/2 
          if (nbndr.ne.0) then
              call iosys ('read integer "largest bound-complex pair '//
     1                    'block" from atomci',1,mxprbc,0,' ')
              maxblk=max(maxblk,2*mxprbc*npts)
              nprbc=nbndc*nbndr
          endif
          do 2 i=1,nbndc
             maxl=max(maxl,lc(i))
    2     continue
      endif
      maxlf=0
      sze=0
      call iosys ('does "free on" exist on atomci',0,0,0,ans)
      nprff=0
      nprbf=0
      nprcf=0
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
          call iosys ('read integer "largest free-free pair '//
     1                'block" from atomci',1,mxprff,0,' ')
          maxblk=max(maxblk,2*mxprff*npts)
          nprff=sze*(sze+1)/2
          if (nbndr.ne.0) then
              call iosys ('read integer "largest bound-free pair '//
     1                'block" from atomci',1,mxprbf,0,' ')
              maxblk=max(maxblk,2*mxprbf*npts)
              nprbf=sze*nbndr
          endif
          if (nbndc.ne.0) then
              call iosys ('read integer "largest complex-free pair '//
     1                'block" from atomci',1,mxprcf,0,' ')
              maxblk=max(maxblk,2*mxprcf*npts)
              nprcf=sze*nbndc
          endif
          cntlm=1
          do 10 i=1,nchan
              call iosys('read integer "channel l values" from '//
     1                   'atomci without rewinding',nlmfre(i),
     2                    lfree(cntlm),0,' ')
              call iosys('read integer "channel m values" from '//
     1                   'atomci without rewinding',nlmfre(i),
     2                    mfree(cntlm),0,' ')
              pnt=cntlm
              do 3 j=1,nlmfre(i)
                 maxlf=max(maxlf,lfree(pnt))
                 pnt=pnt+1
    3         continue
              cntlm=cntmlm+nlmfre(i)
   10     continue     
      endif
      prnt=logkey(ops,'print=m7003',.false.,' ')
      maxfn=max(nbndr,nbndc,sze)
      maxl=max(maxl,maxlf)
      maxk=maxl+maxl
c**********************************************************************c
c            size of largest pair block for any given density          c
c                     has already been read in above.                  c 
c**********************************************************************c
      call iosys ('read maxsiz from rwf',1,maxcor,0,' ')
      avail=maxcor-10000
      if(avail.ge.maxblk) then
         words=min(avail,maxblk)
         write(iout,3) words
      else 
         words=avail
         write(iout,4) words
      endif
      call getscm(words,z,ngot,'m7003',0)
c**********************************************************************c
c                  calculate the bound bound f(k)                      c
c**********************************************************************c
      if (nprbb.ne.0) then
          call iosys ('create real "bound-bound f(k)" on atomci',-1,0,
     1        0,' ')
          wrdcnt=0
          do 20 i=1,                 
          filsiz=npr*npts
              call iosys ('create real "bound-bound densities" on '//
     1                    'atomci',filsiz,0,0,' ')
      call iosys ('write integer maxsiz to rwf',1,maxcor,0,' ')
      call chainx(0)
    3 format(//,5x,'m7003:',/,5x,'      calculation performed in main me
     1mory. use',1x,i6,1x,'words')   
    4 format(//,5x,'m7003:',/,5x,'      calculation cannot be performed
     1in main memory. use',1x,i6,1x,'words')    
  100 format(//,a80)    
      stop
      end






