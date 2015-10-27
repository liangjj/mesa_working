c $Header: ylm.f,v 1.4 92/12/31 14:45:59 bis Exp $
*deck ylm.f
c***begin prologue     m6201
c***date written       921219   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6201, link 6201, spherical harmonics
c***author             schneider, b. i.(lanl/nsf)
c***source             m6201
c***purpose            spherical harmonic generation
c***description        calculates spherical harmonics on grid
c***                   for integral calculation. fully vectorized
c*** 
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       m6201
      program ylm 
      implicit integer (a-z)
      parameter ( dimshl=50 )
      character*2 chvr1, chvr2, itoc
      character*8 cpass, chrkey
      character*128 fillam
      character*1600 card
      character*4096 ops
      character*80 title
      logical logkey, prnt, prntlm, posinp, tstnrm, tstint
      real*8 z, rms
      dimension z(1)
      dimension nr(dimshl), nth(dimshl), nph(dimshl), nsub(dimshl)
      common a(1)
      common /io/ inp, iout
      common /memory / ioff
      equivalence (z,a)
      call drum
      call iosys ('read character options from rwf',-1,0,0,ops)
      prnt=logkey(ops,'print=m6201=ylm',.false.,' ')
      prntlm=logkey(ops,'print=m6201=lm-decomposition',.false.,' ')
      tstnrm=logkey(ops,'m6201=test-normalization',.false.,' ')
      tstint=logkey(ops,'m6201=test-yukawa-integral',.false.,' ')
      write (iout,100)
      mumax=0
      lmax=20
      if ( posinp ('$ylm',cpass) ) then
           call cardin(card)
           lmax=intkey(card,'maximum-l-value',20,' ')
           mumax=intkey(card,'maximum-m-value',0,' ')
      endif
      write (iout,300) lmax, mumax
      call iosys ('read character "linear algebraic filename" from rwf',
     1                 -1,0,0,fillam)
      call iosys ('open lamdat as old',0,0,0,fillam)
      call iosys ('read integer "point buffer" from lamdat',1,
     1            pntbuf,0,' ')
      call iosys ('read integer "number of shells" from lamdat',
     1             1,nshell,0,' ')
      call iosys ('read integer "number radial points per shell" '//
     1            'from lamdat',nshell,nr,0,' ')
      call iosys ('read integer "number theta points per shell" '//
     1            'from lamdat',nshell,nth,0,' ')
      call iosys ('read integer "number phi points per shell" '//
     1            'from lamdat',nshell,nph,0,' ')
      write (iout,400) pntbuf
      call iosys ('write integer "max l in ylm" to lamdat',1,
     1             lmax,0,' ')
      call iosys ('write integer "max m in ylm" to lamdat',1,
     1             mumax,0,' ')
      maxnr=0
      maxth=0
      maxph=0
      do 10 ns=1,nshell
         maxnr=max(maxnr,nr(ns))
         maxth=max(maxth,nth(ns))
         maxph=max(maxph,nph(ns))
   10 continue 
      lmmax=lmax+mumax
      maxfac=lmmax+10
c----------------------------------------------------------------------c
c                     memory allocation                                c
c----------------------------------------------------------------------c
      lplus=lmax+1
      lsq=lplus*lplus
      lpr=2*maxnr*maxth*maxph
      lhold=2*maxth
      lpr1=maxnr*lhold
      lpr2=2*maxnr*lplus
      lpr3=max(lsq,4,lpr1)
      ioff=1
      do 20 i=1,2
          dfct=ioff
          ddfct=dfct+maxfac
          theta=ddfct+maxfac
          plm=theta+maxth
          wtthe=plm+maxth*(lmax+1)
          cphi=wtthe+maxth+maxth
          sphi=cphi+maxph
          phi=sphi+maxph
          phifn=phi+maxph
          wtphi=phifn+(maxph+maxph)
          words=wtphi+maxph+maxph
          if (tstnrm.and..not.tstint) then
              scr=words
              words=scr+max(lsq,4)
          elseif (.not.tstnrm.and.tstint) then
              fio=words
              scr=fio+lpr
              flm=scr+lpr1
              words=flm+lpr2
          elseif (tstnrm.and.tstint) then
              fio=words
              flm=fio+lpr
              scr=flm+lpr2
              words=scr+lpr3
          endif         
          words=wpadti(words)
          if (i.eq.1) then
              call iosys ('read integer maxsiz from rwf',1,
     1                     canget,0,' ')
              if (words.gt.canget) then
                  call lnkerr('not enough memory. will quit')
              endif
              call iosys ('write integer maxsiz to rwf',1,
     1                     words,0,' ')
              call getscm(words,z,ngot,'m6201',0)      
          endif
   20 continue
      call fact(z(dfct),z(ddfct),maxfac)
c----------------------------------------------------------------------c
c                    loop over shells                                  c
c----------------------------------------------------------------------c
      do 30 ns=1,nshell
c----------------------------------------------------------------------c
c                 calculate some necessary functions                   c
c                 used over and over again in plm routine              c
c----------------------------------------------------------------------c 
         call iosys ('read real "theta points shell-'//itoc(ns)//
     1               '" from lamdat',nth(ns),z(theta),0,' ')
         call iosys ('read real "theta weights shell-'//itoc(ns)//
     1               '" from lamdat',nth(ns),z(wtthe),0,' ')
         call copy(z(wtthe),z(wtthe+nth(ns)),nth(ns))
         call iosys ('read real "phi points shell-'//itoc(ns)//
     1               '" from lamdat',nph(ns),z(phi),0,' ')
         call iosys ('read real "phi weights shell-'//itoc(ns)//
     1               '" from lamdat',nph(ns),z(wtphi),0,' ')
         call copy(z(wtphi),z(wtphi+nph(ns)),nph(ns))
         call miscfn(z(phi),z(cphi),z(sphi),nph(ns))
c----------------------------------------------------------------------c
c             calculate the legendre functions                         c
c----------------------------------------------------------------------c
         call rzero(z(fio),2*nr(ns)*nth(ns)*nph(ns))
         do 40 mu=0,mumax
            call legend(z(plm),z(theta),z(dfct),z(ddfct),nth(ns),
     1                  lmax,mu,maxfac)
            wds=nth(ns)*(lmax+1)
            wds1=nph(ns)+nph(ns)
            start=plm+nth(ns)*mu
            if (prnt) then
                title='p(l,m) functions for m = '//itoc(mu)
                lentit=length(title)
                len=lmax-mu+1
                call prntfm(title,z(start),nth(ns),len,nth(ns),len,
     1                      iout)
                write(iout,*) 'writing iosys file ', title(1:lentit)
                write(iout,*) 'words written = ',wds
            endif
            call scmphi(z(cphi),z(sphi),z(phifn),nph(ns),mu)
            num=2
            if (mu.eq.0) then
                wds1=nph(ns)
                num=1
            endif
            if (prnt) then
                title='phi angular functions for m = '//itoc(mu)
                lentit=length(title)
                call prntfm(title,z(phifn),nph(ns),num,nph(ns),
     1                      num,iout)
                write(iout,*) 'writing iosys file ', title(1:lentit)
                write(iout,*) 'words written = ',wds1
            endif
            if (tstnrm) then
c               calculate square root of weights
c
                call sqrtwt(z(wtthe+nth(ns)),z(wtphi+nph(ns)),
     1                      nth(ns),nph(ns))
c               scale functions by the new weights in order to
c               make  integral test easy.
c
                call scalfn(z(start),z(wtthe+nth(ns)),nth(ns),lmax,mu)
                call scalfn(z(phifn),z(wtphi+nph(ns)),nph(ns),num,1)
                call legint(z(start),z(phifn),z(scr),z(scr),
     1                      nth(ns),nph(ns),num,lmax,mu)
c               calculate inverse of new weights
c
                call vinv(z(wtthe+nth(ns)),z(wtthe+nth(ns)),nth(ns))
                call vinv(z(wtphi+nph(ns)),z(wtphi+nph(ns)),nph(ns))
c               scale functions by inverse of weights to return them
c               to original state.
c
                call scalfn(z(start),z(wtthe+nth(ns)),nth(ns),lmax,mu)
                call scalfn(z(phifn),z(wtphi+nph(ns)),nph(ns),num,1)
            endif
c           scale functions by weights to make p(l,m) projections
c           easily vectorized.
c           write out scaled weights for later use.
c 
            call scalfn(z(start),z(wtthe),nth(ns),lmax,mu)
            call scalfn(z(phifn),z(wtphi),nph(ns),num,1)
            chvr1=itoc(ns)
            chvr2=itoc(mu) 
            len=length(chvr2)
            len1=length(chvr1)
            title='"weight scaled p(l,'//chvr2(1:len)//')-'
            lentit=length(title) 
            title=title(1:lentit)//chvr1(1:len1)//'"'
            call iosys ('write real '//title//' to lamdat',
     1                   wds,z(plm),0,' ')
            title='"weight scaled phi('//chvr2(1:len)//')-'
            lentit=length(title)
            title=title(1:lentit)//chvr1(1:len1)//'"'
            call iosys ('write real '//title//' to lamdat',wds1,
     1                   z(phifn),0,' ')  
            if (tstint) then
                call iosys ('read real "multicenter yukawa '//
     1                      'potential" from lamdat without '//
     2                      'rewinding',nr(ns)*nth(ns)*nph(ns),z(fio),
     3                       0,' ')
                call decomp(z(fio),z(start),z(phifn),z(flm),z(scr),lmax,
     1                      mu,nr(ns),nth(ns),nph(ns),num,prntlm)
                call vinv(z(wtthe),z(wtthe),nth(ns))
                call vinv(z(wtphi),z(wtphi),nph(ns))
                call scalfn(z(start),z(wtthe),nth(ns),lmax,mu)
                call scalfn(z(phifn),z(wtphi),nph(ns),num,1)
                call plmtof(z(fio),z(start),z(phifn),z(flm),z(scr),
     1                      lmax,mu,nr(ns),nth(ns),nph(ns),num,prntlm)
            endif
   40    continue
         if (tstint) then
             nprd=nr(ns)*nth(ns)*nph(ns)
             call frms(z(fio),z(fio+nprd),rms,nprd)
             write(iout,500) ns,rms
         endif
   30 continue
      call iosys ('write integer maxsiz to rwf',1,canget,0,' ')
      call chainx(0)
      stop
  100 format (//,20x,'***** m6201: calculate spherical harmonics *****')
  300 format(/,5x,'maximum l value',1x,i3,5x,'maximum m value',1x,i3)
  400 format(/,5x,'point buffer',1x,i8)
  500 format(/,5x,'rms error for the ',i2,' region = ',e15.8)
      end






