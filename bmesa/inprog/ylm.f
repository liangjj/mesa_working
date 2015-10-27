*deck ylm.f
c***begin prologue     m6201
c***date written       921219   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6201, link 6201, spherical harmonics
c***author             schneider, b. i.(lanl/nsf)
c***source             m6201
c***purpose            solve inhomogeneous wave equation by expansion
c***                   in spherical harmonics.
c***description        solution of the equation ( Del**2 + k**2 ) V = Rho
c***                   is calculated by expansion in spherical harmonics.
c***                   we allow k to be zero to handle the poisson equation.
c*** 
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       m6201
      program ylm 
      implicit integer (a-z)
      parameter ( dimcen=10 , dimshl=100, dimtmp=200, dimene=50, acc=30)
      character*4 chvr1
      character*4 itoc
      character*8 cpass
      character*20 radqud, chrkey
      character*128 fillam
      character*1600 card
      character*4096 ops
      character*80 title
      logical logkey, prnt, prntlm, posinp, tstnrm, tstint
      real*8 z, dummy, center, alpha, rr, rshel, rbox, energy
      real*8 tmp
      dimension z(1), dummy(2)
      dimension nrad(dimshl), rshel(dimshl), alpha(dimcen)
      dimension center(3,dimcen), mval(dimtmp), ibuf(dimtmp), nl(dimtmp)
      dimension energy(dimene), k, tmp
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
      write (iout,1)
      call iosys ('read character "linear algebraic filename" from rwf',
     1                 -1,0,0,fillam)
      call iosys ('open lamdat as old',0,0,0,fillam)
c----------------------------------------------------------------------c
c          if m is positive it denotes cosine like phi behavior.       c
c          negative m is sinelike behavior.                            c                         
c----------------------------------------------------------------------c      
      if ( posinp ('$ylm',cpass) ) then
           call cardin(card)
           nm=intkey(card,'number-of-m-values',1,' ')   
           call iosys ('write integer "number m values" to lamdat',
     1                  1,nm,0,' ')
           call intarr(card,'m-values',mval,nm,' ')
           call iosys ('write integer "m values" to lamdat',
     1                  nm,mval,0,' ')
           call intarr(card,'number-l-values-for-each-m',nl,nm,' ')
           call iosys ('write integer "number of l values for each m"'//
     1                 ' to lamdat',nm,nl,' ')
           lmax=0
           mumax=0
           totl=0
           do 10 mmax=1,nm
              mumax=max(mumax,abs(mval(mmax)))
              totl=totl+nl(mmax)
              call intarr(card,'l-values-for-m-'//itoc(mval(mmax)),
     1                    ibuf,nl(mmax),' ')
              call order(ibuf,ibuf,nl(mmax),'integer')
              call iosys ('write integer "l values for m-'//
     1                    itoc(mval(mmax))//'" to lamdat',nl(mmax),
     2                    ibuf,' ')
              do 20 l=1,nl(mmax)
                 lmax=max(lmax,ibuf(l))
   20         continue
   10      continue
      endif
c
c            we may drop lm pairs in a given region by using zero
c            radial quadrature points for that pair in that region
c
      if ( posinp ('$radial',cpass) ) then
           call cardin(card)
           nshell=intkey(card,'number-of-radial-shells',1,' ')
           radqud=chrkey(card,'type-radial-quadrature',
     1                   'newton-cotes',' ')
           call iosys ('write integer "number of radial shells" '//
     1                 'to lamdat',1,nshell,0,' ')
           write (iout,3) nshell
           call intarr(card,'number-of-radial-points-per-shell',nrad,
     1                 nshell,' ')
           call iosys ('write integer "number radial points per '//
     1                 'shell" to lamdat',nshell,nrad,' ') 
           write (iout,4) ( nrad(ns),ns=1,nshell)
           call fparr(card,'shell-radii',rshel,nshell,' ')
           write (iout,5) (rshel(ns),ns=1,nshell)
           maxr=0
           totr=0
           do 30 ns=1,nshell
              maxr=max(maxr,nrad(ns))
              totr=totr+nrad(ns)
   30      continue
           rbox=rshel(nshell)
           lenwt=totr
           if (radqud.eq.'newton-cotes') then
               lenwt=0
               do 40 ns=1,nshell
                  lenwt=lenwt+nrad(ns)*(nrad(ns)-1)
   40          continue
           endif                 
      endif
      if ( posinp ('$energy',cpass) ) then
           call cardin(card)
           nen=intkey(card,'number-of-energies',1,' ')
           call fparr(card,'energies',energy,nen,' ')
      endif                                 
      write (iout,6) lmax, mumax
      call iosys ('write integer "maximum l in ylm" to lamdat',1,
     1             lmax,0,' ')
      call iosys ('write integer "maximum m in ylm" to lamdat',1,
     1             mumax,0,' ')
      write(iout,8) nen, (energy(ii),ii=1,nen)
c       
c             based on the values of lmax and mumax we can compute
c             the maximum number of quadrature points in theta and
c             phi to integrate a polynomial integrand exactly of a
c             given order.
c
      maxfac=lmax+mumax+10
      maxth=lmax+1
      maxph=2*mumax+1
      maxang=max(maxth,maxph)
      maxall=max(maxang,maxr)
      totpt=maxph*maxth*totr
      if (posinp('$yukawa',cpass)) then
          call cardin(card)
          ncen=intkey(card,'number-yukawa-centers',1,' ')
          do 50 cen=1,ncen
             call intarr(card,'position-center-'//itoc(cen),
     1                   center(1,cen),3,' ')
             alpha(cen)=fpkey(card,'yukawa-exponent-center-'
     1                       //itoc(cen),1.d0,' ')
   50     continue           
      endif
c----------------------------------------------------------------------c
c                     memory allocation                                c
c                          RE DO                                       c
c----------------------------------------------------------------------c
      lplus=lmax+1
      lsq=lplus*lplus
      maxscr=max(maxall,totl,lplus*maxth,2*maxr*maxth,
     1           2*maxr*(ltop1+1))
      if (tstnrm) then
          maxscr=max(maxscr,lsq,4)
      endif               
      tmp=lmax*acc
      ltop=lmax+sqrt(tmp)
      ltop=max(strtl,lmax)
      ltop1=ltop+1
      ioff=1
      do 60 i=1,2
          dfct=ioff
          ddfct=dfct+maxfac
          theta=ddfct+maxfac
          wtthe=theta+maxth
          phi=theta+maxth+maxth
          wtphi=phi+maxph
          rad=wtphi+maxph+maxph
          wtrad=rad+totr+totr
          sphi=wtrad+lenwt+lenwt
          cphi=sphi+maxph
          phifn=cphi+maxph
          plm=phifn+nm*maxph
          yuk=plm+maxth*totl
          j=yuk+totpt
          y=j+totr*ltop1
          wron=y+totr*ltop1
          flm=wron+ltop1
          psilm=flm+totr*totl
          lndex=wpadti(psilm+totr*totl)
          jp=iadtwp(lndex+totl)
          yp=jp+maxr*ltop1
          scr=yp+maxr*ltop1
          plmscr=jp
          words=wpadti(scr+maxscr)
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
   60 continue
c----------------------------------------------------------------------c          
c             calculate the gauss-legendre quadrature points for the   c
c             plm functions and the simpson quadrature points for      c 
c             the phifn integration.                                   c
c----------------------------------------------------------------------c
      call gaussq('legendre',maxth,0.d0,0.d0,0,dummy,z(scr),
     1             z(theta),z(wtthe))
      call gaussq('simpson',maxph,0.d0,0.d0,0,dummy,z(scr),
     1             z(phi),z(wtphi))
c----------------------------------------------------------------------c     
c     the nature of the weights and points is very different for       c
c     newton-cotes and gauss quadrature. in the newton-cotes           c
c     formulas the end ranges of the points overlap and the weights    c
c     for each region are a matrix of dimension n*(n-1).               c 
c----------------------------------------------------------------------c 
      rr=0.d0
      locr=rad
      locwtr=wtrad
      do 70 ns=1,nshell
         if (radqud.eq.'legendre') then
             call gaussq('legendre',nrad(ns),0.d0,0.d0,0,z(scr),
     1                    z(locr),z(locwtr))
             call convr(rr,rshel(ns),z(locr),z(locwtr),nrad(ns))
             locwtr=locwtr+nrad(ns)
         elseif (radqud.eq.'newton-cotes') then
             call necote(rr,rshel(ns),z(locr),z(locwtr),nrad(ns))
             locwtr=locwtr+nrad(ns)*(nrad(ns)-1)
         else
             call lnkerr('error in radial quadrature type')
         endif
         locr=locr+nrad(ns)
         rr=rshel(ns)
   70 continue     
c----------------------------------------------------------------------c 
c             calculate factorials and double factorials needed        c
c             for computation of legendre functions.                   c   
c----------------------------------------------------------------------c 
      call fact(z(dfct),z(ddfct),maxfac)
c----------------------------------------------------------------------c
c             begin section of code to calculate the spherical         c
c             functions. by setting appropriate flags the              c
c             normalization integrals can be examined.                 c
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c             calculate cosine and sine of phi.                        c
c             they are used over and over again in plm routine.        c
c----------------------------------------------------------------------c
      call miscfn(z(phi),z(cphi),z(sphi),maxph)
c----------------------------------------------------------------------c
c             calculate the yukawa potential                           c
c----------------------------------------------------------------------c
      locr=rad
      locyuk=yuk
      do 80 ns=1,nshell
         call yukawa(z(locyuk),alpha,center,z(locr),z(theta),z(cphi),
     1               z(sphi),ncen,nrad(ns),nthet,nphi)
         locr=locr+nrad(ns)
         locyuk=locyuk+nrad(ns)*nthet*nphi
  80  continue         
c----------------------------------------------------------------------c
c              need copies of weights                                  c
c----------------------------------------------------------------------c
      call copy(z(wtthe),z(wtthe+maxth),maxth)
      call copy(z(wtphi),z(wtphi+maxph),maxph)
c----------------------------------------------------------------------c
c             compute the legendre and phifn functions.                c
c----------------------------------------------------------------------c 
      locplm=plm
      locphm=phifn
      loclnd=lndex
      do 90 m=1,nm
         mu=mval(m)
         absmu=abs(mu)
         call legend(z(plmscr),z(theta),z(dfct),z(ddfct),maxth,
     1               lmax,absmu,maxfac)
         call strleg(z(plmscr),z(locplm),a(loclnd),maxth,lmax,
     1               absmu,nl(m))
         loclnd=loclnd+nl(m)
         title='p(l,m) functions for m = '//itoc(mu)
         lentit=length(title)
         wds=maxth*nl(m)
         if (prnt) then
             call prntfm(title,z(locplm),maxth,nl(m),maxth,nl(m),iout)
         endif
         call scmphi(z(cphi),z(sphi),z(locphm),maxph,mu)
         title='phi angular function for m = '//itoc(mu)
         lentit=length(title)
         if (prnt) then
             call prntfm(title,z(phifn),maxph,1,maxph,
     1                   1,iout)
         endif
         if (tstnrm) then
c
c               calculate square root of weights
c
             call sqrtwt(z(wtthe+maxth),z(wtphi+maxph),maxth,maxph)
c
c               scale functions by the new weights in order to
c               make  normalization test integral easy.
c
             call scalfn(z(locplm),z(wtthe+maxth),maxth,nl(m),1)
             call scalfn(z(locphm),z(wtphi+maxph),maxph,1,1)
             call legint(z(locplm),z(locphm),z(scr),maxth,maxph,
      1                  nl(m),1)
c
c               calculate inverse of new weights
c
             call vinv(z(wtthe+maxth),z(wtthe+maxth),maxth)
             call vinv(z(wtphi+maxph),z(wtphi+maxph),maxph)
c
c               scale functions by inverse of weights to return them
c               to original state.
c
             call scalfn(z(locplm),z(wtthe+maxth),maxth,nl(m),1)
             call scalfn(z(locphm),z(wtphi+maxph),maxph,1,1)
         endif
c
c           scale functions by weights to make p(l,m) projections
c           easily vectorized.
c           write out scaled weights for later use.
c 
         call scalfn(z(locplm),z(wtthe),maxth,nl(m),1)
         call scalfn(z(locphm),z(wtphi),maxph,1,1)
         chvr1=itoc(mu) 
         len=length(chvr1)
         title='"weight scaled p(l,'//chvr1(1:len)//')"'
         lentit=length(title)
         call iosys ('write real '//title//' to lamdat',
     1                wds,z(locplm),0,' ')
         write(iout,*) 'writing iosys file ', title(1:lentit)
         write(iout,*) 'words written = ',wds
         title='"weight scaled phi('//chvr1(1:len)//')"'
         lentit=length(title)
         call iosys ('write real '//title//' to lamdat',maxph,
     1                z(locphm),0,' ')  
         write(iout,*) 'writing iosys file ', title(1:lentit)
         write(iout,*) 'words written = ',maxph
         locplm=locplm+wds
         locphm=locphm+maxph
   90 continue
c
c----------------------------------------------------------------------c   
c             begin energy loop                                        c
c----------------------------------------------------------------------c
      do 100 ien=1,nen
         k=sqrt(energy(ien))
c----------------------------------------------------------------------c
c             begin loop over shells, decomposition of                 c
c             inhomogeneity into spherical harmonics shell by shell.   c 
c----------------------------------------------------------------------c
         locyuk=yuk
         locflm=flm
         locr=rad
         locj=j
         locy=y           
         do 110 ns=1,nshell
c----------------------------------------------------------------------c
c             calculate bessel's of kr                                 c
c----------------------------------------------------------------------c 
            call vscale(z(locr+totr),z(locr),k,nrad(ns))
            call rcbes(z(locr+totr),z(locj),z(jp),z(locy),z(yp),
     1                 z(wron),scr,nrad(ns),lmax,ltop,
     2                 'derivatives',.false.)
c----------------------------------------------------------------------c
c             the next bit of monkey business converts derivatives     c
c             wrt k*r to derivatives wrt r, recomputes the wronskian   c
c             which is needed for the green's function and then        c       
c             scales the both functions so that the proper factor      c
c             will be incorporated into the radial integral.           c
c----------------------------------------------------------------------c  
            call sscal(ltop1,k,z(wron),1)
            call vinv(z(wron),z(wron),ltop1)
            call vsqrt(z(wron),z(wron),ltop1)
            call mvmul(z(wron),z(locj),z(locj),nrad(ns),ltop1)
            call mvmul(z(wron),z(locy),z(locy),nrad(ns),ltop1)
c----------------------------------------------------------------------c 
c             decompose the function into its angular components       c
c----------------------------------------------------------------------c      
            call lmdcmp(z(locyuk),z(plm),z(phifn),z(locr),z(locflm),
     1                  z(scr),mval,nl,nrad(ns),maxth,maxph,nm,prntlm)
            locr=locr+nrad(ns)
            locyuk=locyuk+nrad(ns)*maxth*maxph
            locflm=locflm+nrad(ns)*totl
            locj=locj+nrad(ns)*ltop1
            locy=locy+nrad(ns)*ltop1    
  110    continue
c----------------------------------------------------------------------c
c             from the radial components of the inhomogeneity          c
c             calculate the radial component of the solution           c
c             using the green's function to solve the integral         c
c             equation.                                                c
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c               do the forward integral                                c
c----------------------------------------------------------------------c
         call drvfrd(z(psilm),z(flm),z(j),z(y),z(wtrad),z(scr),nrad,
     1               mval,a(loclnd),nl,nm,ltop,totl,nshell,radqud)
  100 continue
      call iosys ('write integer maxsiz to rwf',1,canget,0,' ')
      call chainx(0)
      stop
    1 format (//,20x,'***** m6201: solve wave equation *****')
    2 format (/,5x,'the number of (l,m) pairs = ',i3,' their values are'
     1         ,/,(5x,5(:,'(',i3,',',i3,')')))
    3 format (/,5x,'number of radial shells = ',i4)
    4 format (/,5x,'number radial points per shell',(/,10(i4,1x)))
    5 format (/,5x,'shell radii',(/,5(1x,e15.8)))
    6 format (/,5x,'maximum l value',1x,i3,5x,'maximum m value',1x,i3)
    7 format (/,5x,'rms error for the ',i2,' region = ',e15.8)
    8 format (/,5x,'number of energies = ',i3,1x,'and their values',
     1              (/,5x,5(e15.8,1x)))
      end
