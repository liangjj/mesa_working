*deck poissson.f
c***begin prologue     m6203
c***date written       930810  (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6203, link 6203, spherical harmonics
c***author             schneider, b. i.(lanl/nsf)
c***source             m6203
c***purpose            solve poisson equation by expansion
c***                   in spherical harmonics.
c***description        solution of the equation  Del**2  V = Rho
c***                   is calculated by expansion in spherical harmonics.
c*** 
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       m6203
      program poisson
      implicit integer (a-z)
      parameter ( dimcen=10 , dimshl=100)
      character*4 chvr1
      character*3 yn, chra, itoc, defp, defnt, defnrm, defdsk, defspl
      character*3 ny
      character*8 cpass
      character*30 str, rtyp, thtyp, phtyp
      character*128 fillam, str1, str2
      character*4096 ops
      character*80 title
      logical logkey, prnty, prntlm, tstnrm, tstint, ylmon
      logical prntyu, noscal, prntnt, prntg, nonsep, ovride
      logical spln, onlysp, nodisk
      real*8 z, znuc, eta, center
      dimension znuc(dimcen), eta(dimcen), center(3,dimcen), nr(dimshl)
      common /io/ inp, iout
      pointer (pz,z(1)), (pz,iz(1))
      call drum
      call iosys ('read character options from rwf',-1,0,0,ops)
c----------------------------------------------------------------------c
c          read in options to print, test normalization of             c
c          spherical harmonics, spline fit or test integral            c
c          calculation                                                 c
c----------------------------------------------------------------------c
      onlysp=logkey(ops,'m6203=only-spherical-harmonics',.false.,' ')
      prnty=logkey(ops,'print=m6203=ylm',.false.,' ')
      prntlm=logkey(ops,'print=m6203=lm-decomposition',.false.,' ')
      prntyu=logkey(ops,'print=m6203=yukawa',.false.,' ')
      noscal=logkey(ops,'m6203=bare-yukawa',.false.,' ')
      prntg=logkey(ops,'print=m6203=grid',.false.,' ')
      nodisk=logkey(ops,'m6203=no-disk-output',.false.,' ')
      call iosys ('read character "linear algebraic filename" from rwf',
     1                 -1,0,0,fillam)
      call iosys ('open lamdat as old',0,0,0,fillam)
      defp='off'
      prntnt=logkey(ops,'print=m6203=partial-integrals',.false.,' ')
      if (prnty.or.prntlm) then
          defp='on'
      endif    
      tstnrm=logkey(ops,'m6203=test-normalization',.false.,' ')
      defnrm='off'
      if (tstnrm) then
          defnrm='on'
      endif      
      tstint=logkey(ops,'m6203=test-yukawa-integral',.false.,' ')
      defnt='off'
      if (tstint) then
          defnt='on'
      endif    
      ylmon=logkey(ops,'m6203=ylm=on-disk',.false.,' ')
      defdsk='off'
      if (ylmon) then
          defdsk='on'
      endif    
      spln=logkey(ops,'m6203=no-spline',.false.,' ')
      defspl='on'
      if (spln) then
          defspl='off'
      else    
          order=intkey(ops,'m6203=order-of-spline-fit',4,' ')
          call iosys ('write integer "order of spline fit" '//
     1                'to lamdat',1,order,0,' ')
      endif
c----------------------------------------------------------------------c
c            read in basic information on centers, points,             c
c            angular integration type, l and m values, energies        c
c----------------------------------------------------------------------c
      write (iout,1)
      write (iout,2) defp, defnrm, defnt, defdsk, defspl, order
      call iosys ('read integer "max number of shells" from '//
     1            'lamdat',1,maxshl,0,' ')
      if (maxshl.gt.dimshl) then
          call lnkerr('quit. number of shells too big')
      endif    
      call iosys ('read integer "number of atomic centers" from lamdat',
     1            1,ncen,0,' ')
      call iosys('read character "scattering center" from lamdat',-1,0,
     1            0,yn)
      call iosys ('read real "yukawa exponents" from lamdat',ncen,
     1             eta,0,' ')
      call iosys ('read real "nuclear charges" from lamdat',ncen,
     1             znuc,0,' ')
      do 10 i=1,ncen
         call iosys ('read real "atomic center positions" from '//
     1               'lamdat without rewinding',3,center(1,i),0,' ')
   10 continue
      ncplus=ncen 
      if (yn.eq.'yes') then
          ncplus=ncplus+1 
          call iosys ('read real "scattering center position" '//
     1                'from lamdat',3,center(1,ncplus),0,' ')
      endif
      call iosys ('read integer "max theta points" from lamdat',1,
     1             nthmax,0,' ')
      call iosys ('read integer "max phi points" from lamdat',1,
     1             nphmax,0,' ')
      call iosys('read character "non separable quadrature" '//
     1           'from lamdat',-1,0,0,ny)
      nonsep=.false. 
      if (ny.eq.'yes') then
          nonsep=.true.
          call iosys ('read integer "max number lebedev points" '//
     1                'from lamdat',1,maxleb,0,' ')
          if (maxleb.ne.nthmax) then
              call lnkerr('error in angular points:nthmax')
          endif
          if (maxleb.ne.nphmax) then
              call lnkerr('error in angular points:nphmax')
          endif                                                         
      endif
      totpt=0
      nrmax=0
      totr=0    
      if (.not.onlysp) then
           call iosys ('read integer "max grid points" from lamdat',1,
     1                  totpt,0,' ')
           call iosys('read integer "max radial points in shell" '//
     1                'from lamdat',1,nrmax,0,' ')
      call iosys('read integer "max radial points" from lamdat',1,
     1            totr,0,' ')
      endif
      write (iout,3) ncen, ncplus, nrmax, totr, nthmax, nphmax, totpt
      if (nonsep) then
          write(iout,4) maxleb
      endif
      call iosys ('read integer "biggest l" from lamdat',1,
     1             lmax,0,' ')
      call iosys ('read integer "biggest m" from lamdat',1,
     1             mumax,0,' ')
      ovride=logkey(ops,'m6203=override-angular-momentum',.false.,' ')
      if (ovride) then
          lmax=intkey(ops,'m6203=lmax',lmax,' ')
          mumax=intkey(ops,'m6203=mmax',mumax,' ')
      endif    
      lplus=lmax+1
      mplus=mumax+1 
      write (iout,6) mumax, lmax
c----------------------------------------------------------------------c
c                     memory allocation                                c
c----------------------------------------------------------------------c
c
c     set up some integers giving maximum values of dimensions and
c
      lsize=0
      do 20 lstor=0,mumax
         lsize=lsize+lmax-lstor+1
   20 continue          
      tlsize=lsize+lsize
      maxfac=lplus+mplus+10
      maxfcd=maxfac+1
      maxang=max(nthmax,nphmax)
      lsq=lplus*lplus
      wtdim=totr*(totr-1)
      nsplw=(4*totr+1)*order*2
      nbreak=totr-order+1
      maxscr=max(nthmax,2*totr*nthmax,2*totr*lplus,4*totr)
      if (tstnrm) then
          maxscr=max(maxscr,tlsize*tlsize)
      endif               
      ioff=1
      do 30 i=1,2
c
c         memory allocation for factorials, grid points, weights, sine and
c         cosine, phi functions and legendre polynomials.      
          dfct=ioff
          ddfct=dfct+maxfcd
          if (nonsep) then
             angle=ddfct+maxfcd
             theta=angle
             phi=theta+maxleb 
             wtang=phi+maxleb
             sphi=wtang+maxleb+maxleb
             cphi=sphi+maxleb
c            the doubling of the dimension for plm allows the storage 
c            for the second half of the y(l,m) array to overlap p(l,m).  
             ylm=cphi+maxleb
             plm=ylm+maxleb*lsize 
             phifn=plm+maxleb*lsize
             rad=phifn+2*mplus*maxleb
          else
             theta=ddfct+maxfcd
             phi=theta+nthmax
             wtthe=phi+nphmax
             wtang=wtthe
             wtphi=wtthe+nthmax+nthmax
             sphi=wtphi+nphmax+nphmax
             cphi=sphi+nphmax
             plm=cphi+nphmax
             phifn=plm+nthmax*lsize
             rad=phifn+2*mplus*nphmax                      
          endif
          words=rad
          wtrad=rad+totr+totr
c         
c         memory allocation for yukawa potential, bessel functions, legendre
c         decomposed functions and solutions of inhomogeneous wave equation.
c          
          if (.not.onlysp) then
              yuk=wtrad+wtdim+wtdim
              rho=yuk+totpt
              j=rho+totpt
              y=j+totr*lplus
              flm=y+totr*lplus
              psilm=flm+totr*tlsize
              scr=psilm+totr*tlsize
              words=wpadti(scr+maxscr)
              if (.not.spln) then
                  sc=scr+maxscr
                  xtmp=sc+nsplw
                  break=xtmp+totr
                  c=break+nbreak+1
                  ind=wpadti(c+order*nbreak)
                  words=ind+totr
              endif
          endif
          if (i.eq.1) then
              call iosys ('read integer maxsiz from rwf',1,
     1                     canget,0,' ')
              if (words.gt.canget) then
                  call lnkerr('not enough memory. will quit')
              endif
              call iosys ('write integer maxsiz to rwf',1,
     1                     words,0,' ')
              call getscm(words,z,ngot,'m6203',0)
              write (iout,*) 'get ',words,' words of core'      
          endif
   30 continue
c----------------------------------------------------------------------c 
c             calculate factorials and double factorials needed        c
c             for computation of legendre functions.                   c   
c----------------------------------------------------------------------c 
      call fact(z(dfct),z(ddfct),maxfac)
c----------------------------------------------------------------------c
c          loop over the atoms to find the atomic contribution         c
c          to the property.                                            c
c----------------------------------------------------------------------c
      do 40 atom=1,ncplus
c----------------------------------------------------------------------c
c          read in the necessary basic information and then either     c
c          compute or read in the spherical harmonics.                 c
c----------------------------------------------------------------------c
         if (atom.lt.ncplus) then
             cpass='atom'
             chra=itoc(atom)
             len=length(chra)
             str='atom-'//chra(1:len)
         else
             if (yn.ne.'yes') then
                 cpass='atom'
                 chra=itoc(atom)
                 len=length(chra)
                 str='atom-'//chra(1:len)
             else
                 cpass='scattering'
                 str='scattering center'
             endif
         endif           
         write (iout,*) 'beginning calculation on '//str
c
c        we branch here depending on the type of angular quadrature
c
         if (nonsep) then
             call iosys ('read integer "number lebedev points '//
     1                    str//'" from lamdat',1,nang,0,' ')     
             call iosys ('read real "lebedev angular points '//
     1                    str//'" from lamdat',2*nang,
     1                    z(angle),0,' ')
             call iosys ('read real "lebedev angular weights '//
     1                    str//'" from lamdat',nang,
     2                    z(wtang),0,' ')
             nthet=nang
             nphi=nang
         else                
             call iosys('read integer "theta quadrature order '//
     1                   str//'" from lamdat',1,nthet,0,' ')  
             call iosys('read integer "phi quadrature order '//
     1                   str//'" from lamdat',1,nphi,0,' ') 
             call iosys('read character "theta quadrature type '//
     1                   str//'" from lamdat',-1,0,0,thtyp)
             call iosys('read character "phi quadrature type '//
     1                   str//'" from lamdat',-1,0,0,phtyp)
             call iosys ('read real "theta points '//str//
     1                   '" from lamdat',nthet,z(theta),0,' ')
             call iosys ('read real "theta weights '//str//
     1                   '" from lamdat',nthet,z(wtthe),0,' ')
             call iosys ('read real "phi points '//str//
     1                   '" from lamdat',nphi,z(phi),0,' ')
             call iosys ('read real "phi weights '//str//
     1                   '" from lamdat',nphi,z(wtphi),0,' ')
         endif     
         call iosys('read integer "max l value '//
     1               str//'" from lamdat',1,lval,0,' ')
         call iosys ('read integer "max m value '//
     1               str//'" from lamdat',1,mval,0,' ')
         if (ovride) then
             lval=lmax
             mval=mumax 
         endif         
         if (.not.onlysp) then
c        
c        read radial shell information and associated points
c
             call iosys ('read integer "number of shells '//
     1                    str//'" from lamdat',1,nshell,0,' ')
             call iosys ('read integer "number radial points per '//
     1                   'shell '//str//'" from lamdat',nshell,nr,0,' ')
             call iosys('read character "radial quadrature type '//
     1                   str//'" from lamdat',-1,0,0,rtyp)
             locr=rad
             locwt=wtrad
             nrtot=0
             nints=0
             nwgts=0
             locru=xtmp
             do 50 ns=1,nshell
                call iosys ('read real "radial points '//
     1                      str//'" from lamdat without rewinding',
     2                      nr(ns),z(locr),0,' ')
                nwts=nr(ns)
                if (rtyp.eq.'newton-cotes') then
                    nwts=nr(ns)*(nr(ns)-1)
                    nwgts=nwgts+nwts
                    nints=nints+nr(ns)-1
                    if (.not.spln) then
                         call runiqe(z(locr),z(locru),rtyp,nr(ns))
                         locru=locru+nr(ns)-1
                    endif     
                else
                    if (.not.spln) then
                         call runiqe(z(locr),z(locru),rtyp,nr(ns))
                         locru=locru+nr(ns)
                    endif     
                    nints=nints+nr(ns)
                    nwgts=nwgts+nwts
                endif
                call iosys ('read real "unscaled radial weights '//
     1                       str//'" from lamdat without '//
     2                       'rewinding',nwts,z(locwt),0,' ')
                locr=locr+nr(ns)
                locwt=locwt+nwts
                nrtot=nrtot+nr(ns)
   50        continue
         endif
         str1='"y(l,m) '//str//'"'
         str2='"phim '//str//'"'       
         if (.not.ylmon) then
             nwds=(lval+1)*(mval+1)*nthet
             if (.not.nodisk) then
                  call iosys ('create real '//str1//' on lamdat',
     1                         nwds,0,0,' ')
                  nwds=(2*mval+1)*nphi
                  call iosys ('create real '//str2//' on lamdat',
     1                         nwds,0,0,' ')
             endif
         else
             call iosys ('rewind '//str1//' on lamdat '//
     1                   'read-and-write',0,0,0,' ')
             call iosys ('rewind '//str2//' on lamdat '//
     1                   'read-and-write',0,0,0,' ')
         endif
c----------------------------------------------------------------------c
c             begin section of code to calculate the spherical         c
c             functions. by setting appropriate flags the              c
c             normalization integrals can be examined.                 c
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c             calculate cosine and sine of phi.                        c
c             they are used over and over again in plm routine.        c
c----------------------------------------------------------------------c
         call miscfn(z(phi),z(cphi),z(sphi),nphi)
c----------------------------------------------------------------------c
c             compute the legendre and phifn functions.                c
c----------------------------------------------------------------------c 
         lv1=lval+1 
         locplm=plm
         locphm=phifn
         locylm=ylm
         totplm=0
         totylm=0
         totphm=0
         do 60 m=0,mval
            nm=2
            if (m.eq.0) then
                nm=1
            endif
            rowlen=lv1-m
            wds=nthet*rowlen
            wds1=nm*nphi
            if (.not.ylmon) then
c----------------------------------------------------------------------c
c              compute the spherical harmonics                         c
c----------------------------------------------------------------------c     
                call legend(z(locplm),z(theta),z(dfct),z(ddfct),nthet,
     1                      lval,m,maxfac)
                call scmphi(z(phi),z(locphm),nphi,m)
                if (.not.nodisk) then
                     chvr1=itoc(m) 
                     len=length(chvr1)
                     title='"p(l,'//chvr1(1:len)//')"'
                     lentit=length(title)
                     write(iout,*) 'writing iosys file ', 
     1                              title(1:lentit)
                     write(iout,*) 'words written = ',wds
                     call iosys ('write real '//str1//' to lamdat '//
     1                           'without rewinding',wds,
     2                            z(locplm),0,' ')
                     title='"phi('//chvr1(1:len)//')"'
                     lentit=length(title)
                     write(iout,*) 'writing iosys file ', 
     1                              title(1:lentit)
                     write(iout,*) 'words written = ',wds1
                     call iosys ('write real '//str2//' to lamdat '//
     1                      'without rewinding',wds1,z(locphm),0,' ')
                endif
                title='p(l,m) functions for m = '//itoc(m)
                lentit=length(title)
                if (prnty) then
                    call prntfm(title,z(locplm),nthet,rowlen,nthet,
     1                          rowlen,iout)
                endif
                title='phi angular functions for m = '//itoc(m)
                lentit=length(title)
                if (prnty) then
                    call prntfm(title,z(locphm),nphi,nm,nphi,nm,iout)
                endif
c----------------------------------------------------------------------c
c               form spherical harmonics if required                   c
c----------------------------------------------------------------------c
                if (nonsep) then
                    call mkylm(z(locplm),z(locphm),z(locylm),nang,
     1                         lval,m)
                endif                   
            else
c----------------------------------------------------------------------c
c              the spherical harmonics have been computed previously   c
c              read them in from disk.                                 c
c----------------------------------------------------------------------c      
                call iosys ('read real '//str1//' on lamdat without '//
     1                      'rewinding',wds,z(locplm),0,' ') 
                call iosys ('read real '//str2//' on lamdat without '//
     1                      'rewinding',wds1,z(locphm),0,' ')
c----------------------------------------------------------------------c
c               form spherical harmonics if required                   c
c----------------------------------------------------------------------c
                if (nonsep) then
                    call mkylm(z(locplm),z(locphm),z(locylm),nang,
     1                         lval,m)
                endif                 
            endif
            totplm=totplm+rowlen
            totylm=totylm+nm*rowlen
            totphm=totphm+nm
            locplm=locplm+wds
            locphm=locphm+wds1
            locylm=locylm+rowlen*nang*nm
   60    continue
c----------------------------------------------------------------------c
c             spherical harmonics are computed.  calculate the         c
c             normalization integral if required.                      c
c----------------------------------------------------------------------c 
c----------------------------------------------------------------------c
c              need copies of weights                                  c
c----------------------------------------------------------------------c
         if (nonsep) then
             call copy(z(wtang),z(wtang+nang),nang)
         else    
             call copy(z(wtthe),z(wtthe+nthet),nthet)
             call copy(z(wtphi),z(wtphi+nphi),nphi)
         endif    
         if (tstnrm) then
             call nrmylm(z(plm),z(phifn),z(ylm),z(wtthe),
     1                   z(wtphi),z(wtang),z(scr),nthet,nphi,
     2                   nang,lval,mval,totplm,totylm,totphm,nonsep)
         endif
c----------------------------------------------------------------------c
c              now on to the main business.  in what follows we        c
c              actually calculate radial functions which are           c
c              eventually r times the functions we want.  this         c
c              simplifies the calculation of the green's function      c
c              and leads to better numerical behavior of the           c
c              functions which are spline fit.                         c
c                                                                      c
c      G(r;l,m)'' -[ l*(l+1)/(r*r) ] G(r;l,m) = -4*pi*rho(r;l,m)       c
c                   G(r;l,m) = V(r;l,m)*r                              c
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
         if (.not.onlysp) then
c----------------------------------------------------------------------c
c             let us assume that the inhomogeneity is a                c
c             multicenter yukawa potential.                            c  
c----------------------------------------------------------------------c
             if (nonsep) then
                 call yukns(z(yuk),eta,center,z(rad),z(theta),
     1                      z(cphi),z(sphi),z(scr),ncen,nshell,nrtot,
     2                      nr,nthet,nphi,nang,nonsep)
             else
                 call yuks(z(yuk),eta,center,z(rad),z(theta),
     1                     z(cphi),z(sphi),z(scr),ncen,nshell,nrtot,
     2                     nr,nthet,nphi,nang,nonsep)
             endif
c----------------------------------------------------------------------c
c           calculate the spline information only dependent on the     c
c           grid used not on the fitted function.                      c
c----------------------------------------------------------------------c
             if (.not.spln) then
                 call prespl(nints,z(xtmp),order,z(sc))
                 call splmat(nints,z(xtmp),order,z(sc))
                 nbk=nints-order+1
                 nspls=(4*nints+1)*order
                 str1='"number of break points '//str//'"'
                 call iosys ('write integer '//str1//' to lamdat',1,
     1                        nbk,0,' ')
             endif
             if (prntg) then    
                 title='radial grid'
                 call prntfm(title,z(rad),nrtot,1,nrtot,1,iout)
             endif    
c----------------------------------------------------------------------c
c            calculate inverse of weights                              c
c----------------------------------------------------------------------c
             if (nonsep) then
                 call vinv(z(wtang+nang),z(wtang+nang),nang)
             else
                 call vinv(z(wtthe+nthet),z(wtthe+nthet),nthet)
                 call vinv(z(wtphi+nphi),z(wtphi+nphi),nphi)
             endif
c----------------------------------------------------------------------c
c             decompose inhomogeneity into spherical harmonics.        c 
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c             calculate regular and irregular solutions to             c
c             homogeneous equation.                                    c
c----------------------------------------------------------------------c
             call homog(z(rad),z(j),z(y),nrtot,lval)
             ncpy=nrtot*nthet*nphi
             if ( nonsep ) then
                  ncpy=nrtot*nang
             endif
             call copy(z(yuk),z(rho),ncpy)
c----------------------------------------------------------------------c 
c             decompose the function into its angular components       c
c----------------------------------------------------------------------c 
             call lmdcmp(z(rho),z(plm),z(phifn),z(ylm),z(flm),
     1                   z(wtthe),z(wtphi),z(wtang),z(scr),lval,mval,
     2                   nrtot,nthet,nphi,nang,nonsep,prntlm)
             call plmtof(z(rho),z(plm),z(phifn),z(ylm),z(flm),z(scr),
     1                   z(rad),lval,mval,nrtot,nthet,nphi,nang,
     2                   nonsep,prntlm)
c----------------------------------------------------------------------c
c             from the radial components of the inhomogeneity          c
c             calculate the radial component of the solution           c
c             using the green's function to solve the integral         c
c             equation.                                                c
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c               do the integral                                        c
c----------------------------------------------------------------------c
             call drvint(z(flm),z(psilm),z(rad),z(j),z(y),z(wtrad),
     1                   z(scr),nrtot,lval,mval,nshell,nr,nints,nwgts,
     2                   rtyp,prntnt)
             call plmtof(z(rho),z(plm),z(phifn),z(ylm),z(psilm),
     1                   z(scr),z(rad),lval,mval,nints,nthet,nphi,
     2                   nang,nonsep,prntlm)
c----------------------------------------------------------------------c
c               get the spline coefficients                            c
c----------------------------------------------------------------------c
             if (.not.spln) then
                 str1='"spline coefficients '//str//'"'
                 wrds=2*(mval+1)*(lval+1)*order*nbk 
                 call iosys ('create real '//str1//' on lamdat',wrds,
     1                        0,0,' ')
                 call mkspln(z(psilm),z(break),z(c),z(sc),mval,lval,
     1                       nints,order,nbk,nspls,str1)
                 call caspln(z(psilm),z(xtmp),z(break),z(c),a(ind),
     1                       mval,lval,nints,order,nbk,nspls,str1,
     2                       prntlm)
             endif  
         endif
   40 continue
      call iosys ('write integer maxsiz to rwf',1,canget,0,' ')
      call chainx(0)
      stop
    1 format (//,20x,'***** m6203: solve wave equation *****')
    2 format(/,5x,'defaults:',//,5x,'print = ',a3,3x,'test norm ylm = '
     1         ,a3,/,5x,'test yukawa integral = ',a3,3x,
     2                  'ylm on disk = ',a3,/,5x,'spline fit = ',a3,3x,
     3                  'order of spline fit = ',i1)
    3 format(/,5x,'atomic data:',
     1      /,5x, 'no. atomic centers = ',i3,3x,'total no. centers = ',
     2      i3,/,5x,'max. radial points in shell = ',i4,3x,
     3              'max. radial points = ',i4,
     4        /5x,  'max. theta points = ',i4,3x,'max. phi points = ',
     5       i3,/,5x,  'max. grid points = ',i6)         
    4 format (/,5x,'using lebedev quadrature. max. lebedev points = ',
     1              i3)
    5 format (/,5x,'number of energies = ',i3,1x,'and their values',
     1             (/,5x,5(e15.8,1x)))
    6 format (/,5x,'maximum m value',1x,i3,5x,'maximum l value',1x,i3)
      end
