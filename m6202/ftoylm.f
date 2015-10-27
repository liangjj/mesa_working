*deck ftoylm.f
c***begin prologue     m6202
c***date written       930810  (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6203, spherical-harmonics, decomposition
c***author             schneider, b. i.(lanl/nsf)
c***source             m6202
c***purpose            decompose a three dimensional function in 
c***                   spherical harmonics.
c***description        a function is decomposed into spherical harmonics
c***                   using angular quadrature.
c***references         none 
c
c***routines called    iosys, util and mdutil
c***end prologue       m6202
      program ftoylm
      implicit integer (a-z)
      parameter ( dimcen=10 , dimshl=100, ngauss=100, acc=30 )
      character*3 yn, itoc, defnrm
      character*3 nyy
      character*10 cpass
      character*30 str, rtyp, thtyp, phtyp
      character*128 fillam, str1
      character*4096 ops
      character*80 title
      character*8 chrkey, type
      character*800 card
      logical logkey, prnty, prntlm, tstnrm
      logical prntfn, prntg, nonsep, ovride, spln
      real*8 z, znuc, center, alpha, tmp, k, fpkey
      dimension znuc(dimcen), center(3,dimcen), nr(dimshl)
      dimension alpha(ngauss), nx(ngauss), ny(ngauss), nz(ngauss)
      dimension nc(ngauss)
      dimension z(1)
      common a(1)
      common /io/ inp, iout
      common /memory / ioff
      equivalence (z,a)
      call drum
      call iosys ('read character options from rwf',-1,0,0,ops)
c----------------------------------------------------------------------c
c          read in options to print, or test normalization of          c
c          spherical harmonics.                                        c
c----------------------------------------------------------------------c
      prnty=logkey(ops,'print=m6202=ylm',.false.,' ')
      prntlm=logkey(ops,'print=m6202=lm-decomposition',.false.,' ')
      prntg=logkey(ops,'print=m6202=grid',.false.,' ')
      prntfn=logkey(ops,'print=m6202=function',.false.,' ')
      spln=logkey(ops,'m6202=no-spline',.false.,' ')
      order=intkey(ops,'m6202=order-of-spline-fit',4,' ')
      type=chrkey(ops,'m6202=function-type','gaussian',' ')
      tstnrm=logkey(ops,'m6202=test-normalization',.false.,' ')
      call iosys ('read character "linear algebraic filename" from rwf',
     1                 -1,0,0,fillam)
      call iosys ('open lamdat as old',0,0,0,fillam)
      defnrm='off'
      if (tstnrm) then
          defnrm='on'
      endif      
c----------------------------------------------------------------------c
c            read in basic information on centers, points,             c
c            angular integration type, l and m values                  c
c----------------------------------------------------------------------c
      write (iout,1)
      call iosys ('read integer "max number of shells" from '//
     1            'lamdat',1,maxshl,0,' ')
      if (maxshl.gt.dimshl) then
          call lnkerr('quit. number of shells too big')
      endif    
      call iosys ('read integer "number of atomic centers" from lamdat',
     1            1,ncen,0,' ')
      call iosys('read character "scattering center" from lamdat',-1,0,
     1            0,yn)
      call iosys ('read real "nuclear charges" from lamdat',ncen,
     1             znuc,0,' ')
      do 10 i=1,ncen
         call iosys ('read real "atomic center positions" from '//
     1               'lamdat without rewinding',3,center(1,i),0,' ')
   10 continue
      call iosys ('read integer "max theta points" from lamdat',1,
     1             nthmax,0,' ')
      call iosys ('read integer "max phi points" from lamdat',1,
     1             nphmax,0,' ')
      call iosys('read character "non separable quadrature" '//
     1           'from lamdat',-1,0,0,nyy)
      nonsep=.false. 
      if (nyy.eq.'yes') then
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
      call iosys ('read integer "max grid points" from lamdat',1,
     1             totpt,0,' ')
      call iosys('read integer "max radial points in shell" '//
     1           'from lamdat',1,nrmax,0,' ')
      call iosys('read integer "max radial points" from lamdat',1,
     1            totr,0,' ')
      write (iout,2) ncen, ncplus, nrmax, totr, nthmax, nphmax, totpt
      if (nonsep) then
          write(iout,3) maxleb
      endif
      call iosys ('read integer "biggest l" from lamdat',1,
     1             lmax,0,' ')
      call iosys ('read integer "biggest m" from lamdat',1,
     1             mumax,0,' ')
      ovride=logkey(ops,'m6202=override-angular-momentum',.false.,' ')
      if (ovride) then
          lmax=intkey(ops,'m6202=lmax',lmax,' ')
          mumax=intkey(ops,'m6202=mmax',mumax,' ')
      endif    
      lplus=lmax+1
      mplus=mumax+1 
      write (iout,4) mumax, lmax
c----------------------------------------------------------------------c
c                     read in information on functions                 c       c                                                                      c
c----------------------------------------------------------------------c
      if (type.eq.'gaussian') then     
          call posinp('$gaussian',cpass)
          call cardin(card)
          nfun=intkey(card,'number-of-gaussians',1,' ')
          call intarr(card,'x-powers',nx,nfun,' ')
          call intarr(card,'y-powers',ny,nfun,' ')
          call intarr(card,'z-powers',nz,nfun,' ')
          call intarr(card,'position',nc,nfun,' ')
          call fparr(card,'exponents',alpha,nfun,' ')
      elseif (type.eq.'bessel') then
          call posinp('$bessel',cpass)
          call cardin(card)
          nfun=1
          lbes=intkey(card,'l-bessel',0,' ')
          k=fpkey(card,'k-value',.1d0,' ')
          tmp=lbes*acc
          ltop=lbes+sqrt(tmp)
          ltop=max(ltop,lbes)
          ltop1=ltop+1
      else
          call lnkerr('error in function type')
      endif
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
c         memory allocation for gaussian, and legendre decomposed functions.
c          
          fun=wtrad+wtdim+wtdim
          flm=fun+2*totpt
          if (type.eq.'bessel') then
              xval=flm+totr*tlsize
              if (nonsep) then
                  j=xval+maxleb
              else
                  j=xval+nthmax*nphmax
              endif
              add=totr*ltop1
          else
              j=flm
              add=totr*tlsize
          endif
          scr=j+add
          words=wpadti(scr+maxscr)
          if (.not.spln) then
              sc=scr+maxscr
              xtmp=sc+nsplw
              break=xtmp+totr
              c=break+nbreak+1
              ind=wpadti(c+order*nbreak)
              words=ind+totr
          endif
          if (i.eq.1) then
              call iosys ('read integer maxsiz from rwf',1,
     1                     canget,0,' ')
              if (words.gt.canget) then
                  call lnkerr('not enough memory. will quit')
              endif
              call iosys ('write integer maxsiz to rwf',1,
     1                     words,0,' ')
              call getscm(words,z,ngot,'m6202',0)
              write (iout,*) 'get ',words,' words of core'      
          endif
   30 continue
c----------------------------------------------------------------------c 
c             calculate factorials and double factorials needed        c
c             for computation of legendre functions.                   c   
c----------------------------------------------------------------------c 
      call fact(z(dfct),z(ddfct),maxfac)
c----------------------------------------------------------------------c
      cpass='scattering'
      str='scattering center'
c
c        we branch here depending on the type of angular quadrature
c
      if (nonsep) then
          call iosys ('read integer "number lebedev points '//
     1                 str//'" from lamdat',1,nang,0,' ')     
          call iosys ('read real "lebedev angular points '//
     1                 str//'" from lamdat',2*nang,
     1                 z(angle),0,' ')
          call iosys ('read real "lebedev angular weights '//
     1                 str//'" from lamdat',nang,
     2                 z(wtang),0,' ')
          nthet=nang
          nphi=nang
      else                
          call iosys('read integer "theta quadrature order '//
     1                str//'" from lamdat',1,nthet,0,' ')  
          call iosys('read integer "phi quadrature order '//
     1                str//'" from lamdat',1,nphi,0,' ') 
          call iosys('read character "theta quadrature type '//
     1                str//'" from lamdat',-1,0,0,thtyp)
          call iosys('read character "phi quadrature type '//
     1                str//'" from lamdat',-1,0,0,phtyp)
          call iosys ('read real "theta points '//str//
     1                '" from lamdat',nthet,z(theta),0,' ')
          call iosys ('read real "theta weights '//str//
     1                '" from lamdat',nthet,z(wtthe),0,' ')
          call iosys ('read real "phi points '//str//
     1                '" from lamdat',nphi,z(phi),0,' ')
          call iosys ('read real "phi weights '//str//
     1                '" from lamdat',nphi,z(wtphi),0,' ')
      endif     
      call iosys('read integer "max l value '//
     1            str//'" from lamdat',1,lval,0,' ')
      call iosys ('read integer "max m value '//
     1            str//'" from lamdat',1,mval,0,' ')
      if (ovride) then
          lval=lmax
          mval=mumax 
      endif         
c        
c        read radial shell information and associated points
c
      call iosys ('read integer "number of shells '//
     1             str//'" from lamdat',1,nshell,0,' ')
      call iosys ('read integer "number radial points per '//
     1            'shell '//str//'" from lamdat',nshell,nr,0,' ')
      call iosys('read character "radial quadrature type '//
     1            str//'" from lamdat',-1,0,0,rtyp)
      locr=rad
      locwt=wtrad
      nrtot=0
      nwgts=0
      locru=xtmp
      do 50 ns=1,nshell
         call iosys ('read real "radial points '//
     1                str//'" from lamdat without rewinding',
     2                nr(ns),z(locr),0,' ')
        nwts=nr(ns)
        if (rtyp.eq.'newton-cotes') then
            nwts=nr(ns)*(nr(ns)-1)
            nwgts=nwgts+nwts
            if (.not.spln) then
                 call runiqe(z(locr),z(locru),rtyp,nr(ns))
                 locru=locru+nr(ns)-1
            endif     
        else
            if (.not.spln) then
                 call runiqe(z(locr),z(locru),rtyp,nr(ns))
                 locru=locru+nr(ns)
            endif     
            nwgts=nwgts+nwts
        endif
        call iosys ('read real "unscaled radial weights '//
     1               str//'" from lamdat without '//
     2              'rewinding',nwts,z(locwt),0,' ')
        locr=locr+nr(ns)
        locwt=locwt+nwts
        nrtot=nrtot+nr(ns)
   50 continue
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
c----------------------------------------------------------------------c
c              compute the spherical harmonics                         c
c----------------------------------------------------------------------c     
         call legend(z(locplm),z(theta),z(dfct),z(ddfct),nthet,
     1               lval,m,maxfac)
         call scmphi(z(phi),z(locphm),nphi,m)
         title='p(l,m) functions for m = '//itoc(m)
         if (prnty) then
             call prntfm(title,z(locplm),nthet,rowlen,nthet,
     1                   rowlen,iout)
         endif
         title='phi angular functions for m = '//itoc(m)
         if (prnty) then
             call prntfm(title,z(locphm),nphi,nm,nphi,nm,iout)
         endif
c----------------------------------------------------------------------c
c               form spherical harmonics if required                   c
c----------------------------------------------------------------------c
         if (nonsep) then
             call mkylm(z(locplm),z(locphm),z(locylm),nang,
     1                  lval,m)
         endif                   
         totplm=totplm+rowlen
         totylm=totylm+nm*rowlen
         totphm=totphm+nm
         locplm=locplm+wds
         locphm=locphm+wds1
         locylm=locylm+rowlen*nang*nm
   60 continue
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
     1                z(wtphi),z(wtang),z(scr),nthet,nphi,
     2                nang,lval,mval,totplm,totylm,totphm,nonsep)
      endif
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c           calculate the spline information only dependent on the     c
c           grid used not on the fitted function.                      c
c----------------------------------------------------------------------c
      if (.not.spln) then
          call prespl(totr,z(xtmp),order,z(sc))
          call splmat(totr,z(xtmp),order,z(sc))
          nbk=totr-order+1
          nspls=(4*totr+1)*order
          str1='"number of break points '//str//'"'
          call iosys ('write integer '//str1//' to lamdat',1,
     1                 nbk,0,' ')
      endif
      if (prntg) then    
          title='radial grid'
          call prntfm(title,z(rad),nrtot,1,nrtot,1,iout)
      endif    
c----------------------------------------------------------------------c
c             loop over the gaussians                                  c
c----------------------------------------------------------------------c
      do 70 i=1,nfun
         if (type.eq.'gaussian') then
             write(iout,5) nx(i),ny(i),nz(i),alpha(i),
     1                     (center(ii,nc(i)),ii=1,3)
             if (nonsep) then
                 call gnns(z(fun),nx(i),ny(i),nz(i),alpha(i),
     1                     center(1,nc(i)),z(rad),z(theta),z(cphi),
     2                     z(sphi),z(scr),nshell,nrtot,nr,nthet,
     3                     nphi,nang)
             else
                 call gns(z(fun),nx(i),ny(i),nz(i),alpha(i),
     1                    center(1,nc(i)),z(rad),z(theta),z(cphi),
     2                    z(sphi),z(scr),nshell,nrtot,nr,nthet,
     3                    nphi)
             endif
         else
             if (nonsep) then
                 call besns(z(fun),k,lbes,z(j),z(xval),
     1                      center(1,nc(i)),z(rad),z(theta),
     2                      z(cphi),z(sphi),z(scr),nshell,
     3                      nrtot,nr,nthet,nphi,nang,ltop)
             else
                 call bess(z(fun),k,lbes,z(j),z(xval),
     1                     center(1,nc(i)),z(rad),z(theta),
     2                     z(cphi),z(sphi),z(scr),nshell,
     3                     nrtot,nr,nthet,nphi,nang,ltop)
             endif
         endif                                                           
c----------------------------------------------------------------------c
c            calculate inverse of weights                              c
c----------------------------------------------------------------------c
         if (nonsep) then
             call vinv(z(wtang+nang),z(wtang+nang),nang)
             call copy(z(fun),z(fun+totpt),nrtot*nang)
         else
             call vinv(z(wtthe+nthet),z(wtthe+nthet),nthet)
             call vinv(z(wtphi+nphi),z(wtphi+nphi),nphi)
             call copy(z(fun),z(fun+totpt),nrtot*nthet*nphi)
         endif
c----------------------------------------------------------------------c
c             decompose gaussian into spherical harmonics.             c 
c----------------------------------------------------------------------c
         call lmdcmp(z(fun),z(plm),z(phifn),z(ylm),z(flm),
     1               z(wtthe),z(wtphi),z(wtang),z(scr),lval,mval,
     2               nrtot,nthet,nphi,nang,nonsep,prntlm)
         call plmtof(z(fun),z(fun),z(plm),z(phifn),z(ylm),
     1               z(flm),z(scr),z(rad),lval,mval,nrtot,nthet,
     2               nphi,nang,nonsep,prntlm)
         call compar(z(fun+totpt),z(fun+totpt),z(fun),z(fun),nrtot,
     1               nthet,nphi,nang,nonsep)
c----------------------------------------------------------------------c
c               get the spline coefficients                            c
c----------------------------------------------------------------------c
         if (.not.spln) then
             str1='"spline coefficients '//str//'"'
             wrds=(mval+1)*(lval+1)*order*nbk 
             call iosys ('create real '//str1//' on lamdat',wrds,
     1                    0,0,' ')
             call mkspln(z(flm),z(break),z(c),z(sc),mval,lval,
     1                   nrtot,order,nbk,nspls,str1)
             call caspln(z(flm),z(xtmp),z(break),z(c),a(ind),
     1                   mval,lval,nrtot,order,nbk,nspls,str1,
     2                   prntlm)
         endif  
   70 continue
      call iosys ('write integer maxsiz to rwf',1,canget,0,' ')
      call chainx(0)
      stop
    1 format (//,20x,'***** m6202: decompose functions in spherical harm
     1onics *****')
    2 format(/,5x,'atomic data:',
     1      /,5x, 'no. atomic centers = ',i3,3x,'total no. centers = ',
     2      i3,/,5x,'max. radial points in shell = ',i4,3x,
     3              'max. radial points = ',i4,
     4        /5x,  'max. theta points = ',i4,3x,'max. phi points = ',
     5       i3,/,5x,  'max. grid points = ',i6)         
    3 format (/,5x,'using lebedev quadrature. max. lebedev points = ',
     1              i3)
    4 format (/,5x,'maximum m value',1x,i3,5x,'maximum l value',1x,i3)
    5 format(/,1x,'gaussian function with:'/,1x,
     1            'nx = ',i2,1x,'ny = ',i2,1x,'nz = ',i2,/,1x,
     2            'exponent = ',f12.5,1x,'center = (',f12.5,','f12.5,
     3                                              ',',f12.5,')')        
      end
















