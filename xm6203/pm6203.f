*deck @(#)pm6203.f	1.2  10/27/94
      subroutine pm6203(z,a) 
c***begin prologue     pm6203
c***date written       930810  (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6203, link 6203, spherical harmonics
c***author             schneider, b. i.(lanl/nsf)
c***source             m6203
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
c***end prologue       pm6203.f
      implicit integer (a-z)
      parameter ( dimcen=10 , dimshl=100, dimene=50, acc=30)
      character*4 chvr1
      character*3 yn, chra, itoc, defp, defnt, defnrm, defdsk, defspl
      character*3 ny
      character*8 cpass
      character*16 fptoc
      character*30 str, rtyp, thtyp, phtyp
      character*128 fillam, str1, str2
      character*1600 card
      character*4096 ops
      character*80 title
      logical logkey, prnty, prntlm, positn, tstnrm, tstint, ylmon
      logical prntyu, noscal, prntnt, prntg, nonsep
      logical spln, onlysp
      integer cskipb
      integer a(1)
      real*8 z, energy, znuc, eta, center, tmp, k
      dimension znuc(dimcen), eta(dimcen), center(3,dimcen), nr(dimshl)
      dimension z(1), energy(dimene)
      common /io/ inp, iout
c
      call iosys ('read character options from rwf',-1,0,0,ops)
c----------------------------------------------------------------------c
c          read in options to print, test normalization of             c
c          spherical harmonics, spline fit or test integral            c
c          calculation                                                 c
c----------------------------------------------------------------------c
      onlysp=logkey(ops,'only-spherical-harmonics',.false.,' ')
      prnty=logkey(ops,'print=m6203=ylm',.false.,' ')
      prntlm=logkey(ops,'print=m6203=lm-decomposition',.false.,' ')
      prntyu=logkey(ops,'print=m6203=yukawa',.false.,' ')
      noscal=logkey(ops,'bare-yukawa',.false.,' ')
      prntg=logkey(ops,'print=m6203=grid',.false.,' ')
crlm     call iosys ('read character "linear algebraic filename" from rwf',
crlm    1                 -1,0,0,fillam)
      call iosys ('read character "grid filename" from rwf',
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
      ylmon=logkey(ops,'ylm=on-disk',.false.,' ')
      defdsk='off'
      if (ylmon) then
          defdsk='on'
      endif    
      spln=logkey(ops,'no-spline',.false.,' ')
      defspl='on'
      if (spln) then
          defspl='off'
      else    
          order=intkey(ops,'order-of-spline-fit',4,' ')
          call iosys ('write integer "order of spline fit" '//
     1                'to lamdat',1,order,0,' ')
      endif
c----------------------------------------------------------------------c
c            read in basic information on centers, points,             c
c            angular integration type, l and m values, energies        c
c----------------------------------------------------------------------c
      write (iout,1)
      write (iout,2) defp, defnrm, defnt, defdsk, defspl, order 
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
      if (positn('$energy',card,inp)) then
           call cardin(card)
           nen=intkey(card,'number-of-energies',1,' ')
           call fparr(card,'energies',energy,nen,' ')
           write(iout,5) nen, (energy(ii),ii=1,nen)
      endif
      call iosys ('read integer "biggest l" from lamdat',1,
     1             lmax,0,' ')
      lplus=lmax+1    
      call iosys ('read integer "biggest m" from lamdat',1,
     1             mumax,0,' ')
      mplus=mumax+1
      write (iout,6) mumax, lmax
c----------------------------------------------------------------------c
c                     memory allocation                                c
c----------------------------------------------------------------------c
c
c     set up some integers giving maximum values of dimensions and
c     compute the default accuracy for the bessel functions
c
      lsize=0
      do 20 lstor=0,mumax
         lsize=lsize+lmax-lstor+1
   20 continue          
      maxfac=lplus+mplus+10
      maxang=max(nthmax,nphmax)
      lsq=lplus*lplus
      tmp=lmax*acc
      ltop=lmax+sqrt(tmp)
      ltop=max(ltop,lmax)
      ltop1=ltop+1
      wtdim=totr*(totr-1)
      nsplw=(4*totr+1)*order*2
      nbreak=totr-order+1
      maxscr=max(nthmax,2*totr*nthmax,2*totr*lplus)
      if (tstnrm) then
          maxscr=max(maxscr,lsq)
      endif               
      ioff=1
      do 30 i=1,2
c
c         memory allocation for factorials, grid points, weights, sine and
c         cosine, phi functions and legendre polynomials.      
          dfct=ioff
          ddfct=dfct+maxfac+1
          if (nonsep) then
             angle=ddfct+maxfac+1
             theta=angle
             phi=theta+maxleb 
             wtang=angle+maxleb+maxleb
             sphi=wtang+maxleb+maxleb
             cphi=sphi+maxleb
c            the doubling of the dimension for plm allows the storage 
c            for the second half of the y(l,m) array to overlap p(l,m).  
             ylm=cphi+maxleb
crlm         plm=ylm+maxleb*lsize
             plm=ylm+maxleb*lsize*2
             phifn=plm+maxleb*lsize
             rad=phifn+2*mplus*maxleb
          else
             theta=ddfct+maxfac+1
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
              y=j+totr*ltop1
              wron=y+totr*ltop1
              flm=wron+ltop1
              psilm=flm+totr*lsize
c
c         starting at jp storage can be overlapped.
c          
              jp=psilm+totr*lsize
              yp=jp+totr*ltop1
              scr=yp+totr*ltop1
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
              call getscm(words,z,ngot,'m6203',0)
crlm             call iosys ('read integer maxsiz from rwf',1,
crlm     1                     canget,0,' ')
              if (words.gt.ngot) then
                  call lnkerr('not enough memory. will quit')
              endif
crlm              call iosys ('write integer maxsiz to rwf',1,
crlm     1                     words,0,' ')
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
             len=cskipb(chra,' ')
             str='atom-'//chra(1:len)
         else
             if (yn.ne.'yes') then
                 cpass='atom'
                 chra=itoc(atom)
                 len=cskipb(chra,' ')
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
             if(.not.spln) locru=xtmp
             do 50 ns=1,nshell
                call iosys ('read real "radial points '//
     1                      str//'" from lamdat without rewinding',
     2                      nr(ns),z(locr),0,' ')
                nwts=nr(ns)
                if (rtyp.eq.'newton-cotes') then
                    nwts=nr(ns)*(nr(ns)-1)
                    nints=nints+nr(ns)-1
                    if(.not.spln) then
                       call runiqe(z(locr),z(locru),rtyp,nr(ns))
                       locru=locru+nr(ns)-1
                    endif
                else
                    if(.not.spln) then
                       call runiqe(z(locr),z(locru),rtyp,nr(ns))
                       locru=locru+nr(ns)
                    endif
                    nints=nints+nr(ns)
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
             call iosys ('create real '//str1//' on lamdat',
     1                    nwds,0,0,' ')
             nwds=(2*mval+1)*nphi
             call iosys ('create real '//str2//' on lamdat',
     1                    nwds,0,0,' ')
         endif
         call iosys ('rewind '//str1//' on lamdat '//
     1               'read-and-write',0,0,0,' ')
         call iosys ('rewind '//str2//' on lamdat '//
     1               'read-and-write',0,0,0,' ')
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
            wds1=nphi*nm
            if (.not.ylmon) then
c----------------------------------------------------------------------c
c              compute the spherical harmonics                         c
c----------------------------------------------------------------------c     
                call legend(z(locplm),z(theta),z(dfct),z(ddfct),nthet,
     1                      lval,m,maxfac)
                call scmphi(z(cphi),z(sphi),z(locphm),nphi,m)
                call iosys ('write real '//str1//' to lamdat '//
     1                      'without rewinding',wds,z(locplm),0,' ')
                chvr1=itoc(m) 
                len=cskipb(chvr1,' ')
                title='"p(l,'//chvr1(1:len)//')"'
                lentit=cskipb(title,' ')
                write(iout,*) 'writing iosys file ', title(1:lentit)
                write(iout,*) 'words written = ',wds
                title='p(l,m) functions for m = '//itoc(m)
                lentit=cskipb(title,' ')
                if (prnty) then
                    call prntfm(title,z(locplm),nthet,rowlen,nthet,
     1                          rowlen,iout)
                endif
                call iosys ('write real '//str2//' to lamdat '//
     1                      'without rewinding',wds1,z(locphm),0,' ')
                title='"phi('//chvr1(1:len)//')"'
                lentit=cskipb(title,' ')
                write(iout,*) 'writing iosys file ', title(1:lentit)
                write(iout,*) 'words written = ',wds1
                title='phi angular functions for m = '//itoc(m)
                lentit=cskipb(title,' ')
                if (prnty) then
                    call prntfm(title,z(locphm),nphi,nm,nphi,nm,iout)
                endif
c----------------------------------------------------------------------c
c               form spherical harmonics if required                   c
c----------------------------------------------------------------------c
                if (nonsep) then
                    call mkylm(z(locplm),z(locphm),z(locylm),nang,
     1                         lval,m)
                    if(prnty) then
                       call prntfm('ylm',z(locylm),nang,rowlen*nm,
     $                         nang,rowlen*nm,iout)
                    endif
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
            locplm=locplm+ntheta*rowlen
            locphm=locphm+nphi*nm
            locylm=locylm+nang*rowlen*nm
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
     1                    z(wtphi),z(wtang),z(scr),nthet,nphi,
     2                    nang,lval,mval,totplm,totylm,totphm,nonsep)
         endif        
         if (.not.onlysp) then
c----------------------------------------------------------------------c
c             let us assume that the inhomogeneity is a                c
c             product of a multicenter yukawa potential and a          c 
c             bessel function.  this would simulate the first          c
c             iteration of a born series and look like the variation   c
c             iteration method.                                        c
c                                                                      c
c             do the yukawa first.                                     c  
c----------------------------------------------------------------------c
             call yukawa(z(yuk),eta,center,z(rad),z(theta),z(cphi),
     1                   z(sphi),z(scr),ncen,nshell,nrtot,nr,nthet,
     2               nphi,nang,nonsep)
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
c    
c----------------------------------------------------------------------c   
c             begin energy loop                                        c
c----------------------------------------------------------------------c
             do 70 ien=1,nen
                k=sqrt(energy(ien))
                if (k.eq.0.d0) then
                    k=1.d-10
                endif    
c----------------------------------------------------------------------c
c             decompose inhomogeneity into spherical harmonics.        c 
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c             calculate bessel's of kr                                 c
c----------------------------------------------------------------------c 
                call smul(z(rad+nrtot),z(rad),k,nrtot)
                call tstzro(z(rad+nrtot),nrtot)
                call rcbes(z(rad+nrtot),z(j),z(jp),z(y),z(yp),
     1                     z(wron),z(scr),nrtot,lval,ltop,
     2                     'derivatives',.false.)
c----------------------------------------------------------------------c
c          the next bit of monkey business converts derivatives        c
c          wrt k*r to derivatives wrt r, recomputes the wronskian      c
c          which is needed for the green's function and then           c       
c          scales the both functions so that the proper factor         c
c          will be incorporated into the radial integral.              c
c----------------------------------------------------------------------c  
                call sscal(nrtot*lv1,k,z(jp),1)
                call sscal(nrtot*lv1,k,z(yp),1)
                call sscal(lv1,k,z(wron),1)
                call vinv(z(wron),z(wron),lv1)
                call vsqrt(z(wron),z(wron),lv1)
                call mvmul(z(wron),z(j),z(j),nrtot,lv1)
                call mvmul(z(wron),z(y),z(y),nrtot,lv1)
                call mvmul(z(wron),z(jp),z(jp),nrtot,lv1)
                call mvmul(z(wron),z(yp),z(yp),nrtot,lv1)
c----------------------------------------------------------------------c
c             multiply the yukawa in this shell by the zeroth order    c
c             spherical bessel function.                               c
c----------------------------------------------------------------------c
                call scalyu(z(yuk),z(rho),z(j),nrtot,nthet,nphi,nang,
     1                      nonsep,noscal,prntyu)
c----------------------------------------------------------------------c 
c             decompose the function into its angular components       c
c----------------------------------------------------------------------c      
                call lmdcmp(z(rho),z(plm),z(phifn),z(ylm),z(flm),
     1                  z(wtthe),z(wtphi),z(wtang),z(scr),lval,mval,
     2                  nrtot,nthet,nphi,nang,nonsep,prntlm)
                call plmtof(z(rho),z(plm),z(phifn),z(ylm),z(flm),z(scr),
     1                      z(rad),lval,mval,nrtot,nthet,nphi,nang,
     2                      nonsep,prntlm)
c----------------------------------------------------------------------c
c             from the radial components of the inhomogeneity          c
c             calculate the radial component of the solution           c
c             using the green's function to solve the integral         c
c             equation.                                                c
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c               do the integral                                        c
c----------------------------------------------------------------------c
                call drvint(z(flm),z(psilm),z(j),z(y),z(wtrad),z(scr),
     1                      nrtot,lval,mval,nshell,nr,nints,rtyp,prntnt)
                call plmtof(z(rho),z(plm),z(phifn),z(ylm),z(psilm),
     1                      z(scr),z(rad),lval,mval,nints,nthet,nphi,
     2                      nang,nonsep,prntlm)
c----------------------------------------------------------------------c
c               get the spline coefficients                            c
c----------------------------------------------------------------------c
                if (.not.spln) then
                    str1='"spline coefficients '//str//
     1                   ' energy-'//fptoc(energy(ien))//'"'
                    wrds=(mval+1)*(lval+1)*order*nbk 
                    call iosys ('create real '//str1//' on lamdat',wrds,
     1                           0,0,' ')
                    call mkspln(z(psilm),z(break),z(c),z(sc),mval,lval,
     1                          nints,order,nbk,nspls,str1)
                    call caspln(z(psilm),z(xtmp),z(break),z(c),a(ind),
     1                          mval,lval,nints,order,nbk,nspls,str1,
     2                          prntlm)
                endif  
   70        continue
         endif
   40 continue
      call iosys ('write integer maxsiz to rwf',1,canget,0,' ')
      call chainx(0)
      return
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
