*deck orthfn.f 
c***begin prologue     orthfn
c***date written       951230   (yymmdd)
c***revision date               (yymmdd)
c***keywords           orthogonal functions, gauss quadrature
c***author             schneider, b. i.(nsf)
c***source             orthfn
c***purpose            driver for generating orthogonal functionss
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       orthfn
      program orthfn
c
      implicit integer (a-z)
      character*4096 ops
      character*80 cpass
      character*80 title, chrkey
      character*1600 card
      character*16 type
      character*2 itoc
      logical dollar, logkey, prncof, prnply, prnwpt, check
      logical dvrfun, fixed
      logical fixl, fixr 
      common z(1)
      dimension ia(1), endpts(2), der(2), temp(2), lval(maxl)
      dimension pleft(maxl), pright(maxl)
      equivalence (ia(1),z(1))
      common/io/inp, iout      
      common/memory/ioff
      real*8 z, left, right, fpkey, endpts, alpha, beta
      real*8 muzero
c
      call drum
      write(iout,*)
      write(iout,*) '          orthogonal functions'
      call iosys ('read character options from rwf',-1,0,0,ops)
      prncof=logkey(ops,'print=m6300=recursion-coefficients',
     1              .false.,' ')
      prnply=logkey(ops,'print=m6300=functions',.false.,' ')
      prnwpt=logkey(ops,'print=m6300=points/weights',.false.,' ')
      check=logkey(ops,'check-orthogonality',.false.,' ')
      decomp=logkey(ops,'decompose',.false.,' ')
      dvrfun=logkey(ops,'dvr-functions',.false.,' ')
      if ( dollar('$orthfn',card,cpass,inp) ) then
           call cardin(card)
           type=chrkey(card,'type-orthogonal-functions','legendre',' ')
           if(type.eq.'jacobi'.or.type.eq.'laguerre') then
              alpha=fpkey(card,'alpha',0.d0,' ')
              beta=fpkey(card,'beta',0.d0,' ')
           endif              
           left=fpkey(card,'left-boundary',-1.d0,' ')
           right=fpkey(card,'right-boundary',1.d0,' ')
           npts=intkey(card,'quadrature-order',10,' ')
           fixed=logkey(card,'fix-end-points',.false.,' ')
           write(iout,1) type, npts, left, right            
           nfixed=0
           if (fixed) then
               fixl=logkey(card,'fix-left',.false.,' ')
               fixr=logkey(card,'fix-right',.false.,' ')
           endif
           if(fixl) then
              nfixed=nfixed+1
              write(iout,2) 
           endif
           if(fixr) then
              nfixed=nfixed+1
              write(iout,3) 
           endif
      endif
      if(nfixed.eq.1) then
         if (fixl) then
             endpts(1) = left
         endif
         if (fixr) then
             endpts(1) = right
         endif
      endif
      if(nfixed.eq.2) then
         endpts(1)=left
         endpts(2)=right
      endif              
      ioff=1
      do 10 i=1,2
         if(genmat.eq.'full-calculation'.or.
     1      genmat.eq.'matrix-generation-only') then
            x=ioff
            xwts=x+dim
            dxwts=xwts+dim
            ddxwts=dxwts+dim
            wts=ddxwts+dim
            a=wts+dim
            b=a+nmax+1
            p=b+nmax+1
            dp=p+nmax*dim
            ddp=dp+nmax*dim
            scr=ddp+nmax*dim
            wds0=scr
            if(check.or.decomp.or.diffeq) then
               dum=wds0+dim*dim
               f=dum+max(3*dim,dim*dim)
               df=f+dim
               ddf=df+dim
               ham=ddf+dim
               eig=ham+nmax*nmax
               work=eig+dim+1
               dum1=work+max(3*(dim+1),dim*dim)
               ipvt=wpadti(dum1+dim+1)
               wds0=iadtwp(ipvt+dim)
               if(coord) then
                  pn=wds0
                  dpn=pn+nmax*dim
                  ddpn=dpn+nmax*dim
                  newwts=ddpn+nmax*dim
                  mat=newwts+dim
                  work1=mat+dim
                  wds0=work1+dim
               endif
               if(itdiag.eq.'davidson') then
                  vec=wds0
                  hvec=vec+nmax*maxvec
                  root=hvec+nmax*maxvec
                  dvdmat=root+nmax
                  dvdvec=dvdmat+maxvec*maxvec
                  cvn=dvdvec+maxvec*maxvec
                  hmev=cvn+nmax
                  wds0=hmev+nmax*maxvec
               endif               
            endif
         elseif(genmat.eq.'diagonalize') then
            ham=ioff
            eig=ham+nmax*nmax
            work=eig+nmax
            dum1=work+nmax
            wds0=dum1+nmax
            if(itdiag.eq.'davidson') then
               vec=wds0
               hvec=vec+nmax*maxvec
               root=hvec+nmax*maxvec
               dvdmat=root+nmax
               dvdvec=dvdmat+maxvec*maxvec
               cvn=dvdvec+maxvec*maxvec
               hmev=cvn+nmax
               work=hmev+nmax*maxvec
               wds0=work+nmax
            endif
         else
            call lnkerr('error in task assignment')            
         endif                                              
         words=wpadti(wds0)
         if (i.eq.1) then
             call iosys ('read integer maxsiz from rwf',1,
     1                    canget,0,' ')
             if (words.gt.canget) then
                 call lnkerr('not enough memory. will quit')
             endif
             call iosys ('write integer maxsiz to rwf',1,
     1                    words,0,' ')
         else
             call getscm(words,z,ngot,'orthfn',0)
         endif
   10 continue
c
c     get the recursion coefficients for the polynomials.
c
      call class(type,npts,alpha,beta,z(b),z(a),muzero)
c     
c           convert to points and weights on (left,right) if needed
c
               call convt(z(x),z(wts),left,right,norm0,npts,
     1                    type,prnwpt)
c      
            else
               call gaussq(type,npts,alpha,beta,nfixed,endpts,z(scr),
     1                               z(x),z(wts))
            endif
c
c           calculate the recursion coefficients
c
            if (prncof) then
                title='a coefficients'
                call prntrm(title,z(a+1),nmax,1,nmax,1,iout)
                title='b coefficients'
                call prntrm(title,z(b+1),nmax-1,1,nmax-1,1,iout)
            endif
c           calculate the orthonormal functions and their first and 
c                        second derivatives.
c

                  call iosys('destroy "second derivative of '//
     1                       'orthogonal functions for l='
     2                       //itoc(lval(l))//'"',0,0,0,' ')
               endif         
               call iosys('write real "orthogonal functions '//
     1                    'for l='//itoc(lval(l))//'" to lamdat',
     2                     npts*nmax,z(p),0,' ')
               call iosys('write real "first derivative of '//
     1                    'orthogonal functions for l='
     2                     //itoc(lval(l))//'" to lamdat',
     2                     npts*nmax,z(dp),0,' ')
               call iosys('write real "second derivative of '//
     1                    'orthogonal functions for l='
     2                     //itoc(lval(l))//'" to lamdat',
     2                     npts*nmax,z(ddp),0,' ')
            endif
            if(coord) then
               locply=pn
               call diagx(z(p),z(dp),z(ddp),z(a+1),z(b+1),z(pn),
     1                    z(dpn),z(ddpn),z(eig),z(work),z(dum),
     2                    nmax,npts,prnh)
            if (prnply) then
                title='dvr orthfnnomials'
                call prntrm(title,z(pn),npts,nmax,npts,nmax,iout)
                title='first derivative of dvr orthfnnomials'
                call prntrm(title,z(dpn),npts,nmax,npts,nmax,iout)
                title='second derivative of dvr orthfnnomials'
                call prntrm(title,z(ddpn),npts,nmax,npts,nmax,iout)
            endif
               call iosys('does "x-eigenvalues for l='
     1                    //itoc(lval(l))//'" exist on '//
     2                    'lamdat',0,0,0,yn)
               if(yn.eq.'yes') then
                  call iosys('destroy "x-eigenvalues '//
     1                       'for l='//itoc(lval(l))//'"',0,0,0,' ') 
                  call iosys('destroy "x-diagonal functions for l='
     2                       //itoc(lval(l))//'"',0,0,0,' ')
                  call iosys('destroy "first derivative of '//
     1                       'x-diagonal functions for l='
     2                       //itoc(lval(l))//'"',0,0,0,' ')  
                  call iosys('destroy "second derivative of '//
     1                       'x-diagonal functions for l='
     2                       //itoc(lval(l))//'"',0,0,0,' ')
               endif              
               call iosys('write real "x-eigenvalues for l='
     1                     //itoc(lval(l))//'" to lamdat',nmax,
     2                       z(eig),0,' ')
               call iosys('write real "x-diagonal functions for l='
     1                    //itoc(lval(l))//'" to lamdat',npts*nmax,
     2                   z(pn),0,' ')
               call iosys('write real "first derivative of '//
     1                    'x-diagonal functions for l='
     2                    //itoc(lval(l))//'" to lamdat',
     3                     npts*nmax,z(dpn),0,' ')
               call iosys('write real "second derivative of '//
     1                    'x-diagonal functions for l='
     2                    //itoc(lval(l))//'" to lamdat',
     2                     npts*nmax,z(ddpn),0,' ')
               call mkwts(z(pn),z(wts),z(newwts),nmax,npts)
               call matx(z(p),z(pn),z(x),z(wts),z(f),z(mat),z(newwts),
     1                   z(dum),nmax,npts)
               if(check) then
                  call chk(z(pn),z(wts),z(scr),z(dum),nmax,npts)
               endif   
               if(decomp) then
                   call fcoef(z(f),z(df),z(ddf),z(pn),z(dpn),z(ddpn),
     1                       z(x),z(wts),z(scr),z(dum),nmax,
     2                       npts,nocoef)
               endif 
               if(diffeq) then
                  call diff(z(pn),z(dpn),z(ddpn),z(x),z(wts),z(ham),
     1                      z(eig),der,lval(l),nmax,npts,pottyp,prange,
     2                      coord,prnh)
                  if(genmat.ne.'matrix-generation-only') then
                     if(itdiag.eq.'false') then
                        nkept=nmax
                        call diag(z(ham),z(eig),z(work),z(dum1),nmax,
     1                            diatyp,prnh)
                     elseif(itdiag.eq.'lanczos') then
                        call mktrid(z(ham),z(eig),z(work),z(dum1),z(pn),
     1                              nmax,niter,prnh)
                     elseif(itdiag.eq.'davidson') then
                        nkept=nroots
                        call david(z(ham),z(eig),z(vec),z(hvec),z(root),
     1                             z(dvdmat),z(dvdvec),z(cvn),z(hmev),
     2                             z(work),ia(ipvt),thresh,cnverg,
     3                             ops,nmax,nroots,nattim,maxvec,
     4                             niter,prgues,prvec,prtflg,prres)
                     else
                        call lnkerr('error in diagonalization '//
     1                              'technique')
                     endif
                  else
                     write(iout,8)    
                  endif
                  call iosys('does "hamiltonian matrix for l='
     1                       //itoc(lval(l))//'" exist on '//
     2                       'lamdat',0,0,0,yn)
                  if(yn.eq.'yes') then
                     call iosys('destroy "hamiltonian matrix '//
     1                          'for l='//itoc(lval(l))//'"',0,0,0,' ')
                  endif                                          
                  call iosys('write real "hamiltonian matrix for l='
     1                       //itoc(lval(l))//'" to lamdat',
     2                         nmax*nmax,z(ham),0,' ')
               endif
            else
               locply=p      
c           check orthogonality on discrete interval
               if(check) then
                  call chk(z(p),z(wts),z(scr),z(dum),nmax,npts)
               endif
               if(decomp) then
                  call fcoef(z(f),z(df),z(ddf),z(p),z(dp),z(ddp),z(x),
     1                       z(wts),z(scr),z(dum),nmax,npts,nocoef)
               endif
               if(diffeq) then
                  call diff(z(p),z(dp),z(ddp),z(x),z(wts),z(ham),
     1                      z(eig),der,lval(l),nmax,npts,pottyp,prange,
     2                      coord,prnh)
                  if(genmat.ne.'matrix-generation-only') then
                     if(itdiag.eq.'false') then
                        nkept=nmax
                        call diag(z(ham),z(eig),z(work),z(dum1),nmax,
     1                            diatyp,prnh)
                     elseif(itdiag.eq.'lanczos') then
                        call mktrid(z(ham),z(eig),z(work),z(dum1),z(p),
     1                              nmax,niter,prnh)
                     elseif(itdiag.eq.'davidson') then
                        nkept=nroots
                        call david(z(ham),z(eig),z(vec),z(hvec),z(root),
     1                             z(dvdmat),z(dvdvec),z(cvn),z(hmev),
     2                             z(work),ia(ipvt),thresh,cnverg,
     3                             ops,nmax,nroots,nattim,maxvec,
     4                             niter,prgues,prvec,prtflg,prres)
                     else
                        call lnkerr('error in diagonalization '//
     1                              'technique')
                     endif
                  else
                     write(iout,9)   
                  endif
                  call iosys('does "hamiltonian matrix for l='
     1                       //itoc(lval(l))//'" exist on '//
     2                       'lamdat',0,0,0,yn)
                  if(yn.eq.'yes') then
                     call iosys('destroy "hamiltonian matrix '//
     1                       'for l='//itoc(lval(l))//'"',0,0,0,' ')
                  endif
                  call iosys('write real "hamiltonian matrix for l='
     1                       //itoc(lval(l))//'" to lamdat',
     2                       nmax*nmax,z(ham),0,' ')
               endif
            endif
         else
            write(iout,9)
            call iosys('read real "hamiltonian matrix for l='
     1                       //itoc(lval(l))//'" from lamdat',
     2                         nmax*nmax,z(ham),0,' ')
            if(itdiag.eq.'false') then
                   nkept=nmax
                   call diag(z(ham),z(eig),z(work),z(dum1),nmax,
     1                       diatyp,prnh)
            elseif(itdiag.eq.'lanczos') then
                   call mktrid(z(ham),z(eig),z(work),z(dum1),z(vec),
     1                         nmax,niter,prnh)
            elseif(itdiag.eq.'davidson') then
                   nkept=nroots
                   call david(z(ham),z(eig),z(vec),z(hvec),z(root),
     1                        z(dvdmat),z(dvdvec),z(cvn),z(hmev),
     2                        z(work),ia(ipvt),thresh,cnverg,
     3                        ops,nmax,nroots,nattim,maxvec,niter,
     4                        prgues,prvec,prtflg,prres)
            else
                   call lnkerr('error in diagonalization technique')
            endif 
         endif
 500  continue
      if(doscat) then
         if ( posinp('$energy',cpass) ) then         
              call cardin(card)
         endif     
              nen=intkey(card,'number-of-energies',1,' ')
              call trnply(z(locply),z(ham),z(work),z(dum1),npts,
     1                    nmax,nkept)
              do 600 ien=1,nen
                 read(inp,*) energy
                 call phase(z(work),z(eig),energy,right,nkept)
  600         continue
      endif                                         
      call iosys ('write integer maxsiz to rwf',1,canget,0,' ')
      call chainx(0)               
      stop
 1    format(/,5x,'type function        = ',a16,/,5x,
     1            'order of polynomials = ',i3,/,5x,
     2            'left boundary        = ',e15.8,/,5x,
     3            'right boundary       = ',e15.8)
 2    format(/,5x,'fixed left point')
 3    format(/,5x,'fix right point')
 4    format(/,5x,'for orthfnnomial type  = ',a16,/,5x,
     1            'alpha = ',e15.8,/,5x,
     2            'beta = ',e15.8)
 5    format(/,5x,'for orthfnnomial type  = ',a16,/,5x,
     1            'alpha = ',e15.8)
 6    format(/,15x,'using co-ordinate representation')
 7    format(/,15x,'beginning calculation for angular momentum = ',
     1                                                    i3,/,5x,
     2             'order of leading left orthfnnomial  = ',i2,/5x,
     3             'order of leading right orthfnnomial = ',i2,/5x,)
 8    format(/,15,'matrix generated')     
 9    format(/,15x,'matrix already calculated-only diagonalization')
      end

