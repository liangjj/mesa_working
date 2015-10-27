*deck poly.f 
c***begin prologue     poly
c***date written       951230   (yymmdd)
c***revision date               (yymmdd)
c***keywords           orthogonal polynomial, roots
c***author             schneider, b. i.(nsf)
c***source             poly
c***purpose            driver for generating orthogonal polynomials
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       poly
      program poly
c
      implicit integer (a-z)
      parameter ( maxl=20 )
      character*4096 ops
      character*8 cpass, itdiag, prtflg
      character*80 title, chrkey, genmat
      character*1600 card
      character*16 type, pottyp, wtfn
      character*128 fillam
      character*2 itoc
      character*3 yn
      logical posinp, logkey, prncof, prnply, prnwpt, check, decomp
      logical prnh, diffeq, fixed, coord, diatyp, prgues, prvec
      logical prres, doscat
      common z(1)
      dimension ia(1), endpts(2), der(2), temp(2), lval(maxl)
      dimension pleft(maxl), pright(maxl)
      equivalence (ia(1),z(1))
      common/io/inp, iout      
      common/memory/ioff
      real*8 z, left, right, fpkey, endpts, alpha, beta, der
      real*8 norm0, temp, thresh, cnverg, prange, energy
c
      call drum
      write(iout,*)
      write(iout,*) '          orthogonal polynomials'
      call iosys ('read character options from rwf',-1,0,0,ops)
      call iosys ('read character "print flag" from rwf',-1,0,0,prtflg)
      prncof=logkey(ops,'print=m6245=polynomial-coefficients',
     1              .false.,' ')
      prnply=logkey(ops,'print=m6245=polynomials',.false.,' ')
      prnwpt=logkey(ops,'print=m6245=points/weights',.false.,' ')
      prnh=logkey(ops,'print=m6245=hamiltonian',.false.,' ')
      coord=logkey(ops,'to-x-representation',.false.,' ')
      check=logkey(ops,'check-orthogonality',.false.,' ')
      decomp=logkey(ops,'decompose',.false.,' ')
      diffeq=logkey(ops,'solve-differential-equation',.false.,' ')
      diatyp=logkey(ops,'eigenvalues-and-eigenvectors',.false.,' ')
      itdiag=chrkey(ops,'iterative-diagonalization','false',' ')
      genmat=chrkey(ops,'task','full-calculation',' ')
      doscat=logkey(ops,'scattering-calculation',.false.,' ')
      call iosys ('read character "linear algebraic filename" from rwf',
     1                 -1,0,0,fillam)
      call iosys ('open lamdat as unknown',0,0,0,fillam)
      if ( posinp('$poly',cpass) ) then
           call cardin(card)
           type=chrkey(card,'type-polynomials','standard',' ')
           pottyp=chrkey(card,'potential','none',' ')
           left=fpkey(card,'left-boundary',-1.d0,' ')
           call iosys('write real "left boundary" to lamdat',1,
     1                 left,0,' ')
           right=fpkey(card,'right-boundary',1.d0,' ')
           prange=fpkey(card,'potential-range',right,' ')
           call iosys('write real "right boundary" to lamdat',1,
     1                 right,0,' ')
           der(1)=fpkey(card,'left-derivative',0.d0,' ')
           der(2)=fpkey(card,'right-derivative',0.d0,' ')
           nmax=intkey(card,'order-of-polynomials',10,' ')
           nlval=intkey(card,'number-of-angular-momenta',1,' ')
           npts=intkey(card,'number-of-points',nmax,' ')
           fixed=logkey(card,'fix-end-points',.false.,' ')
           niter=intkey(card,'number-of-iterations',nmax,' ')
           if (fixed) then
               nfixed=intkey(card,'number-of-fixed-points',2,' ')
               call fparr(card,'end-points',endpts,nfixed,' ')
           endif
           alpha=0.d0
           beta=0.d0
           if (type.eq.'jacobi'.or.type.eq.'laguerre') then
               alpha=fpkey(card,'alpha',0.d0,' ')
               beta=fpkey(card,'beta',0.d0,' ')
           endif
           wtfn=chrkey(card,'weight-function','one',' ')
           write(iout,1) type, nmax, left, right, der(1), der(2), 
     1                   nlval, pottyp, wtfn
           if (fixed) then
               if(nfixed.eq.1) then
                  write(iout,2) nfixed, endpts(1)
               endif
               if(nfixed.eq.2) then
                  write(iout,3) nfixed, endpts(1), endpts(2)
               endif 
           endif
           if (type.eq.'jacobi') then
               write(iout,4) type, alpha, beta
           endif
           if (type.eq.'laguerre') then
               write(iout,5) type, alpha
           endif
      endif
      call intarr(card,'angular-momenta',lval,nlval,' ')
      call intarr(card,'order-of-leading-left-polynomials',
     1            pleft,nlval,' ')
      call intarr(card,'order-of-leading-right-polynomials',
     1            pright,maxl,' ')
      if(coord) then
         write(iout,6)
      endif
      dim=max(npts,nmax)
      if(itdiag.eq.'davidson') then
         nroots=intkey(ops,'davidson=number-of-roots',nmax,' ')
         thresh=fpkey(ops,'davidson=tolerance',1.0d-06,' ')
         cnverg=fpkey(ops,'davidson=convergence',1.d-08,' ')
         nattim=intkey(ops,'davidson=number-of-roots-at-a-time',
     1                 nroots,' ')
         maxvec=intkey(ops,'davidson=maximum-number-of-vectors',
     1                 nmax,' ')
         niter=intkey(ops,'davidson=maximum-number-of-iterations',
     1                nmax,' ')
         prgues=logkey(ops,'print=davidson=guess',.false.,' ')
         prvec=logkey(ops,'print=davidson=final-vectors',.false.,' ')
         prres=logkey(ops,'print=davidson=residuals',.false.,' ')
         maxvec=min(maxvec,nmax)
         niter=min(niter,nmax)
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
             call getscm(words,z,ngot,'poly',0)
         endif
   10 continue
      call iosys('write integer "number of polynomials" to lamdat',
     1            1,nmax,0,' ')
      call iosys('write integer "size of hamiltonian '//
     1           'matrix" to lamdat',1,nmax,0,' ')
      do 500 l=1,nlval
         if(genmat.ne.'diagonalize') then
            write(iout,7) lval(l), pleft(l), pright(l)
            call iosys('write integer "order of leading left '//
     1                 'polynomial for l='//itoc(lval(l))//'" '//
     2                 'to lamdat',1,pleft(l),0,' ')
            call iosys('write integer "order of leading right '//
     1                 'polynomial for l='//itoc(lval(l))//'" '//
     2                 'to lamdat',1,pright(l),0,' ')
c
c     in order to proceed it is necessary to be able to calculate integrals
c     between the orthogonal polynomials to be generated on the desired
c     integration interval.  the calculation of the points and weights is
c     done on the interval ( -1.,1. ) and then converted to ( left, right ).
c     the fixed points on the standard interval are the two end points.      
c
            if(type.eq.'standard') then
               temp(1)=-1.d0
               temp(2)=1.d0
               if(nfixed.eq.1) then
                  temp(1)=1.d0
               endif
               call gaussq('legendre',npts,alpha,beta,nfixed,temp,
     1                                z(scr),z(x),z(wts))
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
            call wtfun(z(x),z(wts),z(xwts),wtfn,npts)
            call cpoly(z(p),z(x),z(xwts),z(a),z(b),left,right,z(scr),
     1                 nmax,npts,pleft(l),pright(l))
            if (prncof) then
                title='a coefficients'
                call prntrm(title,z(a+1),nmax,1,nmax,1,iout)
                title='b coefficients'
                call prntrm(title,z(b+1),nmax-1,1,nmax-1,1,iout)
            endif
            call iosys('does "polynomial a coefficients for l='//
     1                  itoc(lval(l))//'" exist on lamdat',0,0,0,yn)
            if(yn.eq.'yes') then
               call iosys('destroy "polynomial a coefficients for l='//
     1                     itoc(lval(l))//'"',0,0,0,' ')
               call iosys('destroy "polynomial b coefficients for l='//
     1                     itoc(lval(l))//'"',0,0,0,' ')
            endif
            call iosys('write real "polynomial a coefficients for l='//
     1                  itoc(lval(l))//'" to lamdat',nmax+1,z(a),0,' ')
            call iosys('write real "polynomial b coefficients for l='//
     1                  itoc(lval(l))//'" to lamdat',nmax+1,z(b),0,' ')
c
c           calculate the polynomials and their first and 
c                        second derivatives.
c
            call gpoly(z(p),z(dp),z(ddp),z(x),z(a),z(b),left,right,
     1                 pleft(l),pright(l),nmax,npts,.false.)
            call iosys('does "orthogonal polynomials '//
     1                 'for l='//itoc(lval(l))//'" exist on lamdat',0,
     2                  0,0,yn)
            if(yn.eq.'yes') then
               call iosys('destroy "orthogonal polynomials '//
     1                 'for l='//itoc(lval(l))//'"',0,0,0,' ') 
               call iosys('destroy "first derivative of '//
     1                    'orthogonal polynomials for l='
     2                    //itoc(lval(l))//'"',0,0,0,' ')
               call iosys('destroy "second derivative of '//
     1                    'orthogonal polynomials for l='
     2                    //itoc(lval(l))//'"',0,0,0,' ')
            endif    
            call iosys('write real "orthogonal polynomials '//
     1                 'for l='//itoc(lval(l))//'" to lamdat',
     2                  npts*nmax,z(p),0,' ')
            call iosys('write real "first derivative of '//
     1                 'orthogonal polynomials for l='
     2                 //itoc(lval(l))//'" to lamdat',
     3                   npts*nmax,z(dp),0,' ')
            call iosys('write real "second derivative of '//
     1                 'orthogonal polynomials for l='
     2                 //itoc(lval(l))//'" to lamdat',
     3                   npts*nmax,z(ddp),0,' ')
            if (prnply) then
                title='polynomials'
                call prntrm(title,z(p),npts,nmax,npts,nmax,iout)
                title='first derivative of polynomials'
                call prntrm(title,z(dp),npts,nmax,npts,nmax,iout)
                title='second derivative of polynomials'
                call prntrm(title,z(ddp),npts,nmax,npts,nmax,iout)
            endif
            call toorth(z(p),z(dp),z(ddp),z(x),z(xwts),z(dxwts),
     1                  z(ddxwts),wtfn,nmax,npts,.false.)
            if(wtfn.ne.'one') then
               call iosys('does "orthogonal functions '//
     1                    'for l='//itoc(lval(l))//'" exist on '//
     2                    'lamdat',0,0,0,yn)
               if(yn.eq.'yes') then
                  call iosys('destroy "orthogonal functions '//
     1                 'for l='//itoc(lval(l))//'"',0,0,0,' ') 
                  call iosys('destroy "first derivative of '//
     1                       'orthogonal functions for l='
     2                       //itoc(lval(l))//'"',0,0,0,' ')
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
                title='dvr polynomials'
                call prntrm(title,z(pn),npts,nmax,npts,nmax,iout)
                title='first derivative of dvr polynomials'
                call prntrm(title,z(dpn),npts,nmax,npts,nmax,iout)
                title='second derivative of dvr polynomials'
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
 1    format(/,5x,'type polynomial                   = ',a16,/,5x,
     1            'order of polynomials              = ',i3,/,5x,
     2            'left boundary                     = ',e15.8,/,5x,
     3            'right boundary                    = ',e15.8,/,5x,
     4            'derivative at left boundary       = ',e15.8,/,5x,
     5            'derivative at right boundary      = ',e15.8,/,5x,
     6            'number of angular momenta         = ',i2,/5x,
     7            'potential type                    = ',a16,/5x,
     8            'weight function                   = ',a16)
 2    format(/,5x,'number fixed end points = ',i1,/,5x,     
     1            'end point               = ',e15.8)
 3    format(/,5x,'number fixed end points = ',i1,/,5x,     
     1            'left end point          = ',e15.8,/,5x,
     2            'right end point         = ',e15.8)
 4    format(/,5x,'for polynomial type  = ',a16,/,5x,
     1            'alpha = ',e15.8,/,5x,
     2            'beta = ',e15.8)
 5    format(/,5x,'for polynomial type  = ',a16,/,5x,
     1            'alpha = ',e15.8)
 6    format(/,15x,'using co-ordinate representation')
 7    format(/,15x,'beginning calculation for angular momentum = ',
     1                                                    i3,/,5x,
     2             'order of leading left polynomial  = ',i2,/5x,
     3             'order of leading right polynomial = ',i2,/5x,)
 8    format(/,15,'matrix generated')     
 9    format(/,15x,'matrix already calculated-only diagonalization')
      end

