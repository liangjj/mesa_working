*deck dvr2dvr.f 
c***begin prologue     dvr2dvr
c***date written       960906   (yymmdd)
c***revision date               (yymmdd)
c***keywords           orthogonal polynomial, roots
c***author             schneider, b. i.(nsf)
c***source             dvr2dvr
c***purpose            relate two dvr eigensets to one another
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       dvr2dvr
      program dvr2dvr
c
      implicit integer (a-z)
      parameter ( maxgrd=50 )
      character*4096 ops
      character*8 cpass
      character*2 itoc
      character*80 title, chrkey
      character*1600 card
      character*16 type, wtfn, sets, pottyp, smooth, typges
      character*128 fillam
      logical posinp, logkey, prncof, prnply, prnwpt, prnh, prnhb, check
      logical prgues, prvec, prres, fixed, direct 
      common z(1)
      dimension ia(1), endpts(2), der(2), temp(2), work(5)
      dimension nfin(maxgrd), nfout(maxgrd)
      equivalence (ia(1),z(1))
      common/io/inp, iout      
      common/memory/ioff
      real*8 z, left, right, fpkey, endpts, alpha, beta, der
      real*8 norm0, temp, thresh, cnverg, energy
c
      call drum
      write(iout,*)
      write(iout,*) '       dvr basis generation'
      call iosys ('read character options from rwf',-1,0,0,ops)
      prncof=logkey(ops,'print=m5001=polynomial-coefficients',
     1              .false.,' ')
      prnply=logkey(ops,'print=m5001=polynomials',.false.,' ')
      prnwpt=logkey(ops,'print=m5001=points/weights',.false.,' ')
      prnh=logkey(ops,'print=m5001=hamiltonian',.false.,' ')
      check=logkey(ops,'check-orthogonality',.false.,' ')
      pottyp=chrkey(ops,'potential','none',' ')
      direct=logkey(ops,'iterative-linear-solve',.false.,' ')
      call iosys ('read character "linear algebraic filename" from rwf',
     1                 -1,0,0,fillam)
      call iosys ('open lamdat as unknown',0,0,0,fillam)
      if ( posinp('$dvr2dvr',cpass) ) then
           call cardin(card)
      endif
      energy=fpkey(card,'energy',0.d0,' ')
      type=chrkey(card,'type-polynomials','standard',' ')
      left=fpkey(card,'left-boundary',-1.d0,' ')
      right=fpkey(card,'right-boundary',1.d0,' ')
      der(1)=fpkey(card,'left-derivative',0.d0,' ')
      der(2)=fpkey(card,'right-derivative',0.d0,' ')
      nmin=intkey(card,'minimum-order-of-polynomials',3,' ')
      ngrid=intkey(card,'number-of-grids',2,' ')
      nfin(1)=nmin
      nfout(1)=nfin(1)
      do 10 ng=2,ngrid
         nfin(ng)=2*nfin(ng-1)
         nfout(ng)=nfin(ng)
 10   continue   
      nmax=nfin(ngrid)
      smooth=chrkey(ops,'linslv=smoother','diagonal',' ')
      typges=chrkey(ops,'linslv=guess','solutions',' ')
      lenbuf=intkey(ops,'linslv=buffer',min(1000000,nmax*nmax),' ')
      thresh=fpkey(ops,'linslv=tolerance',1.0d-10,' ')
      cnverg=fpkey(ops,'linslv=convergence',1.d-08,' ')
      nrhs=intkey(ops,'linslv=number-of-right-hand-sides',1,' ')   
      ntrials=intkey(ops,'linslv=number-of-trial-vectors',nrhs,' ')
      niter=intkey(ops,'linslv=maximum-number-of-iterations',
     1             nmax,' ')
      prgues=logkey(ops,'print=linslv=guess',.false.,' ')
      prvec=logkey(ops,'print=linslv=final-solutions',.false.,' ')
      prres=logkey(ops,'print=linslv=residuals',.false.,' ')
      prnhb=logkey(ops,'print=linslv=hamiltonian',.false.,' ')
      maxpt=nmax
      minpt=nmin
      add=0  
      fixed=logkey(card,'fix-end-points',.false.,' ')
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
      write(iout,1) type, nmax, ngrid, left, right, der(1), 
     1              der(2), wtfn
      if (fixed) then
          if(nfixed.eq.1) then
             write(iout,2) nfixed, endpts(1)
             add=1
          endif
          if(nfixed.eq.2) then
             write(iout,3) nfixed, endpts(1), endpts(2)
             add=2
          endif 
      endif
      if (type.eq.'jacobi') then
          write(iout,4) type, alpha, beta
      endif
      if (type.eq.'laguerre') then
          write(iout,5) type, alpha
      endif
      pleft=intkey(card,'order-of-leading-left-polynomial',1,' ')
      pright=intkey(card,'order-of-leading-right-polynomial',0,' ')
      maxpt=maxpt+add
      dim=maxpt
      ioff=1
      do 20 i=1,2
         x=ioff
         wts=x+dim
         a=wts+dim
         b=a+nmax+1
         p=b+nmax+1
         dp=p+nmax*dim
         ddp=dp+nmax*dim
         eig=ddp+nmax*dim
         pn=eig+nmax
         dpn=pn+nmax*dim
         ddpn=dpn+nmax*dim
         hmat=ddpn+nmax*dim
         tmat=hmat+nmax*nmax
         work(1)=tmat+nmax*nmax
         work(2)=work(1)+dim*dim
         work(3)=work(2)+dim*dim
         work(4)=work(3)+2*nmax*dim
         work(5)=work(4)+4*nmax*nmax
         words=work(5)
         if(direct) then
            hbuf=work(5)
            ibuf=wpadti(hbuf+lenbuf)
            diag=iadtwp(ibuf+2*lenbuf)
            dinv=diag+nmax
            rhs=dinv+nmax
            trials=rhs+nmax*nrhs
            pvec=trials+nmax*ntrials
            hpvec=pvec+nmax*niter
            vec=hpvec+nmax*niter
            tvec=vec+nmax*max(nrhs,ntrials)
            b=tvec+nmax*max(nrhs,ntrials)
            btmp=b+niter*niter
            sol=btmp+niter*niter
            soltmp=niter*nrhs
            ipvt=wpadti(soltmp+niter*nrhs)
            ind=ipvt+nmax
            words=iadtwp(ind+nmax)
         endif
         words=wpadti(words)
         if (i.eq.1) then
             call iosys ('read integer maxsiz from rwf',1,
     1                    canget,0,' ')
             if (words.gt.canget) then
                 call lnkerr('not enough memory. will quit')
             endif
             call iosys ('write integer maxsiz to rwf',1,
     1                    words,0,' ')
         else
             call getscm(words,z,ngot,'dvr2dvr',0)
         endif
   20 continue
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
         call gaussq('legendre',maxpt,alpha,beta,nfixed,temp,
     1                z(work(1)),z(x),z(wts))
c     
c         convert to points and weights on (left,right) if needed
c
         call convt(z(x),z(wts),left,right,norm0,maxpt,
     1              type,prnwpt)
      else
         call gaussq(type,maxpt,alpha,beta,nfixed,endpts,z(work(1)),
     1               z(x),z(wts))
      endif
c
c           calculate the recursion coefficients
c
      call cpoly(z(p),z(x),z(wts),z(a),z(b),left,right,z(work(1)),
     1           nmax,maxpt,pleft,pright)
      if (prncof) then
          title='a coefficients'
          call prntrm(title,z(a+1),nmax,1,nmax,1,iout)
          title='b coefficients'
          call prntrm(title,z(b+1),nmax-1,1,nmax-1,1,iout)
      endif
c
c           calculate the polynomials and their first and 
c                        second derivatives.
c
      call gpoly(z(p),z(dp),z(ddp),z(x),z(a),z(b),left,right,
     1           pleft,pright,nmax,maxpt,.false.)
      if (prnply) then
          title='polynomials'
          call prntrm(title,z(p),maxpt,nmax,maxpt,nmax,iout)
          title='first derivative of polynomials'
          call prntrm(title,z(dp),maxpt,nmax,maxpt,nmax,iout)
          title='second derivative of polynomials'
          call prntrm(title,z(ddp),maxpt,nmax,maxpt,nmax,iout)
      endif
      if(check) then
         call chk(z(p),z(wts),z(work(1)),z(work(2)),nmax,maxpt)
      endif   
      call iosys('open ham as scratch',0,0,0,' ')
      do 30 ng=1,ngrid
         write(iout,6) ng, nfin(ng)
         call diagx(z(p),z(dp),z(ddp),z(a+1),z(b+1),z(pn),z(dpn),
     1              z(ddpn),z(eig),z(work(1)),z(tmat),nfin(ng),
     2              maxpt,ng,prnh)
         if(check) then
            call chk(z(pn),z(wts),z(work(1)),z(work(2)),nfin(ng),maxpt)
         endif   
         call ham(z(pn),z(dpn),z(ddpn),z(wts),z(eig),z(hmat),
     1            nfin(ng),maxpt,pottyp,ng,prnh)
 30   continue  
      do 40 ng=1,ngrid
         write(iout,7) ng
         if(direct) then
            call setind(ia(ind),nfin(ng),nfin(ng),1)
            call iham(z(hmat),z(hbuf),ia(ibuf),ia(ind),z(diag),lenbuf,
     1                nfin(ng),prnhb,ntot,incore,ng)
            call mkrhs(z(rhs),nfin(ng),nrhs)
            if(ng.eq.1) then
               nvecs=nfin(ng)
            else
               nvecs=nfin(ng-1)
            endif
            call guesvc(z(vec),z(tmat),z(rhs),z(trials),z(work(2)),
     1                  nvecs,nrhs,ntrials,nfin(ng),prgues,
     2                  typges,ng)
            call drvlst(z(hbuf),ia(ibuf),z(diag),z(dinv),z(rhs),
     1                  z(trials),z(pvec),z(hpvec),z(vec),z(tvec),
     2                  z(b),z(btmp),z(sol),z(soltmp),cnverg,thresh,
     3                  energy,ia(ipvt),nfin(ng),nrhs,ntrials,niter,
     4                  nfinal,lenbuf,ntot,incore,ops,ng)
            call iosys('write real "solutions for grid '//
     1                  itoc(ng)//'" to ham',nfin(ng)*nfinal,z(pvec),
     2                                       0,' ')
            call tr2ply(z(pvec),z(hpvec),z(tmat),nfin(ng),nvecs,ng)
         else
            call rdiag(z(hmat),z(work(1)),z(work(2)),nfin(ng))
         endif
 40   continue 
      call iosys ('write integer maxsiz to rwf',1,canget,0,' ')
      call chainx(0)               
      stop
 1    format(/,5x,'type polynomial                   = ',a16,/,5x,
     1            'maximum order of polynomials      = ',i3,/,5x,
     2            'number of grids                   = ',i3,/,5x,
     3            'left boundary                     = ',e15.8,/,5x,
     4            'right boundary                    = ',e15.8,/,5x,
     5            'derivative at left boundary       = ',e15.8,/,5x,
     6            'derivative at right boundary      = ',e15.8,/,5x,
     7            'weight function                   = ',a16)
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
 6    format(/,5x,'calculation of DVR basis and Hamiltonian for grid'
     1                                            '  = ',i3,/,5x,
     2            'number of polynomials             = ',i3)
 7    format(/,5x,'solving for solutions for grid = ',i3)
      end

