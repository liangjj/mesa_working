*deck lstsq.f
c***begin prologue     lstsq
c***date written       960906   (yymmdd)
c***revision date               (yymmdd)
c***keywords           orthogonal polynomial, roots
c***author             schneider, b. i.(nsf)
c***source             lstsq
c***purpose            
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       lstsq
      program lstsq
c
      implicit integer (a-z)
      character*4096 ops
      character*8 cpass
      character*80 title, chrkey
      character*1600 card
      character*16 type, wtfn, pottyp
      character*128 fillam
      logical posinp, logkey, prncof, prnply, prnwpt, prnh, check
      logical fixed 
      common z(1)
      dimension ia(1), endpts(2), der(2), temp(2), work(5)
      equivalence (ia(1),z(1))
      common/io/inp, iout      
      common/memory/ioff
      real*8 z, left, right, fpkey, endpts, alpha, beta, der
      real*8 norm0, temp, cnverg, energy, thresh
c
      call drum
      write(iout,*)
      write(iout,*) '       least squares solution of linear equation'
      call iosys ('read character options from rwf',-1,0,0,ops)
      prncof=logkey(ops,'print=m5010=polynomial-coefficients',
     1              .false.,' ')
      prnply=logkey(ops,'print=m5010=polynomials',.false.,' ')
      prnwpt=logkey(ops,'print=m5010=points/weights',.false.,' ')
      prnh=logkey(ops,'print=m5010=hamiltonian',.false.,' ')
      check=logkey(ops,'check-orthogonality',.false.,' ')
      pottyp=chrkey(ops,'potential','none',' ')
      call iosys ('read character "linear algebraic filename" from rwf',
     1                 -1,0,0,fillam)
      call iosys ('open lamdat as unknown',0,0,0,fillam)
      if ( posinp('$lstsq',cpass) ) then
           call cardin(card)
      endif
      nmax=intkey(card,'order-of-polynomials',3,' ')
      nrhs=intkey(card,'number-of-right-hand-sides',1,' ')
      ntrials=intkey(card,'number-of-trial-vectors',nrhs,' ')
      cnverg=fpkey(card,'convergence',1.d-08,' ')
      thresh=fpkey(card,'overlap-tolerance',1.d-08,' ')
      niter=intkey(card,'maximum-number-of-iterations',nmax,' ')
      energy=fpkey(card,'energy',0.d0,' ')
      type=chrkey(card,'type-polynomials','standard',' ')
      left=fpkey(card,'left-boundary',-1.d0,' ')
      right=fpkey(card,'right-boundary',1.d0,' ')
      der(1)=fpkey(card,'left-derivative',0.d0,' ')
      der(2)=fpkey(card,'right-derivative',0.d0,' ')
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
      write(iout,1) type, nmax, left, right, der(1), 
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
      npt=nmax+add
      dim=npt
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
         diag=hmat+nmax*nmax 
         eham=diag+nmax
         tmat=eham+nmax*nmax
         work(1)=tmat+nmax*nmax
         work(2)=work(1)+dim*dim
         rhs=work(2)+dim*dim
         trials=rhs+nmax*nrhs
         pvec=trials+nmax*ntrials
         hpvec=pvec+niter*nmax
         vec=hpvec+niter*nmax
         tvec=vec+nmax*nrhs
         b=tvec+nmax*nrhs
         btmp=b+niter*niter
         sol=btmp+niter*niter
         soltmp=sol+niter*nrhs
         ipvt=wpadti(soltmp+niter*nrhs)
         words=ipvt+nmax
         if (i.eq.1) then
             call iosys ('read integer maxsiz from rwf',1,
     1                    canget,0,' ')
             if (words.gt.canget) then
                 call lnkerr('not enough memory. will quit')
             endif
             call iosys ('write integer maxsiz to rwf',1,
     1                    words,0,' ')
         else
             call getscm(words,z,ngot,'lstsq',0)
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
         call gaussq('legendre',npt,alpha,beta,nfixed,temp,
     1                z(work(1)),z(x),z(wts))
c     
c         convert to points and weights on (left,right) if needed
c
         call convt(z(x),z(wts),left,right,norm0,npt,
     1              type,prnwpt)
      else
         call gaussq(type,npt,alpha,beta,nfixed,endpts,z(work(1)),
     1               z(x),z(wts))
      endif
c
c           calculate the recursion coefficients
c
      call cpoly(z(p),z(x),z(wts),z(a),z(b),left,right,z(work(1)),
     1           nmax,npt,pleft,pright)
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
     1           pleft,pright,nmax,npt,.false.)
      if (prnply) then
          title='polynomials'
          call prntrm(title,z(p),npt,nmax,npt,nmax,iout)
          title='first derivative of polynomials'
          call prntrm(title,z(dp),npt,nmax,npt,nmax,iout)
          title='second derivative of polynomials'
          call prntrm(title,z(ddp),npt,nmax,npt,nmax,iout)
      endif
      if(check) then
         call chk(z(p),z(wts),z(work(1)),z(work(2)),nmax,npt)
      endif   
      call diagx(z(p),z(dp),z(ddp),z(a+1),z(b+1),z(pn),z(dpn),
     1           z(ddpn),z(eig),z(work(1)),z(tmat),nmax,
     2           npt,prnh)
      if(check) then
         call chk(z(pn),z(wts),z(work(1)),z(work(2)),nmax,npt)
      endif   
      call ham(z(pn),z(dpn),z(ddpn),z(wts),z(eig),z(hmat),
     1         nmax,npt,pottyp,prnh)
      call mkrhs(z(rhs),nmax,nrhs)
      call copy(z(rhs),z(trials),nmax*nrhs)
      call drvlst(z(hmat),z(diag),z(rhs),z(trials),z(pvec),z(hpvec),
     1            z(vec),z(tvec),z(b),z(btmp),z(sol),z(soltmp),energy,
     2            cnverg,thresh,ia(ipvt),nmax,nrhs,
     3            ntrials,niter,ops)     
      call iosys ('write integer maxsiz to rwf',1,canget,0,' ')
      call chainx(0)               
      stop
 1    format(/,5x,'type polynomial                   = ',a16,/,5x,
     1            'order of polynomials              = ',i3,/,5x,
     2            'left boundary                     = ',e15.8,/,5x,
     3            'right boundary                    = ',e15.8,/,5x,
     4            'derivative at left boundary       = ',e15.8,/,5x,
     5            'derivative at right boundary      = ',e15.8,/,5x,
     6            'weight function                   = ',a16)
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

