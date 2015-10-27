*deck clin.f
      program clin
c
      implicit integer (a-z)
      character*4096 ops
      character*8 prtflg
      character*80 chrkey, cpass, title, pottyp, precon
      character*1600 card
      character*128 fillam
      logical posinp, logkey, prnh, fix, direct 
      common z(1)
      dimension ia(1)
      dimension endpts(2), temp(2), work(2)
      equivalence (ia(1),z(1))
      common/io/inp, iout      
      common/memory/ioff
      real*8 z, fpkey, thresh, cnverg, energy, rmax, endpts
      real*8 temp, norm0
c
      call drum
      write(iout,*)
      write(iout,*) '    test linslv          '
      call iosys ('read character options from rwf',-1,0,0,ops)
      call iosys ('read character "print flag" from rwf',-1,0,0,prtflg)
      call iosys ('read character "linear algebraic filename" from rwf',
     1             -1,0,0,fillam)
      call iosys ('open lamdat as new',0,0,0,fillam)
      if ( posinp('$linslv',cpass) ) then
           call cardin(card)
      endif 
      nmax=intkey(card,'order-of-polynomials',10,' ')
      npt=intkey(card,'number-of-points',nmax,' ')
      dim=max(nmax,npt)
      fix=logkey(card,'fix-end-points',.false.,' ')
      rmax=fpkey(card,'box-size',1.d0,' ')
      direct=logkey(card,'iterative-linear-solve',.false.,' ')
      pottyp=chrkey(card,'potential-type','well',' ')
      precon=chrkey(card,'precondition','diagonal',' ') 
      if (fix) then
          nfix=intkey(card,'number-of-fixed-points',2,' ')
          endpts(1)=0.d0
          endpts(2)=rmax
          call fparr(card,'end-points',endpts,nfix,' ')
      endif
      pleft=intkey(card,'order-of-leading-left-'//
     1                  'polynomials',1,' ')
      pright=intkey(card,'order-of-leading-right-'//
     1                   'polynomials',0,' ')
      nrhs=intkey(card,'number-of-right-hand-sides',1,' ')
      lenbuf=intkey(ops,'linslv=buffer',min(1000000,nmax*nmax),' ')
      thresh=fpkey(ops,'linslv=tolerance',1.0d-10,' ')
      cnverg=fpkey(ops,'linslv=convergence',1.d-08,' ')
      nattim=intkey(ops,'linslv=number-of-solutions-at-a-time',
     1              nrhs,' ')
      maxvec=intkey(ops,'linslv=maximum-number-of-vectors',
     1                 nmax,' ')
      niter=intkey(ops,'linslv=maximum-number-of-iterations',
     1             nmax,' ')
      prgues=logkey(ops,'print=linslv=guess',.false.,' ')
      prvec=logkey(ops,'print=linslv=final-solutions',.false.,' ')
      prres=logkey(ops,'print=linslv=residuals',.false.,' ')
      prnh=logkey(ops,'print=linslv=hamiltonian',.false.,' ')
      write(iout,1) nmax, nrhs, nattim, lenbuf, thresh, cnverg, maxvec,
     1              niter, precon
      dim=max(npt,nmax)
      ioff=1
      do 10 i=1,2
         q=ioff
         wt=q+npt
         a=wt+npt
         b=a+nmax+1
         p=b+nmax+1
         dp=p+nmax*npt
         ddp=dp+nmax*npt
         pn=ddp+nmax*npt
         dpn=pn+nmax*npt
         ddpn=dpn+nmax*npt
         eigc=ddpn+nmax*npt
         ham=eigc+nmax
         v=ham+2*nmax*nmax
         ur=v+2*nmax*nmax
         ul=ur+2*nmax*nmax
         eig=ul+2*nmax*nmax
         diag=eig+2*nmax
         dinv=diag+2*nmax
         rhs=dinv+2*nmax
         vec=rhs+2*nmax*nrhs
         hvec=vec+2*nmax*maxvec
         smat=hvec+2*nmax*maxvec
         svec=smat+4*maxvec*maxvec
         ipvt=wpadti(svec+4*maxvec*nrhs)
         work(1)=iadtwp(ipvt+nmax)
         work(2)=work(1)+2*dim*dim
         words=wpadti(work(2)+2*dim*dim)
         if (i.eq.1) then
             call iosys ('read integer maxsiz from rwf',1,
     1                    canget,0,' ')
             if (words.gt.canget) then
                 call lnkerr('not enough memory. will quit')
             endif
             call iosys ('write integer maxsiz to rwf',1,
     1                    words,0,' ')
         else
             call getscm(words,z,ngot,'linslv',0)
         endif
   10 continue
c
c     calculate the primitive polynomial and then the dvr polynomial basis.
c
      temp(1)=-1.d0
      temp(2)=1.d0
      if(nfix.eq.1) then
         temp(1)=1.d0
      endif
      call gaussq('legendre',npt,0.d0,0.d0,nfix,temp,z(work(1)),
     1             z(q),z(wt))
      call convt(z(q),z(wt),endpts(1),endpts(2),norm0,npt,.false.)
      call cpoly(z(p),z(q),z(wt),z(a),z(b),endpts(1),endpts(2),
     1           z(work(1)),nmax,npt,pleft,pright)
      call gpoly(z(p),z(dp),z(ddp),z(q),z(a),z(b),endpts(1),endpts(2),
     1           pleft,pright,nmax,npt,.false.)
c      title='polynomials'
c      call prntrm(title,z(p),npt,nmax,npt,nmax,iout)
c      title='first derivative of polynomials'
c      call prntrm(title,z(dp),npt,nmax,npt,nmax,iout)
c      title='second derivative of polynomials'
c      call prntrm(title,z(ddp),npt,nmax,npt,nmax,iout)
      call diagx(z(p),z(dp),z(ddp),z(a+1),z(b+1),z(pn),z(dpn),z(ddpn),
     1           z(eigc),z(work(1)),z(work(2)),nmax,npt,.false.)
c      title='dvr polynomials'
c      call prntrm(title,z(pn),npt,nmax,npt,nmax,iout)
c      title='first derivative of dvr polynomials'
c      call prntrm(title,z(dpn),npt,nmax,npt,nmax,iout)
c      title='second derivative of dvr polynomials'
c      call prntrm(title,z(ddpn),npt,nmax,npt,nmax,iout)
c      call chk(z(pn),z(wt),z(work(1)),z(work(2)),nmax,npt)
c     
c     form the zeroth order hamiltonian in the dvr basis.
c
      call ham0(z(pn),z(dpn),z(ddpn),z(wt),z(ham),nmax,npt)
c
c     form the ( here diagonal) potential energy matrix.
c
      call vmat(z(v),z(eigc),nmax,pottyp)
      if(precon.eq.'diagonalize-h0') then
         call cdiag(z(ham),z(ur),z(ul),z(eig),z(work(1)),nmax)
      endif
      energy=fpkey(card,'energy',0.d0,' ')
c
c     calculate the right hand side of the linear equations.
c     here it is simply taken as the driving term to the driven
c     equations formulation in my paper.
c
      call mkpsi0(z(rhs),energy,z(pn),z(q),z(wt),z(work(1)),pottyp,
     1            nmax,npt,nrhs)
      if(direct) then
         call dlin(z(ham),z(v),z(ur),z(ul),z(eig),z(vec),z(hvec),
     1             energy,z(rhs),z(diag),z(dinv),z(smat),z(svec),
     2             ia(ipvt),z(work(1)),z(work(2)),thresh,cnverg,
     3             ops,nmax,nrhs,nattim,maxvec,niter,precon)
      else
          call drslv(z(ham),z(v),z(rhs),energy,ia(ipvt),nmax,nrhs)
      endif                     
      call iosys ('write integer maxsiz to rwf',1,canget,0,' ')
      call chainx(0)               
      stop
 1    format(/,5x,'number of equations                = ',i3,/,5x,
     1            'number of right hand sides         = ',i3,/5x,
     2            'number solved at a time            = ',i3,/,5x,
     3            'buffer length                      = ',i8,/,5x,
     4            'overlap tolerance                  = ',e15.8,/,5x,
     5            'convergence criterion              = ',e15.8,/,5x,
     6            'maximum size of projected subspace = ',i6,/,5x,
     7            'maximum number of iterations       = ',i6,/,5x,
     8            'preconditioner                     = ',a16)
      end




