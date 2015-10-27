*deck lin.f
      program lin
c
      implicit integer (a-z)
      character*4096 ops
      character*8 prtflg
      character*80 cpass, title
      character*1600 card
      character*128 fillam
      logical dollar, logkey, prnh, fix, nfix, direct, diagn 
      common z(1)
      dimension ia(1)
      dimension work(2), nfix(2), energy(100)
      equivalence (ia(1),z(1))
      common/io/inp, iout      
      common/memory/ioff
      real*8 z, fpkey, thresh, cnverg, lftept, rtept, lftder, rtder
      real*8 energy
      call drum
      write(iout,*)
      write(iout,*) '    test complex linslv          '
      call iosys ('read character options from rwf',-1,0,0,ops)
      call iosys ('read character "print flag" from rwf',-1,0,0,prtflg)
      call iosys ('read character "linear algebraic filename" from rwf',
     1             -1,0,0,fillam)
      call iosys ('open lamdat as new',0,0,0,fillam)
      if ( dollar('$clinslv',card,cpass,inp) ) then
         npt=intkey(card,'number-of-points',nmax,' ')
         nen=intkey(card,'number-of-energies',1,' ')
         diagn=logkey(card,'diagonalize',.false.,' ')
         call fparr(card,'energies',energy,nen,' ')
         nmax=npt
         lftept=fpkey(card,'left-end-point',0.d0,' ')
         rtept=fpkey(card,'right-end-point',0.d0,' ')
         fix=logkey(card,'fix-end-points',.false.,' ')
         lftbc=intkey(card,'left-boundary-condition',0,' ')
         rtbc=intkey(card,'right-boundary-condition',0,' ')
         if(lftbc.eq.0) then
            nmax=nmax-1
         endif
         if(rtbc.eq.0) then
            nmax=nmax-1
         endif                    
         lftder=fpkey(card,'left-derivative',0.d0,' ')
         rtder=fpkey(card,'right-derivative',0.d0,' ')      
         direct=logkey(card,'iterative-linear-solve',.false.,' ')
         nfix(1)=.false.
         nfix(2)=.false.
         if (fix) then
             nfix(1)=logkey(card,'fix-left-end-point',.false.,' ')
             nfix(2)=logkey(card,'fix-right-end-point',.false.,' ')
         endif
         nrhs=intkey(card,'number-of-right-hand-sides',1,' ')
      endif 
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
     1              niter
      dim=max(npt,nmax)
      ioff=1
      do 10 i=1,2
         q=ioff
         wt=q+npt
         eigc=wt+npt
         wtc=eigc+npt
         p=wtc+npt
         dp=p+npt*npt
         ddp=dp+npt*npt
         pn=ddp+npt*npt
         dpn=pn+npt*npt
         ddpn=dpn+npt*npt
         v=ddpn+npt*npt
         ham=v+2*npt
         words=ham+2*nmax*nmax
         if(diagn) then
            eig=words
            vecr=eig+2*nmax
            temp=vecr+2*nmax*nmax
            words=temp+3*nmax
         endif            
         rhs=words
         diag=rhs+2*nmax*nrhs
         dinv=diag+2*nmax
         vec=dinv+2*nmax
         hvec=vec+2*nmax*maxvec
         smat=hvec+2*nmax*maxvec
         svec=smat+4*maxvec*maxvec
         hbuf=svec+4*maxvec*nrhs
         ibuf=wpadti(hbuf+2*lenbuf)
         ipvt=ibuf+2*lenbuf
         ind=ipvt+nmax
         work(1)=iadtwp(ind+2*nmax)
         work(2)=work(1)+dim*dim
         words=wpadti(work(2)+dim*dim)
         if (i.eq.1) then
             call iosys ('read integer maxsiz from rwf',1,
     1                    canget,0,' ')
             if (words.gt.canget) then
                 call lnkerr('not enough memory. will quit')
             endif
             call iosys ('write integer maxsiz to rwf',1,
     1                    words,0,' ')
         else
             call getscm(words,z,ngot,'clinslv',0)
         endif
   10 continue
      call getqpt(z(q),z(wt),lftept,rtept,'legendre',
     1            'before',z(work(1)),nfix,npt,npt,1,.false.)
      call lgrply(z(p),z(dp),z(ddp),z(q),z(q),npt,npt,.false.)
      call prepfn(z(p),z(dp),z(ddp),z(q),z(wt),npt,1.d0,'legendre')
      call copy(z(p),z(pn),npt*npt)
      call copy(z(dp),z(dpn),npt*npt)
      call copy(z(ddp),z(ddpn),npt*npt)
c      title='dvr polynomials'
c      call prntfm(title,z(pn),npt,npt,npt,npt,iout)
c      title='first derivative of dvr polynomials '
c      call prntfm(title,z(dpn),npt,npt,npt,npt,iout)
c      title='second derivative of dvr polynomials '
c      call prntfm(title,z(ddpn),npt,npt,npt,npt,iout)
c      call chk(z(pn),z(wt),z(work(1)),z(work(2)),npt,npt)
c      title='overlap matrix dvr polynomials'
c      call prntfm(title,z(work(2)),npt,npt,npt,npt,iout)      
      if(lftbc.eq.0.or.rtbc.eq.0) then
         call dropfn(z(p),z(dp),z(ddp),z(pn),z(dpn),z(ddpn),
     1               z(q),z(wt),z(eigc),z(wtc),lftbc,rtbc,
     2               npt,nmax)
      endif
c      title='dvr polynomials'
c      call prntfm(title,z(pn),nmax,nmax,nmax,nmax,iout)
c      title='first derivative of dvr polynomials '
c      call prntfm(title,z(dpn),nmax,nmax,nmax,nmax,iout)
c      title='second derivative of dvr polynomials '
c      call prntfm(title,z(ddpn),nmax,nmax,nmax,nmax,iout)
      call vmat(z(eigc),z(v),nmax,.true.,ops)      
      call cham(z(pn),z(dpn),z(ddpn),z(eigc),z(wtc),z(ham),z(eig),
     1          z(vecr),z(v),z(temp),nmax,lftbc,rtbc,diagn,prnh)
      if(.not.diagn) then
         call mkrhs(z(pn),z(eigc),z(wtc),z(rhs),nmax)
         if(direct) then
            call setind(ia(ind),nmax,nmax,1)
            call ihamxt(z(ham),z(hbuf),ia(ibuf),ia(ind),z(diag),lenbuf,
     1                  nmax,prnh,ntot,incore)
         endif
         do 20 i=1,nen
            if(direct) then
               call drvlin(z(hbuf),ia(ibuf),z(diag),z(vec),z(hvec),
     1                     z(rhs),z(dinv),z(smat),z(svec),energy(i),
     2                     ia(ipvt),thresh,cnverg,ops,nmax,nrhs,
     3                     nattim,maxvec,niter,lenbuf,ntot,incore)
            else
                call drslv(z(ham),z(rhs),energy(i),ia(ipvt),nmax,nrhs)
            endif
            call compar(z(eigc),z(wtc),z(pn),z(rhs),energy(i),
     1                  rtbc,nmax)
 20      continue   

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
     7            'maximum number of iterations       = ',i6)
      end




