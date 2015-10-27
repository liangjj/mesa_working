*deck rdiis.f 
c***begin prologue     rdiis
c***date written       960718   (yymmdd)
c***revision date               (yymmdd)
c***keywords           
c***                   
c***author             schneider, b. i.(nsf)
c***source             
c***purpose            real version of diis algorithm
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       rdiis
      program rdiis
c
      implicit integer (a-z)
      parameter (nume=100)
      character*4096 ops
      character*8 prtflg
      character*80 cpass, title, chrkey, precon
      character*1600 card
      character*24 pottyp
      character*128 fillam
      real*8 z, fpkey, left, right, endpts, der, norm0, temp, cnverg
      real*8 rsize, todiis, energy
      logical posinp, logkey, fix, doall, drctv, prnt
      logical switch, type, tsteig, dirslv
      common z(1)
      dimension ia(1), endpts(2), der(2), temp(2), work(2)
      dimension energy(nume)
      equivalence (ia(1),z(1))
      common/io/inp, iout      
      common/memory/ioff
      call drum
      write(iout,*)
      write(iout,*) '    test iteration/diis          '
      call iosys ('read character options from rwf',-1,0,0,ops)
      call iosys ('read character "print flag" from rwf',-1,0,0,prtflg)
      doall=logkey(ops,'all-tests-on',.false.,' ')
      drctv=logkey(ops,'diis=on',.false.,' ')
      type=logkey(ops,'diis=on=jacobi',.false.,' ')
      tsteig=logkey(ops,'test-eigenvalues',.false.,' ')
      precon=chrkey(ops,'preconditioner','rhs',' ')
      dirslv=logkey(ops,'direct-solve',.false.,' ')
      prnt=logkey(ops,'print=diis=on',.false.,' ')
      call iosys ('read character "linear algebraic filename" from rwf',
     1             -1,0,0,fillam)
      call iosys ('open lamdat as unknown',0,0,0,fillam)
      if ( posinp('$diis',cpass) ) then
           call cardin(card)
           pottyp=chrkey(card,'potential-type','none',' ')
           nmax=intkey(card,'order-of-polynomials',30,' ')
           npt=intkey(card,'number-of-points',nmax+2,' ')
           cnverg=fpkey(card,'convergence',1.d-08,' ')
           maxit=intkey(card,'maximum-number-of-iterations',10,' ')
           trunc=intkey(card,'restart-value',5,' ') 
           keep=intkey(card,'iterations-retained',4,' ')
           todiis=fpkey(card,'switch-to-diis',1.d0,' ')
           switch=logkey(card,'no-iteration',.false.,' ')
           rsize=fpkey(card,'region-size',5.0d0,' ')
           nrhs=intkey(card,'number-of-right-hand-sides',1,' ')
           nen=intkey(card,'number-of-energies',1,' ')
           call fparr(card,'energies',energy,nen,' ')
           write(iout,1) nrhs, nen
           fix=logkey(card,'fix-end-points',.false.,' ')
           left=0.d0
           right=rsize
           der(1)=fpkey(card,'left-derivative',0.d0,' ')
           der(2)=fpkey(card,'right-derivative',0.d0,' ')
           if (fix) then
               nfix=intkey(card,'number-of-fixed-points',2,' ')
               endpts(1)=left
               endpts(2)=right
               call fparr(card,'end-points',endpts,nfix,' ')
           endif
           pleft=intkey(card,'order-of-leading-left-polynomials',1,' ')
           pright=intkey(card,'order-of-leading-right-polynomials',
     1                   0,' ')
           write(iout,2)  nmax, npt, left, right, der(1), der(2), 
     1                    pleft, pright, nrhs, maxit, nen
           if (fix) then
               if(nfix.eq.1) then
                  write(iout,3) nfix, endpts(1)
               endif
               if(nfix.eq.2) then
                  write(iout,4) nfix, endpts(1), endpts(2)
               endif 
           endif
      endif
      if(drctv) then
         write(iout,5)
      else
         write(iout,6)
      endif
      maxd=max(nmax,npt)
      ioff=1
      do 10 i=1,2
         q=ioff
         wt=q+maxd
         a=wt+maxd
         b=a+nmax+1
         px=b+nmax+1
         dpx=px+nmax*maxd
         ddpx=dpx+nmax*maxd
         pnx=ddpx+nmax*maxd
         dpnx=pnx+nmax*maxd
         ddpnx=dpnx+nmax*maxd
         pdvr=ddpnx+nmax*maxd
         dpdvr=pdvr+nmax*nmax
         ddpdvr=dpdvr+nmax*nmax
         eigc=ddpdvr+nmax*nmax
         eigcm=eigc+nmax
         tham=eigcm+nmax
         matrix=tham+nmax*nmax
         rhs=matrix+nmax*nmax
         xold=rhs+nmax*nrhs
         err=xold+nmax*(maxit+1)
         xcur=err+nmax*(maxit+1)
         errcur=xcur+nmax
         bb=errcur+nmax
         bbtmp=bb+(maxit+1)*(maxit+1)
         sol=bbtmp+(maxit+1)*(maxit+1)
         ipvt=wpadti(sol+maxit+1)
         work(1)=iadtwp(ipvt+max(maxit+1,nmax))
         work(2)=work(1)+max(3*maxd,2*maxd*maxd)
         words=wpadti(work(2)+2*maxd*maxd)                       
         if (i.eq.1) then
             call iosys ('read integer maxsiz from rwf',1,
     1                    canget,0,' ')
             if (words.gt.canget) then
                 call lnkerr('not enough memory. will quit')
             endif
             call iosys ('write integer maxsiz to rwf',1,
     1                    words,0,' ')
         else
             call getscm(words,z,ngot,'diis',0)
         endif
   10 continue
c
c     in order to proceed it is necessary to be able to calculate integrals
c     between the orthogonal polynomials to be generated on the desired
c     integration interval.  the calculation of the points and weights is
c     done on the interval ( -1.,1. ) and then converted to ( left, right ).
c     the fixed points on the standard interval are the two end points.      
c
      temp(1)=-1.d0
      temp(2)=1.d0
      if(nfix.eq.1) then
         temp(1)=1.d0
      endif
      call gaussq('legendre',npt,0.d0,0.d0,nfix,temp,z(work(1)),
     1             z(q),z(wt))
c     
c           convert to points and weights on (left,right) if needed
c
      call convt(z(q),z(wt),endpts(1),endpts(2),norm0,npt,doall)
c
c           calculate the recursion coefficients
c
      call cpoly(z(px),z(q),z(wt),z(a),z(b),endpts(1),endpts(2),
     1           z(work(1)),nmax,npt,pleft,pright)
      if(doall) then
         title='a coefficients'
         call prntrm(title,z(a+1),nmax,1,nmax,1,iout)
         title='b coefficients'
         call prntrm(title,z(b+1),nmax-1,1,nmax-1,1,iout)
      endif   
c
c           calculate the polynomials and their first and 
c                        second derivatives.
c
      call gpoly(z(px),z(dpx),z(ddpx),z(q),z(a),z(b),endpts(1),
     1           endpts(2),pleft,pright,nmax,npt,doall)
      call diagx(z(px),z(dpx),z(ddpx),z(a+1),z(b+1),z(pnx),
     1           z(dpnx),z(ddpnx),z(eigc),z(work(1)),z(tham),
     2           nmax,npt,doall)
      if(doall) then
         title='polynomials at original points'
         call prntrm(title,z(px),npt,nmax,npt,nmax,iout)
         title='first derivative of polynomials at original points'
         call prntrm(title,z(dpx),npt,nmax,npt,nmax,iout)
         title='second derivative of polynomials at original points'
         call prntrm(title,z(ddpx),npt,nmax,npt,nmax,iout)
         title='dvr polynomials at original points'
         call prntrm(title,z(pnx),npt,nmax,npt,nmax,iout)
         title='first derivative of dvr polynomials at original points'
         call prntrm(title,z(dpnx),npt,nmax,npt,nmax,iout)
         title='second derivative of dvr polynomials at original points'
         call prntrm(title,z(ddpnx),npt,nmax,npt,nmax,iout)
         call chk(z(px),z(wt),z(work(1)),z(work(2)),nmax,npt)
         call chk(z(pnx),z(wt),z(work(1)),z(work(2)),nmax,npt)
         call xmtrx(z(px),z(pnx),z(q),z(wt),z(work(1)),nmax,npt)
      endif   
      call newply(z(px),z(dpx),z(ddpx),z(a),z(b),z(pdvr),z(dpdvr),
     1            z(ddpdvr),z(eigc),z(tham),endpts(1),endpts(2),
     2            pleft,pright,nmax)
      if(doall) then
         title='polynomials at dvr points'
         call prntrm(title,z(px),nmax,nmax,nmax,nmax,iout)
         title='first derivative of polynomials at dvr points'
         call prntrm(title,z(dpx),nmax,nmax,nmax,nmax,iout)
         title='second derivative of polynomials at dvr points'
         call prntrm(title,z(ddpx),nmax,nmax,nmax,nmax,iout)
         title='dvr polynomials at dvr points'
         call prntrm(title,z(pdvr),nmax,nmax,nmax,nmax,iout)
         title='first derivative of dvr polynomials at dvr points'
         call prntrm(title,z(dpdvr),nmax,nmax,nmax,nmax,iout)
         title='second derivative of dvr polynomials at dvr points'
         call prntrm(title,z(ddpdvr),nmax,nmax,nmax,nmax,iout)
      endif
      call ham0(z(pnx),z(dpnx),z(ddpnx),z(wt),z(eigc),z(tham),
     1          nmax,npt,pottyp,.false.)
      if(tsteig) then
         call rdiag(z(tham),z(eigcm),z(work(1)),nmax)
         title='eigenvalues of hamiltonian'
         call prntrm(title,z(eigcm),nmax,1,nmax,1,iout)
      endif 
      do 20 ien=1,nen
         call mkrhs(z(rhs),energy(ien),z(pnx),z(q),z(wt),z(work(1)),
     1              pottyp,nmax,npt,nrhs)
         call linsys(z(tham),z(rhs),energy(ien),z(matrix),z(xold),
     1               z(err),z(xcur),z(errcur),z(bb),z(bbtmp),z(sol),
     2               ia(ipvt),todiis,switch,cnverg,nmax,nrhs,maxit,
     3               trunc,prnt,drctv,type,dirslv,precon)
 20   continue   
      call iosys ('write integer maxsiz to rwf',1,canget,0,' ')
      call chainx(0)               
      stop
 1    format(/,5x,'extrapolating on linear system',/,5x,
     1            'number of right hand sides = ',i3,/,5x,
     2            'number of energies         = ',i3)      
 2    format(/,5x,'order of radial polynomials       = ',i3,/,5x,
     1            'number of integration points      = ',i3,/,5x,
     2            'left boundary                     = ',e15.8,/,5x,
     3            'right boundary                    = ',e15.8,/,5x,
     4            'derivative at left boundary       = ',e15.8,/,5x,
     5            'derivative at right boundary      = ',e15.8,/,5x,
     6            'order-of-left-polynomials         = ',i2,/,5x,
     7            'order-of-right-polynomials        = ',i2,/,5x,
     8            'number of right hand sides        = ',i2,/,5x,
     9            'maximum number of diis iterations = ',i3,/,5x,
     x            'number of energies                = ',i2)
 3    format(/,5x,'number fixed radial end points = ',i1,/,5x,     
     1            'end point                      = ',e15.8)
 4    format(/,5x,'number fixed radial end points = ',i1,/,5x,     
     1            'left end point                 = ',e15.8,/,5x,
     2            'right end point                = ',e15.8)
 5    format(/,15x,'diis procedure is turned on')
 6    format(/,15x,'diis procedure is turned off')
      end

