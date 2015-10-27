*deck threed.f 
c***begin prologue     threed
c***date written       960519   (yymmdd)
c***revision date               (yymmdd)
c***keywords           polynomial, cartesian, atom
c***author             schneider, b. i.(nsf)
c***source             mesa
c***purpose            driver for solving three dimensional schroedinger
c***                   equation in a box.
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       threed
      program threed
c
      implicit integer (a-z)
      parameter( ngmax=10 )
      character*4096 ops
      character*8 cpass 
      character*80 chrkey
      character*128 title
      character*1600 card
      character*16 pottyp 
      character*128 fillam, filham
      character*2 itoc
      logical posinp, logkey, prncof, prnply, prnwpt, check, prncoe
      logical prnke, prnkev, prn2h, prn2v, prn3h, prn3v
      logical prnprj, prgues, useone, prnrmt
      logical itdiag, incore, srfst
      common z(1)
      dimension ia(1), endpts(2,3), der(2), temp(2), nply(ngmax)
      dimension npts(ngmax), q(3), wts(3), a(3), b(3), p(3), dp(3)
      dimension ddp(3), eigx(3), t(3), scr(12)
      equivalence (ia(1),z(1))
      common/io/inp, iout      
      common/memory/ioff
      real*8 z, fpkey, endpts, alpha, beta, der
      real*8 norm0, temp, thresh, cnverg,dum
c
      call drum
      write(iout,*)
      write(iout,*) 'atom in three-dimensional box '//
     1              ' using cartesian coordinates'
      call iosys ('read character options from rwf',-1,0,0,ops)
      prncof=logkey(ops,'print=m6275=polynomial-coefficients',
     1              .false.,' ')
      prnply=logkey(ops,'print=m6275=polynomials',.false.,' ')
      prnwpt=logkey(ops,'print=m6275=points/weights',.false.,' ')
      prncoe=logkey(ops,'print=m6275=coordinate-eigenvectors',
     1              .false.,' ')
      prnke=logkey(ops,'print=m6275=kinetic-energy-matrix',.false.,' ')
      prnkev=logkey(ops,'print=m6275=kinetic-energy-eigenvectors',
     1              .false.,' ')
      prn2h=logkey(ops,'print=m6275=2d-hamiltonian',.false.,' ')
      prn2v=logkey(ops,'print=m6275=2d-eigenvectors',.false.,' ')
      prn3h=logkey(ops,'print=m6275=3d-hamiltonian',.false.,' ')
      prn3v=logkey(ops,'print=m6275=3d-eigenvectors',.false.,' ')
      prnprj=logkey(ops,'print=m6275=surface-projections',.false.,' ')
      prnrmt=logkey(ops,'print=m6275=r-matrix',.false.,' ')
      check=logkey(ops,'check-orthogonality',.false.,' ')
      itdiag=logkey(ops,'iterative-diagonalization',.false.,' ')
      call iosys ('read character "linear algebraic filename" from rwf',
     1                 -1,0,0,fillam)
      call iosys ('open lamdat as unknown',0,0,0,fillam)
      if ( posinp('$3-dim',cpass) ) then
           call cardin(card)
      endif           
      pottyp=chrkey(card,'potential','none',' ')
      srfst=logkey(card,'one-surface-set',.false.,' ')
      endpts(1,1)=fpkey(card,'x-left',-1.d0,' ')
      endpts(2,1)=fpkey(card,'x-right',1.d0,' ')
      endpts(1,2)=fpkey(card,'y-left',-1.d0,' ')
      endpts(2,2)=fpkey(card,'y-right',1.d0,' ')
      endpts(1,3)=fpkey(card,'z-left',-1.d0,' ')
      endpts(2,3)=fpkey(card,'z-right',1.d0,' ')
      der(1)=fpkey(card,'left-derivative',0.d0,' ')
      der(2)=fpkey(card,'right-derivative',0.d0,' ')
      pleft=intkey(card,'order-of-leading-left-polynomial',0,' ')
      pright=intkey(card,'order-of-leading-right-polynomial',0,' ')
      nmin=intkey(card,'minimum-order-of-polynomials',3,' ')
      ngrids=intkey(card,'number-of-grids',1,' ')
      nen=intkey(card,'number-of-energies',0,' ')
      nply(1)=nmin
      npts(1)=nply(1)+2
      lenply=nply(1)
      lenpts=npts(1)
      lenpq=lenply*lenpts
      lenpp=lenply*lenply
      lenqq=lenpts*lenpts
      do 10 ng=2,ngrids
         nply(ng)=nply(ng-1)+nply(ng-1)-1
         npts(ng)=nply(ng)+2
         lenply=lenply+nply(ng)
         lenpts=lenpts+npts(ng)
         lenpp=lenpp+nply(ng)*nply(ng)
         lenpq=lenpq+nply(ng)*npts(ng)
         lenqq=lenqq+npts(ng)*npts(ng)
 10   continue
      nmax1=max(nply(ngrids)*nply(ngrids),npts(ngrids)*npts(ngrids),
     1          nply(ngrids)*npts(ngrids))
      nmax2=nply(ngrids)*nply(ngrids)
      nmax3=nmax2*nply(ngrids)
      mmax2=npts(ngrids)*npts(ngrids)
      mmax3=mmax2*npts(ngrids)
      nmax=max(nmax2,mmax2,nmax3,mmax3)
      tri=nply(ngrids)*(nply(ngrids)+1)/2
      alpha=0.d0
      beta=0.d0
      write(iout,1) nmin, ngrids, endpts(1,1), endpts(2,1),
     1              endpts(1,2), endpts(2,2), endpts(1,3), 
     2              endpts(2,3), der(1), der(2), pleft, 
     3              pright, pottyp
      if(itdiag) then
         nrt2d=intkey(ops,'davidson=number-of-2d-roots',nmax2,' ')
         nrt3d=intkey(ops,'davidson=number-of-3d-roots',nmax3,' ')
         nrt2d=min(nrt2d,nmax2)
         nrt3d=min(nrt3d,nmax3)
         nroots=max(nrt2d,nrt3d)
         prgues=logkey(ops,'print=davidson=guess',.false.,' ')
         thresh=fpkey(ops,'davidson=tolerance',1.0d-06,' ')
         cnverg=fpkey(ops,'davidson=convergence',1.d-08,' ')
         nattim=intkey(ops,'davidson=number-of-roots-at-a-time',
     1                 nroots,' ')
         maxvec=intkey(ops,'davidson=maximum-number-of-vectors',
     1                 10*nroots,' ')
         niter=intkey(ops,'davidson=maximum-number-of-iterations',
     1                10*nroots,' ')
         useone=logkey(ops,'davidson=guess=one',.false.,' ')
         lenbuf=intkey(ops,'davidson=hamiltonian-buffer',
     1                 min(1000000,nmax3*nmax3),' ')
c         maxvec=min(maxvec,nmax)
c         niter=min(niter,nmax)
         call iosys ('read character "hamiltonian filename" from rwf',
     1                -1,0,0,filham)
         call iosys('open ham as new',0,0,0,filham)
         call iosys('write integer "buffer size" to ham',1,
     1               lenbuf,0,' ')           
         write(iout,2) nrt2d, nrt3d, nattim, maxvec, niter, thresh, 
     1                 cnverg, lenbuf
      else
         nrt2d=intkey(ops,'direct=number-of-2d-roots',nmax2,' ')
         nrt3d=intkey(ops,'direct=number-of-3d-roots',nmax3,' ')
         nrt2d=min(nrt2d,nmax2)
         nrt3d=min(nrt3d,nmax3)
         nroots=max(nrt2d,nrt3d)
         write(iout,3) nrt2d, nrt3d
      endif
      nchan=6*nrt2d
      ioff=1
      do 20 i=1,2
         q(1)=ioff
         q(2)=q(1)+lenpts
         q(3)=q(2)+lenpts
         wts(1)=q(3)+lenpts
         wts(2)=wts(1)+lenpts
         wts(3)=wts(2)+lenpts
         a(1)=wts(3)+lenpts
         a(2)=a(1)+lenply+ngrids
         a(3)=a(2)+lenply+ngrids
         b(1)=a(3)+lenply+ngrids
         b(2)=b(1)+lenply+ngrids
         b(3)=b(2)+lenply+ngrids
         p(1)=b(3)+lenply+ngrids
         p(2)=p(1)+lenpq
         p(3)=p(2)+lenpq
         dp(1)=p(3)+lenpq
         dp(2)=dp(1)+lenpq
         dp(3)=dp(2)+lenpq
         ddp(1)=dp(3)+lenpq
         ddp(2)=ddp(1)+lenpq
         ddp(3)=ddp(2)+lenpq
         t(1)=ddp(3)+lenpq
         t(2)=t(1)+max(lenpp,lenqq)
         t(3)=t(2)+max(lenpp,lenqq)
         scr(1)=t(3)+max(lenpp,lenqq)
         locsc1=scr(1)
         scr(2)=scr(1)+nmax1
         scr(3)=scr(2)+nmax1
         scr(4)=scr(3)+nmax1
         pn=scr(4)+nmax1
         dpn=pn+npts(ngrids)*nply(ngrids)
         ddpn=dpn+npts(ngrids)*nply(ngrids)
         eigx(1)=ddpn+npts(ngrids)*nply(ngrids)
         eigx(2)=eigx(1)+lenply
         eigx(3)=eigx(2)+lenply
         wds0=eigx(3)+lenply
c        depending on direct or davidson diagonalization, memory is allocated
c        differently.  in the direct approach there are also differences
c        depending on whether some or all of the eigenpairs are computed.
         if(.not.itdiag) then
            ham=wds0
            vec=ham
            eig=ham+nmax3*nmax3
            indx=wpadti(eig+nmax3)
            ipvt=indx+4*nmax3
            enddag=iadtwp(ipvt+3*nmax3)
            if(nrt2d.ge.nmax2.and.nrt3d.ge.nmax3) then
               scr(5)=enddag
               locsc2=scr(5)
               scr(6)=scr(5)+nmax3
               wds0=scr(6)+nmax3
               vec=wds0
            else
               scr(5)=enddag
               locsc2=scr(5)
               scr(6)=scr(5)+nmax3
               scr(7)=scr(6)+nmax3
               scr(8)=scr(7)+nmax3
               scr(9)=scr(8)+nmax3
               scr(10)=scr(9)+nmax3
               scr(11)=scr(10)+nmax3
               scr(12)=scr(11)+nmax3
               vec=scr(12)+nmax3
               wds0=vec+nmax3*nroots
            endif
         else
            vec=wds0
            eig=vec+nmax3*max(maxvec,nrt3d)
            indx=wpadti(eig+nmax3)
            ipvt=indx+4*nmax3
            enddag=iadtwp(ipvt+3*nmax3)
            ttmp=enddag
            hbuf=ttmp+3*nmax2
            diag=hbuf+lenbuf
            ibuf=wpadti(diag+nmax3)
            hvec=ibuf+2*lenbuf
            dvdmat=hvec+nmax3*max(maxvec,nrt3d)
            dvdvec=dvdmat+maxvec*maxvec
            cvn=dvdvec+maxvec*maxvec
            hmev=cvn+nmax3
            eig0=hmev+nmax3*maxvec
            wds0=eig0+nmax3
         endif
         eig2d=enddag
         eig3d=eig2d+nrt2d
         srf=eig3d+nrt3d
         vec2d=srf+6*nrt2d*nrt3d
         energy=vec2d
         vec3d=vec2d+max(nrt2d*nmax2,nen)
         rmt=vec3d 
         prj=vec3d+max(nrt3d*nmax3,nen*nchan*nchan)
         wds1=prj+nmax2*nrt3d
         words=wpadti(wds0,wds1)
         if (i.eq.1) then
             write(iout,99) iadtwp(words)
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
   20 continue
c           
c                      now for all the grids 
c  
c                     set up the orthogonal polynomials
c                            in each dimension
c                     then calculate the one-dimensional
c                           kinetic energy matrix    
c
c     in order to proceed it is necessary to be able to calculate integrals
c     between the orthogonal polynomials to be generated on the desired
c     integration interval.  the calculation of the points and weights is
c     done on the interval ( -1.,1. ) and then converted to ( left, right ).
c     the fixed points on the standard interval are the two end points.      
c
      temp(1)=-1.d0
      temp(2)=1.d0
      do 30 nq=1,3
         xloc=q(nq)
         wtloc=wts(nq)
         aloc=a(nq)
         bloc=b(nq)
         ploc=p(nq)
         dploc=dp(nq)
         ddploc=ddp(nq)
         eigloc=eigx(nq)
         tloc=t(nq)
         do 40 ng=1,ngrids
            call gaussq('legendre',npts(ng),alpha,beta,2,temp,z(scr(1)),
     1                   z(xloc),z(wtloc))
c     
c        convert to points and weights on (left,right)
c
            call convt(z(xloc),z(wtloc),endpts(1,nq),endpts(2,nq),
     1                 norm0,npts(ng),'standard',prnwpt)
c     
c
c           calculate the recursion coefficients
c
            call cpoly(z(ploc),z(xloc),z(wtloc),z(aloc),z(bloc),
     1                 endpts(1,nq),endpts(2,nq),z(scr(1)),nply(ng),
     2                 npts(ng),pleft,pright)
            if (prncof) then
                title='a coefficients'
                call prntrm(title,z(aloc+1),nply(ng),1,nply(ng),
     1                      1,iout)
                title='b coefficients'
                call prntrm(title,z(bloc+1),nply(ng)-1,1,
     1                      nply(ng)-1,1,iout)
            endif
c           calculate the polynomials and their first and 
c                        second derivatives.
c
            call gpoly(z(ploc),z(dploc),z(ddploc),z(xloc),z(aloc),
     1                 z(bloc),endpts(1,nq),endpts(2,nq),pleft,
     2                 pright,nply(ng),npts(ng),.false.)
            if(check) then
               call chk(z(ploc),z(wtloc),z(scr(1)),z(scr(2)),nply(ng),
     1                  npts(ng))
            endif   
            if (prnply) then
                title='primitive polynomials'
                call prntrm(title,z(ploc),npts(ng),nply(ng),npts(ng),
     1                      nply(ng),iout)
                title='first derivative of primitive polynomials'
                call prntrm(title,z(dploc),npts(ng),nply(ng),npts(ng),
     1                      nply(ng),iout)
                title='second derivative of primitive polynomials'
                call prntrm(title,z(ddploc),npts(ng),nply(ng),npts(ng),
     1                      nply(ng),iout)
            endif
c
c                   diagonalize the coordinate operator            
            call diagx(z(ploc),z(dploc),z(ddploc),z(aloc+1),
     1                 z(bloc+1),z(pn),z(dpn),z(ddpn),
     2                 z(eigloc),z(scr(1)),z(scr(2)),nply(ng),
     3                 npts(ng),prncoe,nq,ng)
            if(check) then
               call chk(z(pn),z(wtloc),z(scr(1)),z(scr(2)),
     1                  nply(ng),npts(ng))
            endif
            call copy(z(pn),z(ploc),nply(ng)*npts(ng))
            call copy(z(dpn),z(dploc),nply(ng)*npts(ng))
            call copy(z(ddpn),z(ddploc),nply(ng)*npts(ng))
            if (prnply) then
                title='dvr polynomials'
                call prntrm(title,z(ploc),npts(ng),nply(ng),npts(ng),
     1                      nply(ng),iout)
                title='first derivative of dvr polynomials'
                call prntrm(title,z(dploc),npts(ng),nply(ng),npts(ng),
     1                      nply(ng),iout)
                title='second derivative of dvr polynomials'
                call prntrm(title,z(ddploc),npts(ng),nply(ng),npts(ng),
     1                      nply(ng),iout)
            endif
            call tmat(z(ploc),z(dploc),z(ddploc),z(tloc),z(wtloc),
     1                der,nply(ng),npts(ng),prnke)
            call iosys('write real "kinetic energy plus bloch '//
     1                 'operator for coordinate '//itoc(nq)//' and '//
     2                 'grid '//itoc(ng)//'" to lamdat',
     3                  nply(ng)*nply(ng),z(tloc),0,' ')
            xloc=xloc+npts(ng)
            wtloc=wtloc+npts(ng) 
            aloc=aloc+nply(ng)
            bloc=bloc+nply(ng)+1
            ploc=ploc+npts(ng)*nply(ng)+1
            dploc=dploc+npts(ng)*nply(ng)
            ddploc=ddploc+npts(ng)*nply(ng)
            tloc=tloc+nply(ng)*nply(ng)
            eigloc=eigloc+nply(ng)
 40      continue   
 30   continue
c
c           the remaining code is broken into two sections depending
c           on whether a direct or davidson diagonalization is performed.
c
      cord1=q(1)
      cord2=q(2)
      cord3=q(3)
      wt1=wts(1)
      wt2=wts(2)
      wt3=wts(3)
      p1=p(1)
      p2=p(2)
      p3=p(3)
      t1=t(1)
      t2=t(2)
      t3=t(3)
      eig1=eigx(1)
      eig2=eigx(2)
      eig3=eigx(3)
      do 50 ng=1,ngrids
         nsq=nply(ng)*nply(ng)
         if(.not.itdiag) then
c      
c        calculate all matrices explicitly and diagonalize        
c
c     solve first for the states in two-dimensions.  These are the surface
c     states and can be thought of as the analogs of channel functions.
c     these are calculated by restricting the potential function to the
c     corresponding plane, setting up the hamiltonian and diagonalizing.
c     in general the potential function is not symmetric in (x,y,z) so
c     six diagonalizations are required.       
            call setind(ia(indx),nply(ng),nsq,2)
            title='2d surface functions in right (q1,q2) plane'
            call set2x(z(eig1),z(eig2),endpts(2,3),z(ham),
     1                 z(eig),z(vec),z(t1),z(t2),z(locsc2),
     2                 ia(indx),ia(ipvt),nply(ng),npts(ng),nsq,nrt2d,
     3                 pottyp,prn2h,prn2v,title)
            call iosys('write real "surface eigenvalues for right '//
     1                 '(q1,q2) surface for grid '//itoc(ng)//
     2                 '" to lamdat',nrt2d,z(eig),0,' ')
            call iosys('write real "surface functions for right '//
     1                 '(q1,q2) surface for grid '//itoc(ng)//
     2                 '" to lamdat',nsq*nrt2d,z(ham),0,' ')
            title='2d surface functions in left (q1,q2) plane'
            call set2x(z(eig1),z(eig2),endpts(1,3),z(ham),
     1                 z(eig),z(vec),z(t1),z(t2),z(locsc2),
     2                 ia(indx),ia(ipvt),nply(ng),npts(ng),nsq,nrt2d,
     3                 pottyp,prn2h,prn2v,title)
            call iosys('write real "surface eigenvalues for left '//
     1                 '(q1,q2) surface for grid '//itoc(ng)//
     2                 '" to lamdat',nrt2d,z(eig),0,' ')
            call iosys('write real "surface functions for left '//
     1                 '(q1,q2) surface for grid '//itoc(ng)//
     2                 '" to lamdat',nsq*nrt2d,z(ham),0,' ')
            title='2d surface functions in right (q1,q3) plane'
            call set2x(z(eig1),z(eig3),endpts(2,2),z(ham),
     1                 z(eig),z(vec),z(t1),z(t2),z(locsc2),
     2                 ia(indx),ia(ipvt),nply(ng),npts(ng),nsq,nrt2d,
     3                 pottyp,prn2h,prn2v,title)
            call iosys('write real "surface eigenvalues for right '//
     1                 '(q1,q3) surface for grid '//itoc(ng)//
     2                 '" to lamdat',nrt2d,z(eig),0,' ')
            call iosys('write real "surface functions for right '//
     1                 '(q1,q3) surface for grid '//itoc(ng)//
     2                 '" to lamdat',nsq*nrt2d,z(ham),0,' ')
            title='2d surface functions in left (q1,q3) plane'
            call set2x(z(eig1),z(eig3),endpts(1,2),z(ham),
     1                 z(eig),z(vec),z(t1),z(t2),z(locsc2),
     2                 ia(indx),ia(ipvt),nply(ng),npts(ng),nsq,nrt2d,
     3                 pottyp,prn2h,prn2v,title)
            call iosys('write real "surface eigenvalues for left '//
     1                 '(q1,q3) surface for grid '//itoc(ng)//
     2                 '" to lamdat',nrt2d,z(eig),0,' ')
            call iosys('write real "surface functions for left '//
     1                 '(q1,q3) surface for grid '//itoc(ng)//
     2                 '" to lamdat',nsq*nrt2d,z(ham),0,' ')
            title='2d surface functions in right (q2,q3) plane'
            call set2x(z(eig2),z(eig3),endpts(2,1),z(ham),
     1                 z(eig),z(vec),z(t1),z(t2),z(locsc2),
     2                 ia(indx),ia(ipvt),nply(ng),npts(ng),nsq,nrt2d,
     3                 pottyp,prn2h,prn2v,title)
            call iosys('write real "surface eigenvalues for right '//
     1                 '(q2,q3) surface for grid '//itoc(ng)//
     2                 '" to lamdat',nrt2d,z(eig),0,' ')
            call iosys('write real "surface functions for right '//
     1                 '(q2,q3) surface for grid '//itoc(ng)//
     2                 '" to lamdat',nsq*nrt2d,z(ham),0,' ')
            title='2d surface functions in left (q2,q3) plane'
            call set2x(z(eig2),z(eig3),endpts(1,1),z(ham),
     1                 z(eig),z(vec),z(t1),z(t2),z(locsc2),
     2                 ia(indx),ia(ipvt),nply(ng),npts(ng),nsq,nrt2d,
     3                 pottyp,prn2h,prn2v,title)
            call iosys('write real "surface eigenvalues for left '//
     1                 '(q2,q3) surface for grid '//itoc(ng)//
     2                 '" to lamdat',nrt2d,z(eig),0,' ')
            call iosys('write real "surface functions for left '//
     1                 '(q2,q3) surface for grid '//itoc(ng)//
     2                 '" to lamdat',nsq*nrt2d,z(ham),0,' ')
c     
c     solve for the three dimensional eigenfunctions and eigenvalues
c
            ncube= nply(ng)*nply(ng)*nply(ng)
            call setind(ia(indx),nply(ng),ncube,3)
            call set3x(z(eig1),z(eig2),z(eig3),z(ham),z(eig),
     1                 z(vec),z(t1),z(t2),z(t3),z(locsc2),ia(indx),
     2                 ia(ipvt),nply(ng),npts(ng),ncube,nrt3d,pottyp,
     3                 prn3h,prn3v)
            call iosys('write real "3d eigenvalues for grid '
     1                 //itoc(ng)//'" to lamdat',nrt3d,z(eig),0,' ')
            call iosys('write real "3d functions for grid '
     1                 //itoc(ng)//'" to lamdat',ncube*nrt3d,
     2                 z(ham),0,' ')
         else
            call setind(ia(indx),nply(ng),nsq,2)
            title='right (q1,q2)'
            call hamil(z(eig1),z(eig2),dum,z(t1),z(t2),dum,
     1                 endpts(2,3),z(hbuf),ia(ibuf),ia(indx),
     2                 z(diag),lenbuf,nply(ng),2,nsq,pottyp,
     3                 prn2h,title,ntot,incore)
            if(ng.eq.1.and..not.useone) then
               call icopy(ia(indx),ia(ipvt),2*nsq)
               call zguess(z(eig1),z(eig2),dum,z(t1),z(t2),dum,z(ttmp),
     1                     z(vec),z(hvec),z(eig),ia(ipvt),z(locsc1),
     2                     nply(ng),2,nsq,nrt2d,pottyp,prgues)
            else
               call guess(z(vec),z(eig),nsq,nrt2d)
            endif                    
            call david(z(hbuf),ia(ibuf),z(diag),z(eig),z(vec),
     1                 z(hvec),z(dvdmat),z(dvdvec),
     2                 z(cvn),z(hmev),z(eig0),ia(ipvt),thresh,
     3                 cnverg,ops,nsq,nrt2d,nattim,maxvec,
     4                 niter,lenbuf,ntot,incore)
            call iosys('write real "surface eigenvalues for right '//
     1                 '(q1,q2) surface for grid '//itoc(ng)//
     2                 '" to lamdat',nrt2d,z(eig),0,' ')            
            call iosys('write real "surface functions for right '//
     1                 '(q1,q2) surface for grid '//itoc(ng)//
     2                 '" to lamdat',nsq*nrt2d,z(vec),0,' ')            
            title='left (q1,q2)'
            call hamil(z(eig1),z(eig2),dum,z(t1),z(t2),dum,
     1                 endpts(1,3),z(hbuf),ia(ibuf),ia(indx),
     2                 z(diag),lenbuf,nply(ng),2,nsq,pottyp,
     3                 prn2h,title,ntot,incore)
            if(ng.eq.1.and..not.useone) then
               call icopy(ia(indx),ia(ipvt),2*nsq)
               call zguess(z(eig1),z(eig2),dum,z(t1),z(t2),dum,z(ttmp),
     1                     z(vec),z(hvec),z(eig),ia(ipvt),z(locsc1),
     2                     nply(ng),2,nsq,nrt2d,pottyp,prgues)
            else
               call guess(z(vec),z(eig),nsq,nrt2d)
            endif                    
            call david(z(hbuf),ia(ibuf),z(diag),z(eig),z(vec),
     1                 z(hvec),z(dvdmat),z(dvdvec),
     2                 z(cvn),z(hmev),z(eig0),ia(ipvt),thresh,
     3                 cnverg,ops,nsq,nrt2d,nattim,maxvec,
     4                 niter,lenbuf,ntot,incore)
            call iosys('write real "surface eigenvalues for left '//
     1                 '(q1,q2) surface for grid '//itoc(ng)//
     2                 '" to lamdat',nrt2d,z(eig),0,' ')            
            call iosys('write real "surface functions for left '//
     1                 '(q1,q2) surface for grid '//itoc(ng)//
     2                 '" to lamdat',nsq*nrt2d,z(vec),0,' ')            
            title='right (q1,q3)'
            call hamil(z(eig1),z(eig3),dum,z(t1),z(t3),dum,
     1                 endpts(2,2),z(hbuf),ia(ibuf),ia(indx),
     2                 z(diag),lenbuf,nply(ng),2,nsq,pottyp,
     3                 prn2h,title,ntot,incore)
            if(ng.eq.1.and..not.useone) then
               call icopy(ia(indx),ia(ipvt),2*nsq)
               call zguess(z(eig1),z(eig3),dum,z(t1),z(t3),dum,z(ttmp),
     1                     z(vec),z(hvec),z(eig),ia(ipvt),z(locsc1),
     2                     nply(ng),2,nsq,nrt2d,pottyp,prgues)
            else
               call guess(z(vec),z(eig),nsq,nrt2d)
            endif                    
            call david(z(hbuf),ia(ibuf),z(diag),z(eig),z(vec),
     1                 z(hvec),z(dvdmat),z(dvdvec),
     2                 z(cvn),z(hmev),z(eig0),ia(ipvt),thresh,
     3                 cnverg,ops,nsq,nrt2d,nattim,maxvec,
     4                 niter,lenbuf,ntot,incore)
            call iosys('write real "surface eigenvalues for right '//
     1                 '(q1,q3) surface for grid '//itoc(ng)//
     2                 '" to lamdat',nrt2d,z(eig),0,' ')            
            call iosys('write real "surface functions for right '//
     1                 '(q1,q3) surface for grid '//itoc(ng)//
     2                 '" to lamdat',nsq*nrt2d,z(vec),0,' ')            
            title='left (q1,q3)'
            call hamil(z(eig1),z(eig3),dum,z(t1),z(t3),dum,
     1                 endpts(1,2),z(hbuf),ia(ibuf),ia(indx),
     2                 z(diag),lenbuf,nply(ng),2,nsq,pottyp,
     3                 prn2h,title,ntot,incore)
            if(ng.eq.1.and..not.useone) then
               call icopy(ia(indx),ia(ipvt),2*nsq)
               call zguess(z(eig1),z(eig3),dum,z(t1),z(t3),dum,z(ttmp),
     1                     z(vec),z(hvec),z(eig),ia(ipvt),z(locsc1),
     2                     nply(ng),2,nsq,nrt2d,pottyp,prgues)
            else
               call guess(z(vec),z(eig),nsq,nrt2d)
            endif                    
            call david(z(hbuf),ia(ibuf),z(diag),z(eig),z(vec),
     1                 z(hvec),z(dvdmat),z(dvdvec),
     2                 z(cvn),z(hmev),z(eig0),ia(ipvt),thresh,
     3                 cnverg,ops,nsq,nrt2d,nattim,maxvec,
     4                 niter,lenbuf,ntot,incore)
            call iosys('write real "surface eigenvalues for left '//
     1                 '(q1,q3) surface for grid '//itoc(ng)//
     2                 '" to lamdat',nrt2d,z(eig),0,' ')            
            call iosys('write real "surface functions for left '//
     1                 '(q1,q3) surface for grid '//itoc(ng)//
     2                 '" to lamdat',nsq*nrt2d,z(vec),0,' ')            
            title='right (q2,q3)'
            call hamil(z(eig2),z(eig3),dum,z(t2),z(t3),dum,
     1                 endpts(2,1),z(hbuf),ia(ibuf),ia(indx),
     2                 z(diag),lenbuf,nply(ng),2,nsq,pottyp,
     3                 prn2h,title,ntot,incore)
            if(ng.eq.1.and..not.useone) then
               call icopy(ia(indx),ia(ipvt),2*nsq)
               call zguess(z(eig2),z(eig3),dum,z(t2),z(t3),dum,z(ttmp),
     1                     z(vec),z(hvec),z(eig),ia(ipvt),z(locsc1),
     2                     nply(ng),2,nsq,nrt2d,pottyp,prgues)
            else
               call guess(z(vec),z(eig),nsq,nrt2d)
            endif                    
            call david(z(hbuf),ia(ibuf),z(diag),z(eig),z(vec),
     1                 z(hvec),z(dvdmat),z(dvdvec),
     2                 z(cvn),z(hmev),z(eig0),ia(ipvt),thresh,
     3                 cnverg,ops,nsq,nrt2d,nattim,maxvec,
     4                 niter,lenbuf,ntot,incore)
            call iosys('write real "surface eigenvalues for right '//
     1                 '(q2,q3) surface for grid '//itoc(ng)//
     2                 '" to lamdat',nrt2d,z(eig),0,' ')            
            call iosys('write real "surface functions for right '//
     1                 '(q2,q3) surface for grid '//itoc(ng)//
     2                 '" to lamdat',nsq*nrt2d,z(vec),0,' ')            
            title='left (q2,q3)'
            call hamil(z(eig2),z(eig3),dum,z(t2),z(t3),dum,
     1                 endpts(1,1),z(hbuf),ia(ibuf),ia(indx),
     2                 z(diag),lenbuf,nply(ng),2,nsq,pottyp,
     3                 prn2h,title,ntot,incore)
            if(ng.eq.1.and..not.useone) then
               call icopy(ia(indx),ia(ipvt),2*nsq)
               call zguess(z(eig2),z(eig3),dum,z(t2),z(t3),dum,z(ttmp),
     1                     z(vec),z(hvec),z(eig),ia(ipvt),z(locsc1),
     2                     nply(ng),2,nsq,nrt2d,pottyp,prgues)
            else
               call guess(z(vec),z(eig),nsq,nrt2d)
            endif                    
            call david(z(hbuf),ia(ibuf),z(diag),z(eig),z(vec),
     1                 z(hvec),z(dvdmat),z(dvdvec),
     2                 z(cvn),z(hmev),z(eig0),ia(ipvt),thresh,
     3                 cnverg,ops,nsq,nrt2d,nattim,maxvec,
     4                 niter,lenbuf,ntot,incore)
            call iosys('write real "surface eigenvalues for left '//
     1                 '(q2,q3) surface for grid '//itoc(ng)//
     2                 '" to lamdat',nrt2d,z(eig),0,' ')            
            call iosys('write real "surface functions for left '//
     1                 '(q2,q3) surface for grid '//itoc(ng)//
     2                 '" to lamdat',nsq*nrt2d,z(vec),0,' ')            
            title='three dimensional box'
            ncube= nply(ng)*nply(ng)*nply(ng)
            call setind(ia(indx),nply(ng),ncube,3)
            call hamil(z(eig1),z(eig2),z(eig3),z(t1),z(t2),z(t3),
     1                 dum,z(hbuf),ia(ibuf),ia(indx),
     2                 z(diag),lenbuf,nply(ng),3,ncube,pottyp,
     3                 prn3h,title,ntot,incore)
            if(ng.eq.1.and..not.useone) then
               call icopy(ia(indx),ia(ipvt),3*ncube)
               call zguess(z(eig1),z(eig2),z(eig3),z(t1),z(t2),
     1                     z(t3),z(ttmp),z(vec),z(hvec),z(eig),
     2                     ia(ipvt),z(locsc1),nply(ng),3,ncube,nrt3d,
     3                     pottyp,prgues)
            else
               call guess(z(vec),z(eig),ncube,nrt3d)
            endif                    
            call david(z(hbuf),ia(ibuf),z(diag),z(eig),z(vec),
     1                 z(hvec),z(dvdmat),z(dvdvec),
     2                 z(cvn),z(hmev),z(eig0),ia(ipvt),thresh,
     3                 cnverg,ops,ncube,nrt3d,nattim,maxvec,
     4                 niter,lenbuf,ntot,incore)
            call iosys('write real "3d eigenvalues for grid '
     1                 //itoc(ng)//'" to lamdat',nrt3d,z(eig),0,' ')
            call iosys('write real "3d functions for grid '
     1                 //itoc(ng)//'" to lamdat',ncube*nrt3d,
     2                 z(vec),0,' ')
         endif     
c
c           project the eigenfunctions onto the surface states
c
         call iosys('read real "3d functions for grid '
     1              //itoc(ng)//'" from lamdat',ncube*nrt3d,
     2              z(vec3d),0,' ')
         call iosys('read real "surface functions for right '//
     1              '(q1,q2) surface for grid '//itoc(ng)//
     2              '" from lamdat',nsq*nrt2d,z(vec2d),0,' ')
         offset=0
         call filprj(z(p3),z(locsc1),'r',nply(ng),npts(ng))
         call surf(z(vec3d),z(vec2d),z(locsc1),z(prj),
     1             z(srf+offset),ia(indx),nply(ng),ncube,nrt3d,nsq,
     2             nrt2d,'q3',prnprj,prn2v,prn3v)
         call iosys('write real "surface projections for right '//
     1              '(q1,q2) surface for grid '//itoc(ng)//
     2              '" to lamdat',nrt3d*nrt2d,z(srf+offset),0,' ')
         offset=offset+nrt2d*nrt3d
         call iosys('read real "surface functions for left '//
     1              '(q1,q2) surface for grid '//itoc(ng)//
     2              '" from lamdat',nsq*nrt2d,z(vec2d),0,' ')
         call filprj(z(p3),z(locsc1),'l',nply(ng),npts(ng))
         call surf(z(vec3d),z(vec2d),z(locsc1),z(prj),
     1             z(srf+offset),ia(indx),nply(ng),ncube,nrt3d,nsq,
     2             nrt2d,'q3',prnprj,prn2v,prn3v)
         call iosys('write real "surface projections for left '//
     1              '(q1,q2) surface for grid '//itoc(ng)//
     2              '" to lamdat',nrt3d*nrt2d,z(srf+offset),0,' ')
         offset=offset+nrt2d*nrt3d
         call iosys('read real "surface functions for right '//
     1              '(q1,q3) surface for grid '//itoc(ng)//
     2              '" from lamdat',nsq*nrt2d,z(vec2d),0,' ')
         call filprj(z(p2),z(locsc1),'r',nply(ng),npts(ng))
         call surf(z(vec3d),z(vec2d),z(locsc1),z(prj),
     1             z(srf+offset),ia(indx),nply(ng),ncube,nrt3d,nsq,
     2             nrt2d,'q2',prnprj,prn2v,prn3v)
         call iosys('write real "surface projections for right '//
     1              '(q1,q3) surface for grid '//itoc(ng)//
     2              '" to lamdat',nrt3d*nrt2d,z(srf+offset),0,' ')
         offset=offset+nrt2d*nrt3d
         call iosys('read real "surface functions for left '//
     1              '(q1,q3) surface for grid '//itoc(ng)//
     2              '" from lamdat',nsq*nrt2d,z(vec2d),0,' ')
         call filprj(z(p2),z(locsc1),'l',nply(ng),npts(ng))
         call surf(z(vec3d),z(vec2d),z(locsc1),z(prj),
     1             z(srf+offset),ia(indx),nply(ng),ncube,nrt3d,nsq,
     2             nrt2d,'q2',prnprj,prn2v,prn3v)
         call iosys('write real "surface projections for left '//
     1              '(q1,q3) surface for grid '//itoc(ng)//
     2              '" to lamdat',nrt3d*nrt2d,z(srf+offset),0,' ')
         offset=offset+nrt2d*nrt3d
         call iosys('read real "surface functions for right '//
     1              '(q2,q3) surface for grid '//itoc(ng)//
     2              '" from lamdat',nsq*nrt2d,z(vec2d),0,' ')
         call filprj(z(p1),z(locsc1),'r',nply(ng),npts(ng))
         call surf(z(vec3d),z(vec2d),z(locsc1),z(prj),
     1             z(srf+offset),ia(indx),nply(ng),ncube,nrt3d,nsq,
     2             nrt2d,'q1',prnprj,prn2v,prn3v)
         call iosys('write real "surface projections for right '//
     1              '(q2,q3) surface for grid '//itoc(ng)//
     2              '" to lamdat',nrt3d*nrt2d,z(srf+offset),0,' ')
         offset=offset+nrt2d*nrt3d
         call iosys('read real "surface functions for left '//
     1              '(q2,q3) surface for grid '//itoc(ng)//
     2              '" from lamdat',nsq*nrt2d,z(vec2d),0,' ')
         call filprj(z(p1),z(locsc1),'l',nply(ng),npts(ng))
         call surf(z(vec3d),z(vec2d),z(locsc1),z(prj),
     1             z(srf+offset),ia(indx),nply(ng),ncube,nrt3d,nsq,
     2             nrt2d,'q1',prnprj,prn2v,prn3v)
         call iosys('write real "surface projections for left '//
     1              '(q2,q3) surface for grid '//itoc(ng)//
     2              '" to lamdat',nrt3d*nrt2d,z(srf+offset),0,' ')
         call iosys('read real "3d eigenvalues for grid '
     1              //itoc(ng)//'" from lamdat',nrt3d,z(eig),0,' ')
         if(nen.ne.0) then
            if(posinp('$energy',cpass)) then
               call cardin(card)
            endif
            call fparr(card,'energies',z(energy),nen,' ')
            offsti=0
            do 60 ic=1,6
               offstj=0
               do 70 jc=1,ic
                  call rmat(z(srf+offsti),z(srf+offstj),z(eig),
     1                      z(energy),z(rmt),nrt2d,nrt3d,nen,ic,
     2                      jc,prnrmt)
                  offstj=offstj+nrt2d*nrt3d
 70            continue   
               offsti=offsti+nrt2d*nrt3d
 60         continue   
         endif
         eig1=eig1+nply(ng)
         eig2=eig2+nply(ng)
         eig3=eig3+nply(ng)
         cord1=cord1+npts(ng)
         cord2=cord2+npts(ng)
         cord3=cord3+npts(ng)
         wt1=wt1+npts(ng)
         wt2=wt2+npts(ng)
         wt3=wt3+npts(ng)
         p1=p1+npts(ng)*nply(ng)
         p2=p2+npts(ng)*nply(ng)
         p3=p3+npts(ng)*nply(ng)   
         t1=t1+nply(ng)*nply(ng)
         t2=t2+nply(ng)*nply(ng)
         t3=t3+nply(ng)*nply(ng)
 50   continue      
c               davidson diagonalization      
c
c      else
c             call lnkerr('error in diagonalization technique')
c      endif
      call iosys ('write integer maxsiz to rwf',1,canget,0,' ')
      call chainx(0)               
      stop
 1    format(/,5x,'minimum order of polynomials in each dimension = '
     1                                                       ,i3,/,5x,
     2            'number of grids                                = '
     3                                                       ,i3,/,5x,
     4            'x left                                         = '
     5                                                       ,e15.8,/,
     6                                                             5x,
     7            'x right                                        = '
     8                                                       ,e15.8,/,
     9                                                             5x,
     x            'y left                                         = '
     x                                                       ,e15.8,/,
     x                                                             5x,
     x            'y right                                        = '
     x                                                       ,e15.8,/,
     x                                                             5x,
     x            'z left                                         = '
     x                                                       ,e15.8,/,
     x                                                             5x,
     x            'z right                                        = '
     x                                                       ,e15.8,/,
     x                                                             5x,
     x            'derivative at left boundary                    = '
     x                                                       ,e15.8,/,
     x                                                             5x,
     x            'derivative at right boundary                   = '
     x                                                       ,e15.8,/,
     x                                                             5x,
     x            'order of leading left polynomial               = '
     x                                                       ,i1,/,5x,
     x            'order of leading right polynomial              = '
     x                                                       ,i1,/,5x,
     x            'potential type                                 = '
     x                                                       ,a16,/5x,)
 2    format(/,5x,'davidson iterative diagonalization',/,5x,
     1            'number of 2d-roots                       = ',i3,/,5x,
     2            'number of 3d-roots                       = ',i3,/,5x,
     3            'number calculated at a time              = ',i3,/,5x,
     4            'maximum number of vectors stored in core = ',i3,/,5x,
     5            'maximum number of iterations             = ',i3,/,5x,
     6            'threshold criterion for vectors          = ',e15.8,
     7                                                           /,5x,
     8            'convergence criterion                    = ',e15.8,
     9                                                          /,5x,
     x            'buffer size                              = ',i8)
 3    format(/,5x,'direct diagonalization',/,5x,
     1            'number of 2-d roots = ',i3,/,5x,
     2            'number of 3-d roots = ',i3)
 99   format(/,'requesting ',i8,' working precision words')
      end
