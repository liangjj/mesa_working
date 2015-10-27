*deck pde.f 
c***begin prologue     pde
c***date written       970810   (yymmdd)
c***revision date               (yymmdd)
c***keywords           pde
c***                   
c***author             schneider, b. i.(nsf)
c***source             pde
c***purpose            three dimensional inhomogeneous pde code
c
c***description        the one, two or three dimensional inhomogeneous pde
c***                   problem is solved using a discrete variable
c***                   representation. 
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       pde
      program pde
c
      implicit integer (a-z)
      parameter ( ngmax=10 )
      character*4096 ops
      character*8 prtflg
      character*2 itoc
      character*80 cpass, title, chrkey
      character*24 mattyp
      character*1600 card
      character*128 filham
      character*2 atom
      character*16 coord, fptoc
      character*8 qtyp, qdtyp
      logical dollar, logkey, prnply, prnwpt, prnh 
      logical prnht, prall, chko, check, nfix, prnh0
      logical toau, useau, hamd, incore, cgrid
      logical itsolv, prbufh
      logical prlin, prital, frmdsk
      logical bypass, scatt, typot
      real*8 z, fpkey, amass, omegat, left, right, energy 
      real*8 pi, dleft, dright
      real*8 thresh, cnverg, hbar
      real*8 massau, lenau, timau, dum, scale, sfac
      complex*16 cdum
      common z(1)
      dimension ia(1), nmax(ngmax,3), npt(ngmax,3)
      dimension lftbc(3), rtbc(3), dim(3)
      dimension prnply(3), prnwpt(3), prnh(3), check(3)
      dimension q(ngmax,3), wt(ngmax,3), p(ngmax,3), dp(ngmax,3)
      dimension ddp(ngmax,3), pn(ngmax,3), dpn(ngmax,3), ddpn(ngmax,3)
      dimension eigc(ngmax,3), wtc(ngmax,3), work(4)
      dimension omegat(3), eig(3), left(3), right(3)
      dimension dleft(3), dright(3)
      dimension ham(3), v(3), u(3), vecl(3), vecr(3)
      dimension prlin(11), qtyp(3), qdtyp(3)
      dimension sfac(3), energy(100)
      dimension nfix(2)
      equivalence (ia(1),z(1))
      common/io/inp, iout      
      common/memory/ioff
      data pi/3.14159265358979323844d0/
c     hbar in joule-sec      
      data hbar/1.054592d-34/
c
c          mass in Kilograms length in meters
      data massau, lenau, timau / 9.109558d-31, 5.291771d-11, 
     1                            2.418884d-17 /
      call drum
      write(iout,*)
      write(iout,*) '    three dimensional pde eigenvalue code       '
      call iosys ('read character options from rwf',-1,0,0,ops)
      call iosys ('read character "print flag" from rwf',-1,0,0,prtflg)
c
c                       set options
c
      frmdsk=logkey(ops,'from-disk',.false.,' ')   
      ngrid=intkey(ops,'number-of-grids',1,' ')
      nrhs=intkey(ops,'number-of-right-hand-sides',1,' ')
      prrhs=logkey(ops,'print-right-hand-side',.false.,' ')
      ug=intkey(ops,'grid',1,' ')
      numdim=intkey(ops,'number-of-dimensions',1,' ')   
      coord=chrkey(ops,'coordinate-system','cartesian',' ')
      scatt=logkey(ops,'scattering-calculation',.false.,' ')
      typot=logkey(ops,'coulomb-potential',.false.,' ')
      mattyp=chrkey(ops,'matrix-type','real-symmetric',' ')
      prall=logkey(ops,'print=m7040=all',.false.,' ')
      bypass=logkey(ops,'bypass=on',.false.,' ')
c----------------------------------------------------------------------c
c   1. if task = diagonalize-h0, the code will diagonalize a hamiltonian
c      consisting of the kinetic energy plus one and two body local 
c      interactions.  the interaction terms are not assumed to depend on
c      the wavefunction.                                           
c----------------------------------------------------------------------c 
      write(iout,1) coord
      write(iout,2)
      if(bypass) then
         call genmat(numdim,itsolv,mattyp)
         call chainx(0)               
         stop
      endif                  
      if(prall) then
         prnply(1)=.true.
         prnply(2)=.true.
         prnply(3)=.true.
         prnwpt(1)=.true. 
         prnwpt(2)=.true.
         prnwpt(3)=.true.
         prnh0=.true. 
         prnht=.true. 
         prnh(1)=.true.
         prnh(2)=.true.
         prnh(3)=.true.
      else         
         prnply(1)=logkey(ops,'print=m7040=q1-polynomials',.false.,' ')
         prnply(2)=logkey(ops,'print=m7040=q2-polynomials',.false.,' ')
         prnply(3)=logkey(ops,'print=m7040=q3-polynomials',.false.,' ')
         prnwpt(1)=logkey(ops,'print=m7040=q1-points/weights',
     1                    .false.,' ')
         prnwpt(2)=logkey(ops,'print=m7040=q2-points/weights',
     1                    .false.,' ')
         prnwpt(3)=logkey(ops,'print=m7040=q3-points/weights',
     1                    .false.,' ')
         prnh0=logkey(ops,'print=m7040=h0',.false.,' ')
         prnht=logkey(ops,'print=m7040=h',.false.,' ')
         prnh(1)=logkey(ops,'print=m7040=q1-hamiltonian',.false.,' ')
         prnh(2)=logkey(ops,'print=m7040=q2-hamiltonian',.false.,' ')
         prnh(3)=logkey(ops,'print=m7040=q3-hamiltonian',.false.,' ')
      endif
      chko=logkey(ops,'m7040=check-orthogonality',.false.,' ')
      if(chko) then
         check(1)=.true.
         check(2)=.true.
         check(3)=.true.
      else                  
         check(1)=logkey(ops,'check-q1-orthogonality',.false.,' ')
         check(2)=logkey(ops,'check-q2-orthogonality',.false.,' ')
         check(3)=logkey(ops,'check-q3-orthogonality',.false.,' ')
      endif
      nen=intkey(ops,'number-of-energies',1,' ')
      call fparr(ops,'energies',energy,nen,' ')
      toau=logkey(ops,'to-atomic-units',.false.,' ')
      useau=logkey(ops,'use-atomic-units',.false.,' ')
      if(numdim.eq.1) then
         write(iout,3) prnply(1),
     1                 prnwpt(1), prnh0,
     2                 prnh(1), 
     3                 check(1)
      endif
      if(numdim.eq.2) then
         write(iout,4) prnply(1), prnply(2),
     1                 prnwpt(1), prnwpt(2), prnh0,
     2                 prnh(1), prnh(2), 
     3                 check(1), check(2)
      endif     
      if(numdim.eq.3) then
         write(iout,5) prnply(1), prnply(2), prnply(3),
     1                 prnwpt(1), prnwpt(2), prnwpt(3), prnh0,
     2                 prnh(1), prnh(2), prnh(3), 
     3                 check(1), check(2), check(3)
      endif                     
c
c              set various program options
c      
c
c         set trap configuration, numbers of atoms, various constants
c         and other physical parameters
c       
      if ( dollar('$trap',card,cpass,inp) ) then
           atom=chrkey(card,'atom','Cs',' ')
           call atmdat(atom,amass,omegat)
           if(numdim.eq.1) then
              write(iout,6) amass, omegat(1)
           elseif(numdim.eq.2) then
              write(iout,7) amass, (omegat(i),i=1,2)
           elseif(numdim.eq.3) then
              write(iout,8) amass, (omegat(i),i=1,3)
           endif
      endif
c
c               spatial basis set information
c       
      if ( dollar('$q1-polynomials',card,cpass,inp) ) then
           call bfdat(frmdsk,card,qtyp(1),qdtyp(1),left(1),right(1),
     1                dleft(1),dright(1),lftbc(1),rtbc(1),
     2                npt(1,1),nmax(1,1),nfix,ngrid)

      endif     
      if(numdim.gt.1) then
         if ( dollar('$q2-polynomials',card,cpass,inp) ) then
              call bfdat(frmdsk,card,qtyp(2),qdtyp(2),left(2),right(2),
     1                   dleft(2),dright(2),lftbc(2),rtbc(2),
     2                   npt(1,2),nmax(1,2),nfix,ngrid)
         endif
      endif
      if(numdim.gt.2) then
         if ( dollar('$q3-polynomials',card,cpass,inp) ) then
              call bfdat(frmdsk,card,qtyp(3),qdtyp(3),left(3),right(3),
     1                   dleft(3),dright(3),lftbc(3),rtbc(3),
     2                   npt(1,3),nmax(1,3),nfix,ngrid)
         endif
      endif
      if(toau) then
c
c        converts everything to atomic units
c
         write(iout,*) 'converting to atomic units'
         hbar=1.d0
         scatl=scatl/lenau
         do 10 i=1,numdim
            left(i)=left(i)/lenau
            right(i)=right(i)/lenau
            omegat(i)=omegat(i)*timau
 10      continue
         amass=amass/massau        
      endif
      if(useau) then
c
c        if this option is used then all is entered in atomic units
c
         write(iout,*) 'assuming atomic units'
         hbar=1.d0
         amass=1.d0
         do 20 i=1,numdim
            omegat(i)=1.d0
 20      continue            
      endif
      do 30 i=1,numdim
         dim(i)=0
         do 40 j=1,ngrid
            dim(i)=max(dim(i),npt(j,i))            
 40      continue
 30   continue
      maxd=max(dim(1),dim(2),dim(3))
      n3d=1
      do 50 i=1,numdim
         n3d=n3d*nmax(ug,i)
 50   continue 
c
c             set system solver procedures
c      
      itsolv=logkey(ops,'iterative-linear-solve',.false.,' ')
      if(itsolv) then
         call lindat(card,cpass,n3d,nrhs,nattim,cnverg,thresh,
     1               niter,nvec,lenbuf,cgrid,prbufh,hamd,
     2               prlin,prital,filham)
         write(iout,9) nrhs, thresh, cnverg, niter, nvec
      else
         write(iout,11) n3d, nrhs
      endif
      ioff=1
c
c     note the alignment of variables.  it is important in passing these 
c     arrays later.
c 
      do 70 i=1,2
         call wrdcnt(q(ug,1),wt(ug,1),eigc(ug,1),wtc(ug,1),
     1               p(ug,1),dp(ug,1),ddp(ug,1),
     2               pn(ug,1),dpn(ug,1),ddpn(ug,1),
     3               ham(1),eig(1),v(1),u(1),vecl(1),vecr(1),
     4               mattyp,dim,ngrid,ioff)
         if(numdim.gt.1) then
            call wrdcnt(q(ug,2),wt(ug,2),eigc(ug,2),wtc(ug,2),
     1                  p(ug,2),dp(ug,2),ddp(ug,2),
     2                  pn(ug,2),dpn(ug,2),ddpn(ug,2),
     3                  ham(2),eig(2),v(2),u(2),vecl(2),vecr(2),
     4                  mattyp,dim,ngrid,ioff)
         endif
         if(numdim.gt.2) then
            call wrdcnt(q(ug,3),wt(ug,3),eigc(ug,3),wtc(ug,3),
     1                  p(ug,3),dp(ug,3),ddp(ug,3),
     2                  pn(ug,3),dpn(ug,3),ddpn(ug,3),
     3                  ham(3),eig(3),v(3),u(3),vecl(3),vecr(3),
     4                  mattyp,dim,ngrid,ioff)
         endif
         ind=wpadti(ioff)
         idum=ind+(numdim+1)*n3d
         words=iadtwp(idum+numdim*n3d)            
         space=words    
         work(1)=words
         if(mattyp.eq.'complex'.or.
     1             mattyp.eq.'real-unsymmetric') then
            lwork=20*maxd
            work(2)=work(1)+max(4*maxd*maxd,6*n3d,lwork)
            words=work(2)+4*maxd*maxd
         elseif(mattyp.eq.'real-symmetric') then
            lwork=0
            work(2)=work(1)+max(2*maxd*maxd,3*n3d)
            words=work(2)+2*maxd*maxd
         endif                  
         if(.not.itsolv) then
            hamtot=words
            rhs=hamtot+n3d*n3d*str
            vtot=rhs+n3d*str*nrhs
            ipvt=wpadti(vtot+n3d*str)
            words=ipvt+n3d
         else
             hbuf=words
             ibuf=wpadti(hbuf+str*lenbuf)
             diag=iadtwp(ibuf+2*lenbuf)
             vtot=diag+n3d*str
             trials=vtot+n3d*str
             pvec=trials+n3d*nrhs*str
             hpvec=pvec+n3d*nvec*str
             vec=hpvec+nvec*n3d*str
             h=vec+n3d*nvec*str
             htmp=h+nvec*nvec*str
             b=htmp+nvec*nvec*str
             btmp=b+nvec*nvec*str
             svec=btmp+nvec*nvec*str
             sveca=svec+nvec*nvec*str
             eigtot=sveca+nvec*nvec*str
             etmp=eigtot+n3d*str
             work(3)=etmp+n3d*str
             mwork=10*nvec*str
             work(4)=work(3)+max(n3d*str*nvec,mwork)
             words=work(4)+max(n3d*str*nvec,2*str*nvec)
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
             call getscm(words,z,ngot,'poly',0)
         endif
   70 continue
      lwork=max(lwork,mwork)
c
c     in order to proceed it is necessary to be able to calculate integrals
c     between the orthogonal polynomials to be generated on the desired
c     integration interval.  the calculation of the points and weights is
c     done on the interval ( -1.,1. ) and then converted to ( left, right ).
c     the fixed points on the standard interval are the two end points.      
c
      do 80 i=1,numdim
c      
         write(iout,12) i, ug
         if(.not.frmdsk) then
             call getqpt(z(q(ug,i)),z(wt(ug,i)),left(i),right(i),
     1                   qdtyp(i),'before',z(work(1)),
     2                   nfix,npt(ug,i),npt(ug,i),1,.false.)
             if(qdtyp(i).eq.'hermite') then
                call smul(z(q(ug,i)),z(q(ug,i)),sfac(i),npt(ug,i))
                call smul(z(wt(ug,i)),z(wt(ug,i)),sfac(i),npt(ug,i))
             endif   
             call lgrply(z(p(ug,i)),z(dp(ug,i)),z(ddp(ug,i)),z(q(ug,i)),
     1                   z(q(ug,i)),npt(ug,i),npt(ug,i),.false.)
             call prepfn(z(p(ug,i)),z(dp(ug,i)),z(ddp(ug,i)),z(q(ug,i)),
     1                   z(wt(ug,i)),npt(ug,i),sfac(i),qdtyp(i))
             call copy(z(p(ug,i)),z(pn(ug,i)),npt(ug,i)*npt(ug,i))
             call copy(z(dp(ug,i)),z(dpn(ug,i)),npt(ug,i)*npt(ug,i))
             call copy(z(ddp(ug,i)),z(ddpn(ug,i)),npt(ug,i)*npt(ug,i))
             if(prnply(i)) then                  
                title='dvr polynomials '//itoc(i)//' dimension'
                call prntfm(title,z(pn(ug,i)),npt(ug,i),npt(ug,i),
     1                      npt(ug,i),npt(ug,i),iout)
                title='first derivative of dvr polynomials '
     1                 //itoc(i)//' dimension'
                call prntfm(title,z(dpn(ug,i)),npt(ug,i),npt(ug,i),
     1                      npt(ug,i),npt(ug,i),iout)
                title='second derivative of dvr polynomials '
     1                 //itoc(i)//' dimension'
                call prntfm(title,z(ddpn(ug,i)),npt(ug,i),npt(ug,i),
     1                      npt(ug,i),npt(ug,i),iout)
             endif
             if(check(i)) then
                call chk(z(pn(ug,i)),z(wt(ug,i)),z(work(1)),z(work(2)),
     1                   npt(ug,i),npt(ug,i))
                title='overlap matrix dvr polynomials '//itoc(i)
     1                 //' dimension'
                call prntfm(title,z(work(2)),npt(ug,i),npt(ug,i),
     1                      npt(ug,i),npt(ug,i),iout)
             endif
             if(lftbc(i).eq.0.or.rtbc(i).eq.0) then
                call dropfn(z(p(ug,i)),z(dp(ug,i)),z(ddp(ug,i)),
     1                      z(pn(ug,i)),z(dpn(ug,i)),
     2                      z(ddpn(ug,i)),z(q(ug,i)),
     3                      z(wt(ug,i)),z(eigc(ug,i)),z(wtc(ug,i)),
     4                      lftbc(i),rtbc(i),npt(ug,i),nmax(ug,i))
             endif     
         else
             call iosys('read real "points grid = '//itoc(ug)
     1                  //'" from rwf',nmax(ug,i),z(eigc(ug,i)),0,' ')
             call iosys('read real "weights grid = '//itoc(ug)
     1                 //'" from rwf',nmax(ug,i),z(wtc(ug,i)),0,' ')
             call iosys('read real "lgrbf = '//itoc(ug)//
     1                  ' grid = '//itoc(ug)//'" from rwf',
     2                   nmax(ug,i)*nmax(ug,i),z(pn(ug,i)),0,' ')
             call iosys('read real "dlgrbf = '//itoc(ug)//
     1                 ' grid = '//itoc(ug)//'" from rwf',
     2                   nmax(ug,i)*nmax(ug,i),z(dpn(ug,i)),0,' ')
             call iosys('read real "ddlgrbf = '//itoc(ug)//
     3                  ' grid = '//itoc(ug)//'" from rwf',
     4                   nmax(ug,i)*nmax(ug,i),z(ddpn(ug,i)),0,' ')
         endif              
         if(prnply(i)) then                  
            title='dvr polynomials '//itoc(i)//' dimension'
            call prntfm(title,z(pn(ug,i)),nmax(ug,i),nmax(ug,i),
     1                  nmax(ug,i),nmax(ug,i),iout)
            title='first derivative of dvr polynomials '
     1             //itoc(i)//' dimension'
            call prntfm(title,z(dpn(ug,i)),nmax(ug,i),nmax(ug,i),
     1                  nmax(ug,i),nmax(ug,i),iout)
            title='second derivative of dvr polynomials '
     1             //itoc(i)//' dimension'
            call prntfm(title,z(ddpn(ug,i)),nmax(ug,i),nmax(ug,i),
     1                  nmax(ug,i),nmax(ug,i),iout)
         endif
         if(check(i)) then
            call chk(z(pn(ug,i)),z(wtc(ug,i)),z(work(1)),z(work(2)),
     1               nmax(ug,i),nmax(ug,i),)
            title='overlap matrix dvr polynomials '//itoc(i)
     1             //' dimension'
            call prntfm(title,z(work(2)),nmax(ug,i),nmax(ug,i),
     1                  nmax(ug,i),nmax(ug,i),iout)
         endif
         call vmat1d(z(eigc(ug,i)),z(v(i)),hbar,amass,scale,nmax(ug,i),
     1               i,prnh0,ops)
c
c        the hamiltonian matrix contains the kinetic energy and any
c        one body potential.
c
         call ham0(z(pn(ug,i)),z(dpn(ug,i)),z(ddpn(ug,i)),z(eigc(ug,i)),
     1             z(wtc(ug,i)),z(ham(i)),z(ham(i)),z(v(i)),
     2             z(u(i)),z(u(i)),z(eig(i)),z(eig(i)),
     3             z(vecl(i)),z(vecl(i)),z(vecr(i)),z(vecl(i)),
     4             z(work(1)),lwork,z(work(2)),hbar,amass,nmax(ug,i),
     5             lftbc(i),rtbc(i),coord,qtyp(i),prnh0,mattyp)
         if(mattyp.eq.'complex') then
            call cvscal(z(eig(i)),z(eig(i)),scale,nmax(ug,i))
c            if(prnh0) then
               title='eigenvalues of unperturbed hamiltonian '//itoc(i)
     1                //' dimension'
               call prntcm(title,z(eig(i)),nmax(ug,i),1,nmax(ug,i),
     1                     1,iout)
c            endif
            call cvscal(z(eig(i)),z(eig(i)),1.d0/scale,nmax(ug,i))
         else
           call vscale(z(eig(i)),z(eig(i)),scale,nmax(ug,i))
c           if(prnh0) then
               title='eigenvalues of unperturbed hamiltonian '//itoc(i)
     1                //' dimension'
               call prntfm(title,z(eig(i)),nmax(ug,i),1,
     1                     nmax(ug,i),1,iout)
c            endif 
            call vscale(z(eig(i)),z(eig(i)),1.d0/scale,nmax(ug,i))
         endif
 80   continue
      call setind(ia(ind),nmax,n3d,ug,numdim,ngmax)
      if(numdim.eq.1) then
         call vmat(z(eigc(ug,1)),dum,dum,z(vtot),z(work(1)),
     1             hbar,amass,scale,n3d,nmax(ug,1),nmax(ug,2),
     2             nmax(ug,3),numdim,prnh0,ops)   
      elseif(numdim.eq.2) then
         call vmat(z(eigc(ug,1)),z(eigc(ug,2)),dum,z(vtot),z(work(1)),
     1             hbar,amass,scale,n3d,nmax(ug,1),nmax(ug,2),
     2             nmax(ug,3),numdim,prnh0,ops)   
      elseif(numdim.eq.3) then
         call vmat(z(eigc(ug,1)),z(eigc(ug,2)),z(eigc(ug,3)),z(vtot),
     1             z(work(1)),hbar,amass,scale,n3d,
     2             nmax(ug,1),nmax(ug,2),nmax(ug,3),
     3             numdim,prnh0,ops)     
      endif
c
c     subtract the one body interations included in the zeroth
c     order hamiltonian.
c
      call vpert(z(vtot),z(v(1)),z(v(2)),z(v(3)),n3d,
     1           nmax(ug,1),nmax(ug,2),nmax(ug,3),numdim)
c
c     make the right hand side.  in this example we use the R-matrix
c     equations in the DVR representation as an example.
c     note that this is only appropriate for a one-dimensional problem
c     in what follows.  the generalization is trivial but tedious for
c     more than one dimension.
c
      call calrhs(z(pn(ug,1)),z(rhs),z(rhs),z(eigc(ug,1)),z(wtc(ug,1)),
     1            n3d,mattyp,'delta',prrhs)
      if(itsolv) then
         call iosys('write integer "length of vector" '//
     1              'to ham',1,n3d,0,' ')
         if(numdim.eq.1) then
            title='buffered hamiltonian for '//itoc(numdim)
     1             //'d dimensional hamiltonian'
            call hamil(z(ham(1)),z(ham(1)),
     1                 dum,cdum,
     2                 dum,cdum,
     3                 z(vtot),z(hbuf),z(hbuf),ia(ibuf),
     4                 ia(ind),z(diag),z(diag), lenbuf,
     5                 nmax(ug,1),nmax(ug,2),nmax(ug,3),
     6                 numdim,n3d,prbufh,nel,incore,title,mattyp)
            do 100 ien=1,nen
               if(mattyp.eq.'complex') then
                  call cindvd(z(hbuf),ia(ibuf),z(diag),energy(ien),
     1                        z(eig(1)),dum,dum,
     2                        z(vecl(1)),dum,dum,
     3                        z(vecr(1)),dum,dum,     
     4                        z(pvec),z(hpvec),z(vec),z(rhs),
     5                        z(h),z(htmp),z(b),z(btmp),
     6                        z(vec),z(work(3)),z(work(4)),
     7                        scale,cnverg,thresh,
     8                        nmax(ug,1),nmax(ug,2),nmax(ug,3),
     9                        numdim,n3d,nrhs,niter,nvec,
     x                        lenbuf,hamd,incore,nel,prlin,mattyp)
                   call tocord(z(rhs),z(rhs),z(pn(ug,1)),n3d,mattyp)
                   title='solution vectors energy = '
     1                    //fptoc(energy(ien))
                   call prntcm(title,z(rhs),n3d,nrhs,n3d,nrhs,iout)
               else
                   call lindvd(z(hbuf),ia(ibuf),z(diag),energy(ien),
     1                         z(eig(1)),dum,dum,
     2                         z(vecl(1)),dum,dum,
     3                         z(pvec),z(hpvec),z(vec),z(rhs),
     4                         z(h),z(htmp),z(b),z(btmp),
     5                         z(vec),z(work(3)),z(work(4)),
     6                         scale,cnverg,thresh,
     7                         nmax(ug,1),nmax(ug,2),nmax(ug,3),
     8                         numdim,n3d,nrhs,niter,nvec,
     9                         lenbuf,hamd,incore,nel,prlin,mattyp)
                   call tocord(z(rhs),z(rhs),z(pn(ug,1)),n3d,mattyp)
                   title='solution vectors energy = '
     1                   //fptoc(energy(ien))
                   call prntrm(title,z(rhs),n3d,nrhs,n3d,nrhs,iout)
                   title='exact solution vectors energy = '
     1                   //fptoc(energy(ien))
                   call exact(z(rhs),z(rhs),z(eigc(ug,1)),n3d,
     1                        mattyp,'delta') 
                   call prntrm(title,z(rhs),n3d,nrhs,n3d,nrhs,iout)
               endif
 100        continue                         
         elseif(numdim.eq.2) then
            title='buffered hamiltonian for '//itoc(numdim)
     1             //'d dimensional hamiltonian'
            call hamil(z(ham(1)),z(ham(1)),
     1                 z(ham(2)),z(ham(2)),
     2                 dum,cdum,
     3                 z(vtot),z(hbuf),z(hbuf),ia(ibuf),
     4                 ia(ind),z(diag),z(diag), lenbuf,
     5                 nmax(ug,1),nmax(ug,2),nmax(ug,3),numdim,n3d,
     6                 prbufh,nel,incore,title,mattyp)
            do 200 ien=1,nen
               if(mattyp.eq.'complex') then
c                  call clindvd(z(hbuf),ia(ibuf),z(diag),energy(ien),
c     1                         z(eig(1)),z(eig(2)),dum,
c     2                         z(vecl(1)),z(vecl(2)),dum,
c     3                         z(vecr(1)),z(vecr(2)),dum,     
c     4                         z(pvec),z(hpvec),z(vec),z(rhs),
c     5                         z(h),z(htmp),z(b),z(btmp),
c     6                         z(vec),z(work(3)),z(work(4)),
c     7                         scale,cnverg,thresh,
c     8                         nmax(ug,1),nmax(ug,2),nmax(ug,3),
c     9                         numdim,n3d,nrhs,niter,nvec,
c     x                         lenbuf,hamd,incore,nel,prlin,mattyp)
                  title='solution vectors energy = '
     1                   //fptoc(energy(ien))
                  call prntcm(title,z(rhs),n3d,nrhs,n3d,nrhs,iout)
               else
                  call lindvd(z(hbuf),ia(ibuf),z(diag),energy(ien),
     1                        z(eig(1)),z(eig(2)),dum,
     2                        z(vecl(1)),z(vecl(2)),dum,
     3                        z(pvec),z(hpvec),z(vec),z(rhs),
     4                        z(h),z(htmp),z(b),z(btmp),
     5                        z(vec),z(work(3)),z(work(4)),
     6                        scale,cnverg,thresh,
     7                        nmax(ug,1),nmax(ug,2),nmax(ug,3),
     8                        numdim,n3d,nrhs,ntrials,niter,nvec,
     9                        lenbuf,hamd,incore,nel,prlin,mattyp)
                  title='solution vectors energy = '
     1                    //fptoc(energy(ien))
                  call prntrm(title,z(rhs),n3d,nrhs,n3d,nrhs,iout)     
               endif
 200        continue                                     
         elseif(numdim.eq.3) then
            title='buffered hamiltonian for '//itoc(numdim)
     1             //'d dimensional hamiltonian'
            call hamil(z(ham(1)),z(ham(1)),
     1                 z(ham(2)),z(ham(2)),
     2                 z(ham(3)),z(ham(3)),
     3                 z(vtot),z(hbuf),z(hbuf),ia(ibuf),
     4                 ia(ind),z(diag),z(diag), lenbuf,
     5                 nmax(ug,1),nmax(ug,2),nmax(ug,3),numdim,n3d,
     6                 prbufh,nel,incore,title,mattyp)
            do 300 ien=1,nen
               if(mattyp.eq.'complex') then
c                  call clindvd(z(hbuf),ia(ibuf),z(diag),energy(ien),
c     1                         z(eig(1)),z(eig(2)),z(eig(3)),
c     2                         z(vecl(1)),z(vecl(2)),z(vecl(3)),
c     3                         z(vecr(1)),z(vecr(2)),z(vecr(3)),     
c     4                         z(pvec),z(hpvec),z(vec),z(rhs),
c     5                         z(h),z(htmp),z(b),z(btmp),
c     6                         z(vec),z(work(3)),z(work(4)),
c     7                         scale,cnverg,thresh,
c     8                         nmax(ug,1),nmax(ug,2),nmax(ug,3),
c     9                         numdim,n3d,nrhs,niter,nvec,
c     x                         lenbuf,hamd,incore,nel,prlin,mattyp)
                  title='solution vectors energy = '
     1                    //fptoc(energy(ien))
                  call prntcm(title,z(rhs),n3d,nrhs,n3d,nrhs,iout)
               else
                  call lindvd(z(hbuf),ia(ibuf),z(diag),energy(ien),
     1                        z(eig(1)),z(eig(2)),z(eig(3)),
     2                        z(vecl(1)),z(vecl(2)),z(vecl(3)),
     3                        z(pvec),z(hpvec),z(vec),z(rhs),
     4                        z(h),z(htmp),z(b),z(btmp),
     5                        z(vec),z(work(3)),z(work(4)),
     6                        scale,cnverg,thresh,
     7                        nmax(ug,1),nmax(ug,2),nmax(ug,3),
     8                        numdim,n3d,nrhs,niter,nvec,
     9                        lenbuf,hamd,incore,nel,prlin,mattyp)
                  title='solution vectors energy = '
     1                    //fptoc(energy(ien))
                  call prntrm(title,z(rhs),n3d,nrhs,n3d,nrhs,iout)     
               endif               
 300        continue
         endif
      else    
        if(numdim.eq.1) then
           do 400 ien=1,nen
              call lschr(z(ham(1)),z(ham(1)),
     1                   dum,cdum,
     2                   dum,cdum,
     3                   z(vtot),z(hamtot),z(hamtot),
     4                   z(rhs),z(rhs),energy(ien),ia(ind),ia(ipvt),
     5                   n3d,nrhs,nmax(ug,1),nmax(ug,2),nmax(ug,3),
     6                   numdim,.true.,.true.,prnht,mattyp)
              call tocord(z(rhs),z(rhs),z(pn(ug,1)),n3d,mattyp)
              if(mattyp.eq.'complex') then
                 title='solution vectors energy = '
     1                 //fptoc(energy(ien))
                 call prntcm(title,z(rhs),n3d,nrhs,n3d,nrhs,iout)
                 title='exact solution vectors energy = '
     1                 //fptoc(energy(ien))
                 call exact(z(rhs),z(rhs),z(eigc(ug,1)),n3d,
     1                      mattyp,'delta') 
                 call prntcm(title,z(rhs),n3d,nrhs,n3d,nrhs,iout)
              else
                 title='solution vectors energy = '
     1                 //fptoc(energy(ien))
                 call prntrm(title,z(rhs),n3d,nrhs,n3d,nrhs,iout)
                 title='exact solution vectors energy = '
     1                 //fptoc(energy(ien))
                 call exact(z(rhs),z(rhs),z(eigc(ug,1)),n3d,
     1                      mattyp,'delta') 
                 call prntrm(title,z(rhs),n3d,nrhs,n3d,nrhs,iout)
              endif
 400       continue     
         elseif(numdim.eq.2) then
            do 500 ien=1,nen
               call lschr(z(ham(1)),z(ham(1)),
     1                    z(ham(2)),z(ham(2)),
     2                    dum,cdum,
     3                    z(vtot),z(hamtot),z(hamtot),
     4                    z(rhs),z(rhs),energy(ien),ia(ind),ia(ipvt),
     5                    n3d,nrhs,nmax(ug,1),nmax(ug,2),nmax(ug,3),
     6                    numdim,.true.,.true.,prnht,mattyp) 
 500        continue                         
         elseif(numdim.eq.3) then
            do 600 ien=1,nen
               call lschr(z(ham(1)),z(ham(1)),
     1                    z(ham(2)),z(ham(2)),
     2                    z(ham(3)),z(ham(3)),
     3                    z(vtot),z(hamtot),z(hamtot),
     4                    z(rhs),z(rhs),energy(ien),ia(ind),ia(ipvt),
     5                    n3d,nrhs,nmax(ug,1),nmax(ug,2),nmax(ug,3),
     6                    numdim,.true.,.true.,prnht,mattyp)
 600        continue              
         endif
         if(scatt) then
c            do 700 i=1,nen
c               call rmtrx(z(pn(ug,1)),z(eigc(ug,1)),z(eigtot),
c     1                    z(bigvec),energy(i),typot,n3d)
c  700       continue         
         endif
      endif
      call iosys ('write integer maxsiz to rwf',1,canget,0,' ')
      call chainx(0)               
      stop
 1    format(/,15x,'coordinate system = ',a16)
 2    format(/,15x,'program options',/,/,5x,
     1             'diagonalize zeroth-order hamiltonian')  
 3    format(/,15x,'code options',/,
     1       /,5x,'print q1 polynomials                   = ',l1,
     2       /,5x,'print q1 points/weights                = ',l1,
     3       /,5x,'print spatial hamiltonian information  = ',l1, 
     4       /,5x,'print q1 eigenvectors                  = ',l1, 
     5       /,5x,'check q1 orthonormality                = ',l1)
 4    format(/,15x,'code options',/,
     1       /,5x,'print q1 polynomials                   = ',l1,
     2       /,5x,'print q2 polynomials                   = ',l1,
     3       /,5x,'print q1 points/weights                = ',l1,
     4       /,5x,'print q2 points/weights                = ',l1,
     5       /,5x,'print spatial hamiltonian information  = ',l1, 
     6       /,5x,'print q1 eigenvectors                  = ',l1, 
     7       /,5x,'print q2 eigenvectors                  = ',l1, 
     8       /,5x,'check q1 orthonormality                = ',l1, 
     9       /,5x,'check q2 orthonormality                = ',l1) 
 5    format(/,15x,'code options',/,
     1       /,5x,'print q1 polynomials                   = ',l1,
     2       /,5x,'print q2 polynomials                   = ',l1,
     3       /,5x,'print q3 polynomials                   = ',l1,
     4       /,5x,'print q1 points/weights                = ',l1,
     5       /,5x,'print q2 points/weights                = ',l1,
     6       /,5x,'print q3 points/weights                = ',l1,
     7       /,5x,'print spatial hamiltonian information  = ',l1, 
     8       /,5x,'print q1 eigenvectors                  = ',l1, 
     9       /,5x,'print q2 eigenvectors                  = ',l1, 
     x       /,5x,'print q3 eigenvectors                  = ',l1, 
     x       /,5x,'check q1 orthonormality                = ',l1, 
     x       /,5x,'check q2 orthonormality                = ',l1, 
     x       /,5x,'check q3 orthonormality                = ',l1) 
 6    format(/,15x,'trap configuration and atomic parameters',/,/,5x,
     1             'atomic mass                      = ',e15.8,/,5x,
     2             'trap q1 frequency                = ',e15.8)
 7    format(/,15x,'trap configuration and atomic parameters',/,/,5x,
     1             'atomic mass                       = ',e15.8,/,5x,
     2             'trap q1 frequency                 = ',e15.8,/,5x,
     3             'trap q2 frequency                 = ',e15.8)
 8    format(/,15x,'trap configuration and atomic parameters',/,/,5x,
     1             'atomic mass                       = ',e15.8,/,5x,
     2             'trap q1 frequency                 = ',e15.8,/,5x,
     3             'trap q2 frequency                 = ',e15.8,/,5x,
     4             'trap q3 frequency                 = ',e15.8)
 9    format(/,15x,'iterative linear system information',/,/,5x,
     1             'number of right hand sides         = ',i3,/,5x,
     2             'overlap tolerance                  = ',e15.8,/,5x,
     3             'convergence criterion              = ',e15.8,/,5x,
     4             'maximum number of iterations       = ',i6,/,5x,
     5             'maximum number of vectors          = ',i6)
 11   format(/,15x,'direct solve',/,/,5x,
     1             'size of matrix             = ',i3,/,5x,
     2             'number of right hand sides = ',i3)
 12   format(/,5x,'dimension = ',i2,2x,'grid = ',i3)
      end



