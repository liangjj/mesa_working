*deck colloc.f
c***begin prologue     colloc
c***date written       950704   (yymmdd)
c***revision date               (yymmdd)
c***keywords           collocation, scattering
c***author             schneider, barry (nsf)
c***source             mesa
c***purpose            solve one-dimensional bound or scattering equations
c***                   by combination of collocation and finite difference.
c***
c***description        one-dimensional schroedinger equation is solved using
c***                   a combination of collocation and finite difference.
c***                   the equation is treated in spherical co-ordinates.
c***references         
c
c***routines called 
c***                   
      program colloc
      implicit integer(a-z)
      parameter ( ngmax=10, nrmax=100 )
      common ia(1)
      real*8 z(1), fpkey, energy, rbndry
      real*8 rbeg, rend
      dimension nreg(3), npts(nrmax,ngmax,3), nwts(nrmax,ngmax,3)
      dimension rbeg(nrmax,3), rend(nrmax,3), qtype(3), rbndry(2,3)
      dimension ncolfn(ngmax), ptot(ngmax,3), wtot(ngmax,3)
      dimension ptotu(ngmax,3)
      character*4096 ops
      character*128 fillam
      character*16 pottyp, type, qtype
      character*24 cpass
      character*3 itoc
      character*1600 card
      character*80  chrkey, title, dirct, inv
      logical logkey, prntbf, prntiv, prntfv, prntdr, prntsl
      logical prntih, prntv, prntev, prntg, prnth, fixed, nobc
      logical zroinh, posinp
      equivalence (ia(1),z(1))
      common /memory/ ioff
      common/io/inp,iout
c
      call drum
      call iosys ('read character options from rwf',-1,0,0,ops)
      call iosys ('read character "linear algebraic filename" from rwf',
     1                                                  -1,0,0,fillam)
      call iosys ('open lamdat as unknown',0,0,0,fillam)
      write(iout,*)
      title='solve simple schroedinger equation via collocation/finite d 
     1ifference'
      write(iout,*) title
      write(iout,*)
      type=chrkey(ops,'type-fitting-function','sine',' ')
      prntbf=logkey(ops,'print=m6232=basis-information',.false.,' ')
      prntiv=logkey(ops,'print=m6232=inverse',.false.,' ')
      prntfv=logkey(ops,'print=m6232=fitting-vector',.false.,' ')
      prntdr=logkey(ops,'print=m6232=derivative-matrix',.false.,' ')
      prntv=logkey(ops,'print=m6232=potential',.false.,' ')
      prntg=logkey(ops,'print=m6232=grid',.false.,' ')
      prntev=logkey(ops,'print=m6232=effective-potential',.false.,' ')
      prntih=logkey(ops,'print=m6232=inhomogeneity',.false.,' ')
      prntsl=logkey(ops,'print=m6232=solution',.false.,' ')
      prnth=logkey(ops,'print=m6232=hamiltonian',.false.,' ')
      nobc=logkey(ops,'no-boundary-condition',.false.,' ')
      zroinh=logkey(ops,'zero-inhomogeneity',.false.,' ')
      l=intkey(ops,'angular-momentum',0,' ')
      pottyp=chrkey(ops,'potential-type','exponential',' ')
      inv=chrkey(ops,'inversion-type','pseudo',' ')
      order=intkey(ops,'order-numerov',3,' ')
      ng=intkey(ops,'number-of-grids',1,' ')
      frezc=intkey(ops,'freeze-collocation-grid',0,' ')
      nener=intkey(ops,'number-of-energies',1,' ')
c     basic coarse grid information and region definition
      ptcnt=0
      wtcnt=0
      if (posinp('$collocation',cpass) ) then
          call cardin(card)
          call dtacol(card,qtype(1),rbeg(1,1),rend(1,1),nreg(1),
     1                npts(1,1,1),nwts(1,1,1),ptot(1,1),wtot(1,1),
     2                ptcl,wtcl,ncorfn,ng,nrmax,ngmax)
          ncolfn(1)=ncorfn
          ptcnt=ptcnt+ptcl
          wtcnt=wtcnt+wtcl
          rbndry(1,1)=rbeg(1,1)
          rbndry(2,1)=rend(nreg(1),1)
      endif
      if (posinp('$finite-difference',cpass) ) then
          call cardin(card)
          call dtafd(card,qtype(2),rbeg(1,2),rend(1,2),nreg(2),
     1               npts(1,1,2),nwts(1,1,2),ptot(1,2),wtot(1,2),
     2               ptfd,wtfd,ng,nrmax,ngmax)
          fixed=.false.
          if(qtype(2).eq.'default') then
             fixed=.true.
          endif
          ptcnt=ptcnt+ptfd
          wtcnt=wtcnt+wtfd
          rbndry(1,2)=rbeg(1,2)
          rbndry(2,2)=rend(nreg(2),2)
      endif
      if (posinp('$fitting',cpass) ) then
          call cardin(card)
          call dtaft(card,qtype(3),rbeg(1,3),rend(1,3),nreg(3),
     1               npts(1,1,3),nwts(1,1,3),ptot(1,3),wtot(1,3),
     2               ptft,wtft,ng,nrmax,ngmax)
          ptcnt=ptcnt+ptft
          wtcnt=wtcnt+wtft
          rbndry(1,3)=rbeg(1,3)
          rbndry(2,3)=rend(nreg(3),3)
      endif
      bw=1
      if (order.eq.5) then
          bw=2
      endif
      if (ptot(1,2).eq.0.and.ptot(1,1).eq.0) then
          call lnkerr('error in point allocation')
      endif
      if (ptot(1,2).eq.0) then
          dirct='collocation' 
      elseif(ptot(1,1).eq.0) then
          dirct='numerov'
      else
          dirct='collocation-and-numerov'
      endif          
      write(iout,*)
      write(iout,*) '  type calculation = ',dirct
      write(iout,*) title                    
      write(iout,1) ng, ptot(1,1), ptot(1,2), ptot(1,3)
      write(iout,2) ncorfn, nener, pottyp, l
      if (ptot(1,1).ne.0) then
          write(iout,3) rbndry(1,1), rbndry(2,1)
      endif
      if (ptot(1,2).ne.0) then
          write(iout,4) rbndry(1,2), rbndry(2,2)
      endif
      if (ptot(1,3).ne.0) then
          write(iout,5) rbndry(1,3), rbndry(2,3)
      endif
c     generate basic grid information
      maxcol=ptot(ng,1)      
      maxnum=ptot(ng,2)
      maxfn=ncolfn(1)
      nptsq=maxcol*maxcol
      ndim=maxcol*maxfn
      ioff=1
      do 50 i=1,2
c        storage for all grids, potentials, effective potentials
c        and inhomogeneities are outlayed
         rpt=ioff
         wt=rpt+ptcnt
         runiqe=wt+wtcnt 
         v=runiqe+ptcnt
         vff=v+ptcnt
         g=vff+ptcnt
         rhs=g+ptcnt
         wdsusd=rhs+ptcnt
         if (ptot(1,1).ne.0.and.ptot(1,2).eq.0) then
c            all of the functions, their first and second derivatives
c            are outlayed.
             fns=wdsusd
             dfns=fns+ndim
             ddfns=dfns+ndim
             dmat=ddfns+ndim
             ddmat=dmat+nptsq
             hpp=ddmat+nptsq
             ainv=hpp+maxcol*maxcol
             workp=ainv+maxcol*maxfn
             ipvtp=wpadti(workp+maxcol)
c            the inverse matrices and scratch arrays are outlayed
c            for the biggest grid.
             t=iadtwp(ipvtp+maxcol)
             e=t+maxcol*maxfn
             wdsusd=e+maxfn
             if(inv.ne.'inverse') then
                s=e+maxfn
                u=s+maxfn+1
                vv=u+maxcol*maxfn
                wdsusd=vv+maxfn*maxfn
             endif                
         endif
         if (ptot(1,2).ne.0.and.ptot(1,1).eq.0) then
             bandy=wdsusd
             bandg=bandy+order*maxnum
             workq=bandg+order*maxnum
             ipvtq=wpadti(workq+maxnum)
             wdsusd=iadtwp(ipvtq+maxnum)
         endif
         if (ptot(1,1).ne.0.and.ptot(1,2).ne.0) then
             fns=wdsusd
             dfns=fns+ndim
             ddfns=dfns+ndim
             dmat=ddfns+ndim
             ddmat=dmat+nptsq
             ham=ddmat+nptsq
             ainv=ham+(maxcol+maxnum)*(maxcol+maxnum)
c            the inverse matrices and scratch arrays are outlayed
c            for the biggest grid.
             t=ainv+maxcol*maxfn
             workp=t+maxcol*maxfn
             e=workp+maxcol
             wdsusd=e+maxfn
             if(inv.ne.'inverse') then
                s=e+maxfn
                u=s+maxfn+1
                vv=u+maxcol*maxfn
                ipvtp=wpadti(vv+maxfn*maxfn)
                wdsusd=iadtwp(ipvtp+maxcol)
             endif                
             bandy=ainv
             bandg=bandy+order*maxnum
             work=bandg+order*maxnum
             workp=work
             ipvt=wpadti(work+maxcol+maxnum)
             ipvtp=ipvt        
             wdsusd=iadtwp(ipvt+maxcol+maxnum)   
         endif
         if (i.eq.1) then
             need=wpadti(wdsusd)
             call getscm(need,z,maxcor,'colloc',0)
         endif
 50   continue
      ir=rpt
      iru=runiqe
      iwt=wt
      iv=v
      do 60 i=1,ng
c        if there is more than one region, the last point in region (i-1)
c        is the first point in region i.
c        make the point and potential arrays         
         if (ptot(1,1).ne.0) then
             write(iout,*) 
             title='generating grid and potential for collocation'
             len=length(title)
             title=title(1:len)//' region'
             write(iout,*) title
             call mkgrd(qtype(1),z(ir),z(iru),z(iwt),rbeg(1,1),
     1                  rend(1,1),nreg(1),npts(1,i,1),ptotu(i,1),
     2                  nrmax,ngmax,prntg)
             call potntl(z(iv),z(iru),ptotu(i,1),pottyp,prntv)
             ir=ir+ptot(i,1)
             iwt=iwt+wtot(i,1)
             iv=iv+ptot(i,1)
             iru=iru+ptotu(i,1)
             rlst=z(ir-1)
         endif
         if(ptot(1,2).ne.0) then
             write(iout,*) 
             title='generating grid and potential for finite'
             len=length(title)
             title=title(1:len)//' difference region'
             write(iout,*) title
             call mkgrd(qtype(2),z(ir),z(iru),z(iwt),rbeg(1,2),
     1                  rend(1,2),nreg(2),npts(1,i,2),ptotu(i,2),
     2                  nrmax,ngmax,prntg)
             call potntl(z(iv),z(iru),ptotu(i,2),pottyp,prntv)
             ir=ir+ptot(i,2)
             iru=iru+ptotu(i,2)
             iwt=iwt+wtot(i,2)
             iv=iv+ptot(i,2)
             rlst=z(ir-1)
         endif
         if(ptot(1,3).ne.0) then
             write(iout,*) 
             title='generating grid and potential for fitting region'
             write(iout,*) title
             call mkgrd(qtype(3),z(ir),z(iru),z(iwt),rbeg(1,3),
     1                  rend(1,3),nreg(3),npts(1,i,3),ptotu(i,3),
     2                  nrmax,ngmax,prntg)
             call potntl(z(iv),z(iru),ptotu(i,3),pottyp,prntv)
             ir=ir+ptot(i,3)
             iru=iru+ptotu(i,3)
             iwt=iwt+wtot(i,3)
             iv=iv+ptot(i,3)
             rlst=z(ir-1)
         endif
c     add a few points at the end for formulas requiring points outside their
c     range.
      title='generating end points and zeroing potential'
      write(iout,*) title
c      rend=rlst+stplst+stplst
c      call mkgrd(z(ir),rlst,rend,stplst,3,prntg)
c      call potntl(z(iv),z(ir),3,pottyp,prntv)
c      call rzero(z(iv),3)
c      ir=ir+3
c      iv=iv+3         
 60   continue   
      if (ptot(1,1).ne.0) then
          call posinp('$funct',cpass)
          ir=rpt
          iru=runiqe
          if=fns
          iv=v
          idf=dfns
          iddf=ddfns
          idmat=dmat
          iddmat=ddmat
          do 70 i=1,ng
c            calculate the collocation basis functions and/or their
c            derivatives at the collocation points.             
             call funct(z(iru),z(if),z(idf),z(iddf),ncolfn(i),
     1                  ptotu(i,1),prntbf)
             call copy(z(if),z(ainv),ptotu(i,1)*ncolfn(i))
c            make the first and second derivative matrices 
             call derivs(z(idmat),z(iddmat),z(idf),z(iddf),z(ainv),
     1                   z(t),z(s),z(e),z(u),z(vv),z(workp),
     2                   ptotu(i,1),ncolfn(i),prntiv,prntdr,inv)
             ir=ir+ptot(i,1)
             iru=iru+ptotu(i,1)
             iv=iv+ptot(i,1)
             if=if+ncolfn(i)*ptot(i,1)
             idf=idf+ncolfn(i)*ptot(i,1)
             iddf=iddf+ncolfn(i)*ptot(i,1)
             idmat=idmat+ptot(i,1)*ptot(i,1)
             iddmat=iddmat+ptot(i,1)*ptot(i,1)
 70       continue
      endif
c      call lnkerr('quit')
      do 80 ene=1,nener
c
c        energy is input in Rydbergs
c      
         call posinp('$energy-'//itoc(ene),cpass)
         call cardin(card)         
         energy=fpkey(card,'energy',1.d0,' ')
         write (iout,1000) energy
         ir=rpt
         iru=runiqe
         iv=v
         ivff=vff
         ig=g
         do 90 i=1,ng
c        compute the effective potential and then make the right hand side 
c        array.  the latter is needed because we remove the inhomogeneous term 
c        from the solution and bring it to the right hand side in order to
c        impose "outgoing" wave boundary conditions.   
            if (ptot(1,1).ne.0) then
                call veff(z(iv),z(ivff),z(iru),energy,l,ptotu(i,1),
     1                    prntev)
                call mkg(z(ig),z(iv),z(iru),energy,l,ptotu(i,1),
     1                   zroinh,prntg)
                ir=ir+ptot(i,1)
                iru=iru+ptotu(i,1)
                iv=iv+ptot(i,1)
                ivff=ivff+ptot(i,1)
                ig=ig+ptot(i,1)
            endif
            if(ptot(1,2).ne.0) then
                call veff(z(iv),z(ivff),z(iru),energy,l,ptotu(i,2),
     1                    prntev)
                call mkg(z(ig),z(iv),z(iru),energy,l,ptotu(i,2),zroinh,
     1                   prntg)
                ir=ir+ptot(i,2)
                iru=iru+ptotu(i,2)
                iv=iv+ptot(i,2)
                ivff=ivff+ptot(i,2)       
                ig=ig+ptot(i,2)
            endif
            if(ptot(1,3).ne.0) then
                call veff(z(iv),z(ivff),z(iru),energy,l,ptotu(i,3),
     1                    prntev)
                call mkg(z(ig),z(iv),z(iru),energy,l,ptotu(i,3),
     1                   zroinh,prntg)
                ir=ir+ptot(i,3)
                iru=iru+ptotu(i,3)
                iv=iv+ptot(i,3)
                ivff=ivff+ptot(i,3)       
                ig=ig+ptot(i,3)
            endif
c           add the end points            
c            call veff(z(iv),z(ivff),z(ir),energy,l,3,prntev)
c            call mkg(z(ig),z(iv),z(ir),energy,l,3,zroinh,prntg)
c            ir=ir+3
c            iv=iv+3
c            ivff=ivff+3      
c            ig=ig+3
 90      continue
c----------------------------------------------------------------------c
c     the structure of the matrix is very interesting.  if we consider
c     the collocation block as a p-space and the numerov block as q-space
c     the hpp block will in general be a full matrix, the hpq block will
c     be zero, the hqp block will have only one non-zero entry, which will
c     be in the first row of that block as the last element and the hqq
c     block will be tridiagonal.  this is tailor made for special matrix
c     methods and we will exploit that later.
c----------------------------------------------------------------------c
      if(ptot(1,1).ne.0.and.ptot(1,2).eq.0) then
c     set up collocation block of matrix if required
         np=ptotu(1,1)
         call matpp(z(hpp),z(dmat),z(ddmat),z(vff),z(g),z(rhs),
     1              z(workp),energy,z(runiqe),rbndry(2,1),ia(ipvtp),l,
     2              np,np,prnth,prntsl,dirct,nobc)
      elseif(ptot(1,1).eq.0.and.ptot(1,2).ne.0) then
             nq=ptotu(1,2)
             call numerv(z(vff),z(g),z(rhs),z(bandy),z(bandg),z(workq),
     1                   energy,z(rpt),rbndry(2,2),ia(ipvtq),l,
     2                   nq,order,bw,fixed,prnth,prntsl)
      elseif(ptot(1,1).ne.0.and.ptot(1,2).ne.0) then
         np=ptotu(1,1)
         nq=ptotu(1,2)-1
         ntot=np+nq
         call rzero(z(ham),ntot*ntot)
         call matpp(z(ham),z(dmat),z(ddmat),z(vff),z(g),z(rhs),
     1              z(work),energy,z(runiqe),rbndry(2,1),ia(ipvt),l,
     2              np,ntot,prnth,prntsl,dirct)
c         to complete the problem it is necessary to calculate the
c         numerov portion of the discretization, the coupling of
c         collocation to the numerov block and the coupling of the numerov
c         to the collocation block.  the coupling of the collocation to 
c         the numerov block is the result of replacing the last collocation
c         equation by one which insures that the derivative of the 
c         wavefunction is continuous at the interface point.
          hpq=ham+np*ntot  
          call matpq(z(ham),z(hpq),z(dmat),z(rpt+np),z(vff),
     1               z(g),z(rhs),np,nq,ntot,prnth)
c         solve formally for the p-space part of the wavefunction
          call solvep(z(ham),z(hpq),z(rhs),ia(ipvtp),
     1                np,nq,ntot,prntsl)
          hqp=ham+np
          call matqp(z(hqp),z(rpt+np),z(vff+np),z(g+np),
     1               z(rhs+np),3,np,nq,ntot,prnth)
          hqq=hpq+np
          call matqq(z(hqq),z(rpt+np),z(vff+np),z(g+np),z(rhs+np),
     1               z(bandy),z(bandg),z(work),energy,rbndry(2,2),
     2               l,3,np,nq,ntot,order,bw,fixed,prnth)
          call newqq(z(hqq),z(hpq),z(rhs),z(hqp),z(rhs+np),z(bandy),
     1               np,nq,ntot,bw,prnth) 
          call matslv(z(bandy),z(rhs+np),z(work),energy,nq,order,
     1                bw,prnth,prntsl)
      else     
             call lnkerr('error')
      endif
 80   continue
c
      call chainx(0)
c
 1    format(/,1x,'number of grids = ',i2,1x,
     1            'number of coarse grid collocation points = ',
     2                                                          i3,/,1x,
     3            'number of coarse grid finite difference points = '
     4                                                         ,i4,/,1x,
     5            'number of coarse grid fitting  points          = '
     6                                                              ,i3)
 2    format(/,1x,'number of coarse grid collocation functions = ',i3,
     1                                                             /,1x,
     2            'number of energies                          = ',i3,
     3                                                             /,1x,
     4            'potential type                              = ',a16,
     5                                                            /,1x,
     6            'angular momentum                            = ',i2)
 3    format(/,1x,'collocation region begins at = ',e15.8,/,1x,
     1            'collocation region ends at   = ',e15.8)
 4    format(/,1x,'finite difference region begins at = ',e15.8,/,1x,
     1            'finite difference region ends at   = ',e15.8)
 5    format(/,1x,'fitting difference region begins at = ',e15.8,/,1x,
     1            'fitting difference region ends at   = ',e15.8)

 9    format(7x,i4,10x,i4,19x,f10.5)
 1000 format(/,1x,'calculation at energy = ',e15.8)
      stop
      end

