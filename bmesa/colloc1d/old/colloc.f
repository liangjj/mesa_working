*deck colloc.f
c***begin prologue     collloc
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
      parameter ( ngrids=10 )
      common ia(1)
      real*8 z(1), fpkey,energy, rbndry
      real*8 stp, delr
      dimension stp(ngrids,3), npts(ngrids,3), delr(3), rbndry(2,3)
      dimension ncolfn(ngrids), nsum(3)
      character*4096 ops
      character*128 fillam
      character*16 cpass, pottyp, type
      character*3 itoc
      character*800 card
      character*80  chrkey, title, dirct
      logical logkey, prntbf, prntiv, prntfv, prntdr, prntsl
      logical prntih, prntv, prntev, prntg, prnth
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
      call posinp('$colloc',cpass)
      call cardin(card)
c     basic information on angular momentum, potential and energies
      l=intkey(ops,'angular-momentum',0,' ')
      pottyp=chrkey(ops,'potential-type','exponential',' ')
      nener=intkey(card,'number-of-energies',1,' ')
c     basic coarse grid information and region definition
      ng=intkey(card,'number-of-grids',1,' ')
      ncorcl=intkey(card,'number-of-coarse-grid-collocation-points',
     1                    3,' ')
c     make ncorcl odd and use a minimum of three points
      if (ncorcl.ne.0) then
          frezc=intkey(card,'freeze-collocation-grid',0,' ')
          call fparr(card,'collocation-region-boundaries',
     1               rbndry(1,1),2, ' ')
          ntst=mod(ncorcl,2)
          if(ntst.eq.0) then
             ncorcl=ncorcl+1
          endif
          ncorcl=max(ncorcl,3)
          ncorfn=intkey(card,'number-of-coarse-grid-'//
     1                       'collocation-functions',ncorcl, ' ')
          ncolfn(1)=ncorfn
          npts(1,1)=ncorcl
          nsum(1)=npts(1,1)
          delr(1) =  rbndry(2,1)-rbndry(1,1)
          stp(1,1) = delr(1)/( npts(1,1)-1 )
      endif
      ncorfd=intkey(card,'number-of-coarse-grid-finite-'//
     1                   'difference-points',26,' ')
      order=0
      if (ncorfd.ne.0) then     
          call fparr(card,'finite-difference-region-boundaries',
     1               rbndry(1,2),2,' ')
c         make ncorfd odd and use a minimum of three points
          ntst=mod(ncorfd,2)
          if(ntst.eq.0) then
             ncorfd=ncorfd+1
          endif
          ncorfd=max(ncorfd,3)
          npts(1,2)=ncorfd
          nsum(2)=npts(1,2)
          delr(2) = rbndry(2,2)-rbndry(1,2)
          stp(1,2) = delr(2)/( npts(1,2)-1 )
          order=intkey(ops,'order-numerov',3,' ')
          bw=1
          if (order.eq.5) then
              bw=2
          endif
      endif
      ncorft=intkey(card,'number-of-coarse-grid-fitting-points',
     1                    3,' ')
      if (ncorft.ne.0) then     
          call fparr(card,'fitting-region-boundaries',rbndry(1,3),
     1               2,' ')
c         use a minimum of two points
          if (ncorft.lt.2) then
              ncorft=2
          endif
          npts(1,3)=ncorft
          nsum(3)=npts(1,3)
          delr(3) = rbndry(2,3)-rbndry(1,3)
          stp(1,3) = delr(3)/( npts(1,3)-1 )
      endif
      type=chrkey(ops,'type-fitting-function','sine',' ')
      prntbf=logkey(ops,'print=m1204=basis-information',.false.,' ')
      prntiv=logkey(ops,'print=m1204=inverse',.false.,' ')
      prntfv=logkey(ops,'print=m1204=fitting-vector',.false.,' ')
      prntdr=logkey(ops,'print=m1204=derivative-matrix',.false.,' ')
      prntv=logkey(ops,'print=m1204=potential',.false.,' ')
      prntg=logkey(ops,'print=m1204=grid',.false.,' ')
      prntev=logkey(ops,'print=m1204=effective-potential',.false.,' ')
      prntih=logkey(ops,'print=m1204=inhomogeneity',.false.,' ')
      prntsl=logkey(ops,'print=m1204=solution',.false.,' ')
      prnth=logkey(ops,'print=m1204=hamiltonian',.false.,' ')
      if (ncorfd.eq.0.and.ncorcl.eq.0) then
          call lnkerr('error in point allocation')
      endif
      if (ncorfd.eq.0) then
          dirct='collocation' 
      elseif(ncorcl.eq.0) then
          dirct='numerov'
      else
          dirct='collocation-and-numerov'
      endif          
      write(iout,*)
      write(iout,*) '  type calculation = ',dirct
      write(iout,*) title                    
      write(iout,1) ng, ncorcl, ncorfd, ncorft
      write(iout,2) ncorfn, nener, pottyp, l
      if (ncorfd.ne.0) then
          write(iout,3) rbndry(1,1), rbndry(2,1)
      endif
      if (ncorfd.ne.0) then
          write(iout,4) rbndry(1,2), rbndry(2,2)
      endif
      if (ncorft.ne.0) then
          write(iout,5) rbndry(1,3), rbndry(2,3)
      endif
c     generate basic grid information
      maxcol=npts(1,1)      
      maxnum=npts(1,2)
      maxfn=ncolfn(1)
      nptsq=npts(1,1)*npts(1,1)
      ndim=npts(1,1)*ncolfn(1)
      do 10 i=2,ng
         if (ncorcl.ne.0) then
             npts(i,1) = npts(i-1,1) + npts(i-1,1) -1
             stp(i,1) = .5d0*stp(i-1,1)
             ncolfn(i)=npts(i,1)
             nsum(1)=nsum(1)+npts(i,1)
             ndim=ndim+npts(i,1)*ncolfn(i)
             nptsq=nptsq+npts(i,1)*npts(i,1)
             maxcol=max(maxcol,npts(i,1))
             maxfn=max(maxfn,ncolfn(i))
         endif
         if (ncorfd.ne.0) then
             npts(i,2) = npts(i-1,2) + npts(i-1,2) -1
             stp(i,2) = .5d0*stp(i-1,2)
             nsum(2)=nsum(2)+npts(i,2)
             maxnum=max(maxnum,npts(i,2))
         endif
         if (ncorft.ne.0) then
             npts(i,3) = npts(i-1,3) + npts(i-1,3) -1
             stp(i,3) = .5d0*stp(i-1,3)
             nsum(3)=nsum(3)+npts(i,3)
         endif
 10   continue
      if (ncorcl.ne.0) then
          write(iout,6)
          do 20 i=1,ng
             write(iout,9) i, npts(i,1), stp(i,1)
 20       continue
      endif   
      if (ncorfd.ne.0) then
          write(iout,7)
          do 30 i=1,ng
             write(iout,9) i, npts(i,2), stp(i,2)
 30       continue
      endif   
      if (ncorft.ne.0) then
          write(iout,8)
          do 40 i=1,ng
             write(iout,9) i, npts(i,3), stp(i,3)
 40       continue
      endif   
      nsumt=nsum(1)+nsum(2)+nsum(3)+10*ng
      ioff=1
      do 50 i=1,2
c        storage for all grids, potentials, effective potentials
c        and inhomogeneities are outlayed
         rpt=ioff
         v=rpt+nsumt
         vff=v+nsumt
         g=vff+nsumt
         rhs=g+nsumt
         wdsusd=rhs+nsumt
         if (ncorcl.ne.0) then
c            all of the functions, their first and second derivatives
c            are outlayed.
             fns=wdsusd
             dfns=fns+ndim
             ddfns=dfns+ndim
             dmat=ddfns+ndim
             ddmat=dmat+nptsq
             hampp=ddmat+nptsq
             ainv=hampp+maxcol*maxcol
             wdsusd=ainv
c            the inverse matrices, work and ipvt arrays are outlayed
c            for the biggest grid.
             workp=ainv+maxcol*maxcol
             ipvtp=wpadti(workp+maxcol+maxcol)           
         endif
         if (ncorfd.ne.0) then
             band=wdsusd
             workq=band+order*maxnum
             ipvtq=wpadti(workq+maxnum)
             ham=iadtwp(ipvtq+maxnum)
             wdsusd=ham+(maxnum+maxcol)*(maxnum+maxcol)
         endif
         if (i.eq.1) then
             need=wpadti(wdsusd)
             call getscm(need,z,maxcor,'colloc',0)
         endif
 50   continue
      ir=rpt
      iv=v
      do 60 i=1,ng
c        if there is more than one region, the last point in region (i-1)
c        is the first point in region i.
c        make the point and potential arrays         
         if (ncorcl.ne.0) then
             title='generating grid and potential for collocation'
             len=length(title)
             title=title(1:len)//' region'
             write(iout,*) title
             call mkgrd(z(ir),rbndry(1,1),rbndry(2,1),
     1                  stp(i,1),npts(i,1),0,prntg)
             call potntl(z(iv),z(ir),npts(i,1),0,pottyp,prntv)
             ir=ir+npts(i,1)
             iv=iv+npts(i,1)
         endif
         if(ncorfd.ne.0) then
             title='generating grid and potential for finite'
             len=length(title)
             title=title(1:len)//' difference region'
             write(iout,*) title
             call mkgrd(z(ir),rbndry(1,2),rbndry(2,2),
     1                  stp(i,2),npts(i,2),0,prntg)
             call potntl(z(iv),z(ir),npts(i,2),0,pottyp,prntv)
             ir=ir+npts(i,2)
             iv=iv+npts(i,2)
         endif
         if(ncorft.ne.0) then
             title='generating grid and potential for fitting region'
             write(iout,*) title
             call mkgrd(z(ir),rbndry(1,3),rbndry(2,3),
     1                  stp(i,3),npts(i,3),0,prntg)
             call potntl(z(iv),z(ir),npts(i,3),0,pottyp,prntv)
             ir=ir+npts(i,3)
             iv=iv+npts(i,3)
         endif
 60   continue   
      if (ncorcl.ne.0) then
          call posinp('$funct',cpass)
          ir=rpt
          if=fns
          iv=v
          idf=dfns
          iddf=ddfns
          idmat=dmat
          iddmat=ddmat
          do 70 i=1,ng
c            calculate the collocation basis functions and/or their
c            derivatives at the collocation points.             
             call funct(z(ir),z(if),z(idf),z(iddf),ncolfn(i),
     1                  npts(i,1),prntbf)
             call copy(z(if),z(ainv),npts(i,1)*npts(i,1))
c            make the first and second derivative matrices 
             call derivs(z(idmat),z(iddmat),z(idf),z(iddf),z(ainv),
     1                   z(workp),ia(ipvtp),npts(i,1),ncolfn(i),
     2                   prntiv,prntdr)
             ir=ir+npts(i,1)
             iv=iv+npts(i,1)
             if=if+ncolfn(i)*npts(i,1)
             idf=idf+ncolfn(i)*npts(i,1)
             iddf=iddf+ncolfn(i)*npts(i,1)
             idmat=idmat+npts(i,1)*npts(i,1)
             iddmat=iddmat+npts(i,1)*npts(i,1)
 70       continue
      endif
      do 80 ene=1,nener
c
c        energy is input in Rydbergs
c      
         call posinp('$energy-'//itoc(ene),cpass)
         call cardin(card)         
         energy=fpkey(card,'energy',1.d0,' ')
         write (iout,1000) energy
         ir=rpt
         iv=v
         ivff=vff
         ig=g
         do 90 i=1,ng
c        compute the effective potential and then make the right hand side 
c        array.  the latter is needed because we remove the inhomogeneous term 
c        from the solution and bring it to the right hand side in order to
c        impose "outgoing" wave boundary conditions.   
            if (ncorcl.ne.0) then
                call veff(z(iv),z(ivff),z(ir),energy,l,npts(i,1),prntev)
                call mkg(z(ig),z(iv),z(ir),energy,l,npts(i,1),prntg)
                ir=ir+npts(i,1)
                iv=iv+npts(i,1)
                ivff=ivff+npts(i,1)
                ig=ig+npts(i,1)
            endif
            if(ncorfd.ne.0) then
                call veff(z(iv),z(ivff),z(ir),energy,l,npts(i,2),prntev)
                call mkg(z(ig),z(iv),z(ir),energy,l,npts(i,2),prntg)
                ir=ir+npts(i,2)
                iv=iv+npts(i,2)
                ivff=ivff+npts(i,2)       
            endif
            if(ncorft.ne.0) then
                call veff(z(iv),z(ivff),z(ir),energy,l,npts(i,2),prntev)
                call mkg(z(ig),z(iv),z(ir),energy,l,npts(i,2),prntg)
                ir=ir+npts(i,3)
                iv=iv+npts(i,3)
                ivff=ivff+npts(i,3)       
            endif
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
      n=maxcol+maxnum
      np=npts(1,1)
      nq=npts(1,2)
      if(ncorcl.ne.0) then
c     set up collocation block of matrix if required
         call matpp(z(hampp),z(dmat),z(ddmat),z(vff),z(g),z(rhs),
     1              z(workp),energy,z(rpt),rbndry(2,1),ia(ipvtp),l,
     2              np,prnth,prntsl,dirct)
      elseif(ncorcl.eq.0.and.ncorfd.ne.0) then
             call numerv(z(vff),z(g),z(rhs),z(band),z(workq),
     1                   stp(1,2),energy,z(rpt),rbndry(2,2),
     2                   ia(ipvtq),l,nq,order,bw,prnth,prntsl)
      elseif(ncorcl.ne.0.and.ncorfd.ne.0) then
c         add to the collocation matrix the numerov discretization of the
c         rest of the problem.
             ntot=np+nq
             call filham(z(ham),z(hampp),ntot,np)
             call matqp(z(ham),z(rpt),z(vff),z(g),z(rhs),z(band),
     1                  3,np,nq,ntot)
             lochqq=ham+n*np+np
             call matqq(z(lochqq),z(vff+np),z(g+np),z(locrhq),z(band),
     1                  stp(1,2),3,np,nq,bw)
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
 6    format(/,6x,'   grid     collocation points           step size')
 7    format(/,6x,'   grid     finite difference points     step size')
 8    format(/,6x,'   grid     matching points              step size')
 9    format(7x,i4,10x,i4,19x,f10.5)
 1000 format(/,1x,'calculation at energy = ',e15.8)
      stop
      end

