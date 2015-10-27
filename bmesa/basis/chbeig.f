*deck chbeig.f
c***begin prologue     chbeig
c***date written       950704   (yymmdd)
c***revision date               (yymmdd)
c***keywords           r-matrix, diagonalization
c***author             schneider, barry (nsf)
c***source             mesa
c***purpose            calculate eigenvalues and eigenfunctions of 1-d
c***                   schroedinger equation by expansion in chebychev
c***                   polynomials.         
c***   
c***                   
c***references         
c
c***routines called 
c***end prologue       chbeig                   
      program chbeig
      implicit integer(a-z)
      common ia(1)
      real*8 z(1)
      character*4096 ops
      character*128 fillam
      character*16 pottyp, bc
      character*24 cpass
      character*1600 card
      character*80  chrkey, title
      logical logkey, prntv, prntg, prntsl, prnth
      logical posinp
      equivalence (ia(1),z(1))
      common /memory/ ioff
      common/io/inp,iout
c
      call drum
      call iosys ('read character options from rwf',-1,0,0,ops)
      call iosys ('read character "linear algebraic filename" from rwf',
     1                                                  -1,0,0,fillam)
      call iosys ('open lamdat as unknown',0,0,0,fillam)
      prntv=logkey(ops,'print=m6237=potential',.false.,' ')
      prntg=logkey(ops,'print=m6237=grid',.false.,' ')
      prntsl=logkey(ops,'print=m6237=solution',.false.,' ')
      prnth=logkey(ops,'print=m6237=hamiltonian',.false.,' ')
      prntvc=logkey(ops,'print=m6237=vectors',.false.,' ')
      bc=chrkey(ops,'boundary-condition','zero-function',' ')
      l=intkey(ops,'angular-momentum',0,' ')
      pottyp=chrkey(ops,'potential-type','exponential',' ')
      evecs=logkey(ops,'eigenvectors',.false.,' ')
      drct=chrkey(ops,'diagonalize','all',' ')
      if(drct.ne.'all') then
         evecs=.true.
      endif         
      extrap=logkey(ops,'extrapolate',.false.,' ')
c     basic coarse grid information and region definition
      if (posinp('$chbeig',cpass) ) then
          call cardin(card)
      endif
      call iosys('read integer "order of chebyschev polynomials" '//
     1           'from lamdat',1,n,0,' ')
      if(drct.ne.'all') then
         nroots=intkey(ops,'number-of-roots',5,' ')
      endif
      dim=2
      if(.not.extrap) then
         write(iout,3)
         dim=1
      endif
c     generate basic grid information
      ioff=1
      do 10 i=1,2
c        storage for all grids, potentials, effective potentials
c        and inhomogeneities are outlayed
         roots=ioff
         t=roots+n
         dt=t+n*n
         ddt=dt+n*n
         v=ddt+n*n
         if (i.eq.1) then
             need=wpadti(wdsusd)
             call getscm(need,z,maxcor,'finite-difference',0)
         endif
 10   continue
      call iosys('read real "chebyschev roots" from lamdat',n,
     1            z(roots),0,' ')
      call iosys('read real "chebychev polynomials" from lamdat',
     1            n*n,z(t),0,' ')
      call iosys('read real "first derivative of chebychev '//
     1           'polynomials" from lamdat',n*n,z(dt),0,' ')
      call iosys('read real "second derivative of chebychev '//
     1           'polynomials" from lamdat',n*n,z(ddt),0,' ')
      call potntl(z(v),z(roots),n,pottyp,l,prntv)
      
         call bndcnd(z(band),points(i),order,bw,bc)
         call hdiag(z(band),z(ei),z(ev),z(work),z(dum),z(s),ia(ind),
     1              points(i),order,bw,evecs,nroots,prntvc)
         if(evecs) then
            call addend(z(ev),nroots,points(i),bc)
            call orth(z(s),z(ev),z(wts),z(runiqe),z(work),
     1                nroots,points(i))
         endif
         ei=ei+pt(i)
         ev=ev+pt(i)*pt(i)
 20   continue   
c     use richardson extrapolation technique
      if(extrap) then
         call richdn(z(eig),z(eig+pt(1)),z(evec),z(evec+pt(1)*pt(1)),
     1               stp,points(1),points(2),nroots)
      else
c
c     find better eigenvalues and eigenfunctions using these as guesses
c     and numerov method.
c
         call better(z(band),z(runiqe),z(v),z(eig),z(evec),
     1               z(s),z(wts),ia(ind),z(work),stp(1),nroots,
     2               points(1),order,bw,iter,dele,err,bc,
     3               typefd,prntsl)
      endif
      call chainx(0)
c
 1    format(/,1x,'number of finite difference points = ',i4,/,1x,
     1            'potential type                     = ',a16,/,1x,
     2            'angular momentum                   = ',i2,/1x,
     3            'type boundary condition            = ',a16,/,1x,
     4            'number of roots                    = ',i3)
 2    format(/,1x,'finite difference region begins at = ',e15.8,/,1x,
     1            'finite difference region ends at   = ',e15.8)
 3    format(/,1x,'do not perform richardson extrapolation')
      stop
      end

