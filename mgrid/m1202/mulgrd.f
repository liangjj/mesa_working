*deck mulgrd 
c***begin prologue     mulgrd
c***date written       941208   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m1201, link 1200, multigrid
c***author             schneider, b. i.(nsf)
c***source             m1201
c***purpose            simple multigrid code for poisson or
c***                   inhomogeneous wave equation.
c***description        full multigrid method. solves equation
c***                   in an (x,y) region bounded by the variables
c***                   left and right in both dimensions.
c***
c***                   the equation under consideration is;
c***                   - ( del**2 + k*k ) u = f in 2-d
c***references         multigrid methods, siam, s. mccormick ed. 
c
c***routines called    iosys, util and mdutil
c***end prologue       mulgrd
c
c      (-left,right)  _ _ _ _ _ _ _ _ _ _ (right,right)
c                     _                 _
c                     _                 _
c                     _                 _         
c                     _                 _
c         (-left,0)   _      (0,0)      _ (right,0)
c                     _                 _
c                     _                 _
c                     _                 _
c      (-left,-left)  _ _ _ _ _ _ _ _ _ _ (right,-left)
c
      program mulgrd
c
      implicit integer (a-z)
      parameter ( ngmax=50 )
      character*4096 ops
      character*16 cpass, chrkey, model, type, ittyp, restr
      character*1600 card
      logical posinp, logkey, prnrho, prnsol, ecrse, tocom
      common z(1)
      real*8 z, left, right, stpg, tol, fpkey, cstep, hf, k, rms
      real*8 hsq, hsqi, ksq, fac, ifac
      dimension ia(1), u(ngmax), rhs(ngmax), rho(ngmax), res(ngmax)
      dimension nptg(ngmax), stpg(ngmax)
      equivalence (ia(1),z(1))
      common/io/inp, iout      
      common/memory/ioff
c
      call drum
      write(iout,*)
      write(iout,*) '                         m1200:multigrid program'
      call iosys ('read character options from rwf',-1,0,0,ops)
      prnrho=logkey(ops,'print=m1202=input-rho',.false.,' ')
      prnsol=logkey(ops,'print=m1202=coarse-solution',.false.,' ')
      if (posinp('$multigrid',cpass) ) then
          call cardin(card)
      endif          
      left=fpkey(card,'left-boundary',-1.d0,' ')
      right=fpkey(card,'right-boundary',1.d0,' ')
      ncorse=intkey(card,'number-of-coarse-grid-points',3,' ')
      cstep=(right-left)/(ncorse-1)
      ng=intkey(card,'number-of-grids',5,' ')
      npre=intkey(card,'number-of-pre-smoothing-relaxations',
     1                       1,' ')
      npost=intkey(card,'number-of-post-smoothing-relaxations',
     1                        1,' ')
      ncycle=intkey(card,'number-of-v-cycles',1,' ')
      maxit=intkey(card,'maximum-number-of-gauss-seidel-iterations',
     1             100,' ')
      tol=fpkey(card,'convergence-criterion',1.d-07,' ')
      model=chrkey(card,'type-rhs','model',' ')
      type=chrkey(card,'equation-type','poisson',' ')
      ittyp=chrkey(card,'iteration-method','gauss-seidel',' ')
      tocom=logkey(card,'iterate-to-convergence',.false.,' ')
      restr=chrkey(card,'restriction-method','half-weighting',' ')
      k=0.d0
      if (type.eq.'wave') then
          k=fpkey(card,'wave-number',1.d0,' ')
      endif
      ecrse=logkey(card,'exact-coarse-grid-solution',.false.,' ')
      write(iout,1) type, model
      if(type.eq.'wave') then
         write(iout,2) k
      endif
c     calculate number of points and stepsizes for each grid      
      n=ncorse
      nptg(1)=n
      stpg(1)=cstep
      do 10 i=2,ng
         n=n+n-1
         nptg(i)=n
         stpg(i)=.5d0*stpg(i-1)
 10   continue   
      write(iout,3) ng,npre,npost,ncycle,ittyp
      do 20 i=1,ng
         write(iout,4) i, nptg(i), stpg(i)
 20   continue
c     allocate memory        
      ioff=1
      do 30 i=1,2
         ibeg=ioff
         do 40 ig=1,ng
            nptgsq=nptg(ig)*nptg(ig)
            u(ig)=ibeg
            rhs(ig)=u(ig)+nptgsq
            rho(ig)=rhs(ig)+nptgsq
            res(ig)=rho(ig)+nptgsq
            ibeg=res(ig)+nptgsq
 40      continue
            uhold=ibeg
            ibeg=ibeg+nptgsq
         words=wpadti(ibeg)   
         if (i.eq.1) then
             call iosys ('read integer maxsiz from rwf',1,
     1                    canget,0,' ')
             if (words.gt.canget) then
                 call lnkerr('not enough memory. will quit')
             endif
             call iosys ('write integer maxsiz to rwf',1,
     1                    words,0,' ')
         else
             call getscm(words,z,ngot,'mgrid',0)
         endif
 30   continue
c     input exact solution and right hand sides for all grids.
c     note right hand sides are actually scaled by the square
c     of the stepsize.
      do 50 ig=1,ng
         call rhoin(z(rho(ig)),stpg(ig),ig,nptg(ig),model,prnrho)
 50   continue
      nsq=nptg(ng)*nptg(ng)
      call copy(z(rho(ng)),z(u(ng)),nsq)
      call copy(z(rho(ng)),z(uhold),nsq)
      if (ecrse) then
          write(iout,5)
          call slvsml(z(u(1)),z(rho(1)),z(res(1)),stpg(1),k,tol,
     1                nptg(1),maxit,ittyp,prnsol)
      else
          call slvsml(z(u(1)),z(rho(1)),z(res(1)),stpg(1),k,tol,
     1                nptg(1),maxit,ittyp,prnsol)
          do 60 j=2,ng
             write(iout,6) j
             call interp(z(u(j)),z(u(j-1)),nptg(j),nptg(j-1))               
             if(j.ne.ng) then
                 call copy(z(rho(j)),z(rhs(j)),nptg(j)*nptg(j)) 
             else
                 call copy(z(uhold),z(rhs(j)),nptg(j)*nptg(j))
             endif
             do 70 jcycle=1,ncycle
                write(iout,7) jcycle
                do 80 jj=j,2,-1
c               perform some relaxations on the fine grid then inject
c               the result to the next lowest coarse grid.
                   nf=nptg(jj)
                   hf=stpg(jj)
                   hsq=hf*hf
                   hsqi=1.d0/hsq
                   ksq=k*k*hsq
                   ifac=1.d0/(4.d0-ksq)
                   do 90 jpre=1,npre
                      call relax(z(u(jj)),z(rhs(jj)),hsq,ifac,nf)
 90                continue
                   fac=4.d0*hsqi-k*k
                   call resid(z(res(jj)),z(u(jj)),z(rhs(jj)),hf,
     1                        hsqi,fac,nf,rms)
                   write(iout,8) jj, rms
                   nc=nptg(jj-1) 
                   call rstrct(z(rhs(jj-1)),z(res(jj)),nc,nf,restr) 
                   call rzero(z(u(jj-1)),nc*nc) 
 80             continue
                call rzero(z(u(1)),nptg(1)*nptg(1))
                call slvsml(z(u(1)),z(rhs(1)),z(res(1)),stpg(1),k,
     1                      tol,nptg(1),maxit,ittyp,prnsol)
c        now ascend from the coarse grid to the fine grid
                write(iout,9) jcycle
                do 100 jj=2,j
                   nf=nptg(jj)
                   nc=nptg(jj-1)
                   hf=stpg(jj)
                   hsq=hf*hf
                   hsqi=1.d0/hsq
                   ksq=k*k*hsq
                   ifac=1.d0/(4.d0-ksq)
                   fac=4.d0*hsqi-k*k
                   call addint(z(u(jj)),z(u(jj-1)),z(res(jj)),nf,nc)
c           perform some relaxations for further smoothing
                   do 200 jpost=1,npost
                      call relax(z(u(jj)),z(rhs(jj)),hsq,ifac,nf)
 200               continue
                   call resid(z(res(jj)),z(u(jj)),z(rhs(jj)),hf,
     1                        hsqi,fac,nf,rms)
                   write(iout,8) jj, rms     
 100            continue
 70          continue
 60       continue
          if (tocom) then
              prnsol=.true.
              call slvsml(z(u(ng)),z(rho(ng)),z(res(ng)),stpg(ng),k,
     1                    tol,nptg(ng),maxit,ittyp,prnsol)
          endif
      endif
      call iosys ('write integer maxsiz to rwf',1,canget,0,' ')
      call chainx(0)               
      stop
 1    format(//,1x,'equation being solved = ',a16,
     1       /1x,  'right hand side = ',a16)
 2    format(/,1x,'the wave number = ',e15.8)   
 3    format(/,20x,'various grid parameters',/,5x,
     1             'number of grids = ',i3,1x,
     2             'number of pre-smoothing steps  = ',i2,/,5x,
     3             'number of post-smoothing steps = ',i2,1x,
     4             'number of v cycles = ',i2,/,5x,
     5             'iteration method = ',a16 )
 4    format(/,1x,'grid number = ',i2,1x,'number of points = ',i6,1x,
     1            'step size = ',e15.8)        
 5    format(/,1x,'iterating solution to convergence on one grid')
 6    format(/,25x,'beginning grid = ',i2)
 7    format(/,10x,'downward stroke of v cycle = ',i2)
 8    format(/,5x,'residual for grid = ',i2,2x,'rms = ',e15.8)
 9    format(/,10x,'upward stroke of v cycle = ',i2)      
      end





