*deck mulgrd 
c***begin prologue     mulgrd
c***date written       941208   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m1200, link 1200, multigrid
c***author             schneider, b. i.(nsf)
c***source             m1200
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
      character*16 cpass, chrkey, model, type, ittyp
      character*1600 card
      logical posinp, logkey, prnrho, prnsol, ecrse
      common z(1)
      real*8 z, left, right, stpg, tol, fpkey, cstep, h, k
      dimension ia(1), u(ngmax), f(ngmax), ex(ngmax), res(ngmax)
      dimension nptg(ngmax), stpg(ngmax)
      equivalence (ia(1),z(1))
      common/io/inp, iout      
      common/memory/ioff
c
      call drum
      write(iout,*)
      write(iout,*) '                         m1200:multigrid program'
      call iosys ('read character options from rwf',-1,0,0,ops)
      prnrho=logkey(ops,'print=m1200=input-rho',.false.,' ')
      prnsol=logkey(ops,'print=m1200=coarse-solution',.false.,' ')
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
            nn=nptg(ig)*nptg(ig)
            u(ig)=ibeg
            f(ig)=u(ig)+nn
            ex(ig)=f(ig)+nn
            res(ig)=ex(ig)+nn
            ibeg=res(ig)+nn
 40      continue
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
      do 50 ig=1,ng
         nn=nptg(ig)*nptg(ig)
         call rzero(z(u(ig)),nn)
         call rzero(z(f(ig)),nn)
         call rzero(z(ex(ig)),nn)
         call rzero(z(res(ig)),nn)
 50   continue
c     input exact solution and right hand sides for all grids.
c     note right hand sides are actually scaled by the square
c     of the stepsize.
      do 60 ig=1,ng
         call rhoin(z(f(ig)),z(ex(ig)),stpg(ig),ig,nptg(ig),
     1              model,prnrho)
 60   continue
      if (ecrse) then
          write(iout,5)
          call rzero(z(u(1)),nptg(1)*nptg(1))
          call slvsml(z(u(1)),z(f(1)),stpg(1),k,tol,nptg(1),
     1                maxit,ittyp,prnsol)
          call cmpare(z(u(1)),z(ex(1)),z(f(1)),stpg(1),k,nptg(1))
      else
          do 70 ig=1,ng
c        descend from finest grid to coarse grid level two.
c        if outer loop is at first grid just skip over this and
c        relax the solution to convergence if desired.
             if (ig.gt.1) then
                 do 80 jg=ig,2,-1
                    nf=nptg(jg)
                    nc=nptg(jg-1)
                    h=stpg(jg)
c               perform some relaxations on the fine grid then inject
c               the result to the next lowest coarse grid.
                    do 90 cycle=1,npre
                       call relax(z(u(jg)),z(f(jg)),h,k,nf,ittyp)
 90                 continue
                    call resinj(z(u(jg)),z(f(jg)),z(res(jg)),
     1                          z(f(jg-1)),h,k,nf,nc)
 80              continue
             endif
c        now ascend from the coarse grid to the fine grid
             do 100 jg=1,ig
                nc=nptg(jg)
                h=stpg(jg)
c           perform some relaxations for further smoothing
                do 200 cycle=1,npost
                   call relax(z(u(jg)),z(f(jg)),h,k,nc,ittyp)
 200            continue
c           check for various errors and norms
                if(jg.eq.ig) then
                   call cmpare(z(u(jg)),z(ex(jg)),z(f(jg)),h,k,nc)
                endif                
                if (jg.lt.ng) then
                    nf=nptg(jg+1)
                    call intadd(z(u(jg)),z(u(jg+1)),nc,nf)
                endif
 100         continue
 70       continue
          call slvsml(z(u(ng)),z(f(ng)),stpg(ng),k,tol,
     1                nptg(ng),maxit,ittyp,prnsol)
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
      end





