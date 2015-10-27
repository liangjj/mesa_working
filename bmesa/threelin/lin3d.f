*deck lin3d.f
c***begin prologue     lin3d
c***date written       xxxxxx   (yymmdd)
c***revision date      940209   (yymmdd)
c***keywords           variation-iteration, scattering
c***author             schneider, barry (nsf)
c***source             lin3d
c***purpose            solve three-dimensional scattering equations
c***                   using finite-difference/variation-iteration 
c***                   method based on the partial differential equation.
c***
c***description        potential scattering in three-dimensions is solved
c***                   using the variation-iteration method.
c***                   the inhomogeneous partial differential equation
c***                   is solved using using various approaches.
c
c***references         
c
c***routines called 
c***                   
      program lin3d
      implicit integer(a-z)
      parameter ( ngmax=50 )
      common ia(1)
      real*8 z(1), minval, maxval, cstep, tol, energy, convg, ovtol
      character*4096 ops
      character*16 cpass, chrkey, pottyp
      character*3 itoc
      character*800 card
      character*80 title
      logical logkey, varit
      equivalence (ia(1),z(1))
      dimension nptg(ngmax), stpg(ngmax)
      common /memory/ ioff
      common/io/inp,iout
c
      call drum
      call iosys ('read character options from rwf',-1,0,0,ops)
      call posinp('$lin3d',cpass)
      call cardin(card)
      minval=fpkey(card,'minimum-value-of-coordinate',1.d-20,' ')          
      maxval=fpkey(card,'maximum-value-of-coordinate',10.d+00,' ')
      ncorse=intkey(card,'number-of-coarse-grid-points',3,' ')
      cstep=(maxval-minval)/(ncorse-1)
      ng=intkey(card,'number-of-grids',6,' ')
      npre=intkey(card,'number-of-pre-smoothing-relaxations',
     1                       1,' ')
      npost=intkey(card,'number-of-post-smoothing-relaxations',
     1                        1,' ')
      ncycle=intkey(card,'number-of-v-cycles',1,' ')
      maxit=intkey(card,'maximum-number-of-gauss-seidel-iterations',
     1             100,' ')
      tol=fpkey(card,'convergence-criterion',1.d-07,' ')
      write(iout,1)
      write(iout,2) ng,npre,npost,ncycle 
c     calculate number of points and stepsizes for each grid      
      n=ncorse
      nptg(1)=n
      stpg(1)=cstep
      if (ng.gt.1) then
          do 10 i=2,ng
             n=n+n-1
             nptg(i)=n
             stpg(i)=.5d0*stpg(i-1)
 10       continue
      endif 
      do 20 i=1,ng
         write(iout,3) i, nptg(i), stpg(i)
   20 continue         
      prnt=logkey(card,'print',.false.,' ')
      varit=logkey(card,'variation-iteration',.false.,' ')
      nener=intkey(card,'number-of-energies',1,' ')
      pottyp=chrkey(card,'potential-type','exponential',' ')
      write(iout,4) nener, pottyp
      iter=intkey(card,'number-of-iterations',50,' ')
      maxvec=intkey(card,'maximum-number-of-vectors',50,' ')
      convg=fpkey(card,'convergence-criterion',1.d-09,' ')
      ovtol=fpkey(card,'overlap-tolerance',1.d-08,' ')
      if (varit ) then
          write(iout,5) iter, maxvec, convg, ovtol
      else
          write(iout,6)
      endif
      m=points+10
      dim=m+1
      ioff=1
      do 30 i=1,2
         x=ioff
         v=x+dim
         diag=v+dim
         sudiag=diag+dim
         spdiag=sudiag+dim
         g=spdiag+dim
         rhs=g+dim
         f=rhs+dim
         words=f+dim
         if (varit) then
             guess=words
             temp=guess+dim
             exvc=temp+6*dim
             exiter=exvc+dim*maxvec
             aold=exiter+dim*maxvec
             anew=aold+iter*iter
             bold=anew+iter*iter
             bnew=bold+iter
             ipvt=wpadti(bnew+iter)
             words=ipvt+iter
             words=iadtwp(words)
         endif             
         if (i.eq.1) then
             need=wpadti(words)
             call getscm(need,z,maxcor,'onelin',0)
         endif
   30 continue
      do 40 ene=1,nener
c
c         energy is input in Rydbergs
c      
         call posinp('$energy-'//itoc(ene),cpass)
         call cardin(card)         
         energy=fpkey(card,'energy',1.d0,' ')
         if (varit) then
             call vdrive(z(v),z(x),z(g),z(f),z(diag),z(sudiag),
     1                   z(spdiag),z(rhs),z(guess),z(exvc),z(exiter),
     2                   z(aold),z(anew),z(bold),z(bnew),z(temp),
     3                   ia(ipvt),rmin,rmax,rdel,energy,refe,convg,
     4                   ovtol,value,points,m,iter,pottyp,drive,bcond,
     5                   gtype,mgrid,slndir,urefe,ops)
         else
             call ddrive(z(v),z(x),z(g),z(f),z(diag),z(sudiag),
     1                   z(spdiag),z(rhs),rmin,rmax,rdel,energy,
     1                   value,m,points,pottyp,drive,bcond)
         endif
   30 continue
c
      call chainx(0)
c
    1 format(//,10x,'solution of three dimensional schroedinger'/,2x,
     1             ' equation using variation-iteration/multigrid')
    2 format(/,20x,'various grid parameters',/,5x,
     1             'number of grids = ',i3,1x,
     2             'number of pre-smoothing steps  = ',i2,/,5x,
     3             'number of post-smoothing steps = ',i2,1x,
     4             'number of v cycles = ',i2)
    3 format(/,20x,'grid number = ',i2,1x,'number of points = ',i6,
     1       /20x, 'step size = ',e15.8)           
    4 format (/,20x,'number of energies = ',i3,1x,'potential type = '
     1             ,a16)
    5 format(/,20x,'number of iterations = ',i4,/,20x,
     1             'maximum number of vectors = ',i3,/20x,
     2             'convergence criterion = ',e15.8,/20x,
     3             'overlap tolerance = ',e15.8)
    6 format(/'equations solved by direct lu decomposition')                               
      stop
      end
