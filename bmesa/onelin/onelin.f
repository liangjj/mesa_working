*deck onelin.f
c***begin prologue     onelin
c***date written       xxxxxx   (yymmdd)
c***revision date      940209   (yymmdd)
c***keywords           variation-iteration, scattering
c***author             schneider, barry (nsf)
c***source             mesa
c***purpose            solve one-dimensional scattering equations
c***                   by finite-difference/numerov method or by
c***                   variation-iteration method based on the differential
c***                   equation.
c***
c***description        potential scattering in one-dimension is solved using
c***                   a numerov procedure either directly or by the variation-
c***                   iteration method.  in both cases the problem is
c***                   recast in the inhomogeneous form 
c***                   ( E - H ) Y = V Y0 where Y0 is a free wave, V is
c***                   the potential and Y the scattered wave part of
c***                   the solution.  this apparently trivial modification
c***                   enables us to impose log-derivative boundary
c***                   conditions which do not depend on the phase shift
c***                   or k-matrix.
c***
c***                   when the iteration-variation method is used, the
c***                   numerov method is used to solve for the iterates
c***                   based on a free particle zeroth-order H.
c***
c***                   this code grew out of a bound state finite 
c***                   difference code and was tested using boundary 
c***                   conditions of specifiying a function value or 
c***                   its derivative at large values of the radial 
c***                   variable.  this is appropriate for a bound state
c***                   problem.  for scattering the formulation above
c***                   must be used. 
c
c***references         
c
c***routines called 
c***                   
      program onelin
      implicit integer(a-z)
      parameter( noen=200 )
      dimension energy(noen)
      common ia(1)
      real*8 z(1), fpkey, rmin, rmax, rdel, energy, value 
      real*8 convg, ovtol, refe
      character*4096 ops
      character*16 cpass, chrkey, pottyp, drive, bcond, gtype, mgrid
      character*16 slndir
      character*3 itoc
      character*800 card
      character*80 title
      logical logkey, varit, urefe, ink
      equivalence (ia(1),z(1))
      common /memory/ ioff
      common/io/inp,iout
c
      call drum
      call iosys ('read character options from rwf',-1,0,0,ops)
      call posinp('$diffit',cpass)
      call cardin(card)
      rmin=fpkey(card,'minimum-r-value',1.d-20,' ')          
      rmax=fpkey(card,'maximum-r-value',10.d+00,' ')          
      points=intkey(card,'number-of-points',100,' ')
c**********************************************************************c
c**********************************************************************c
c              ignore anything involving reference energy              c
c                                or                                    c
c                            multi-grid.                               c
c              some ideas were tried that did not work.                c
c**********************************************************************c
c**********************************************************************c
c**********************************************************************c
c               if points is not an even number make it even and       c
c               define rdel consistently. this is required for the     c
c               multi-grid formula to work easily and does not really  c
c               require much more work.                                c
c**********************************************************************c
      test=points-2*(points/2)
      if (test.ne.0) then
          points=points+1
      endif
      nener=intkey(card,'number-of-energies',1,' ')
      ink=logkey(card,'k-values-entered',.false.,' ')
      if(ink) then
         call fparr(card,'k-values',energy,nener,' ')
         do 5 i=1,nener
            energy(i)=energy(i)*energy(i)
 5       continue   
      else
         call fparr(card,'energies',energy,nener,' ')
      endif
      rdel=(rmax-rmin)/points
      pottyp=chrkey(card,'potential-type','exponential',' ')
      value=fpkey(card,'value',1.d0,' ')
      drive=chrkey(card,'driving-term','one',' ')
c*********************************************************************c
c         we can;                                                     c
c                1. specify the function at the last point            c
c                2. set the derivative to zero                        c
c                3. set the log derivative to zero                    c
c      note: if 3. is used then the driving term must be kohn         c 
c*********************************************************************c
      bcond=chrkey(card,'boundary-condition','function',' ')
      if (bcond.eq.'log-derivative') then
          drive='kohn'
      endif
      iter=intkey(card,'number-of-iterations',50,' ')
      maxvec=intkey(card,'maximum-number-of-vectors',50,' ')
      convg=fpkey(card,'convergence-criterion',1.d-09,' ')
      ovtol=fpkey(card,'overlap-tolerance',1.d-08,' ')
      prnt=logkey(card,'print',.false.,' ')
      varit=logkey(card,'variation-iteration',.false.,' ')
      gtype=chrkey(card,'variation-iteration=precondition',
     1             'unit-matrix',' ')
      mgrid=chrkey(card,'variation-iteration=grid','single-grid',' ')
      slndir=chrkey(card,'variation-iteration=method','matrix',' ')
      urefe=logkey(card,'reference-energy-method',.false.,' ')
      refe=fpkey(card,'reference-energy',0.d0,' ')
c**********************************************************************c
c                            set up mesh                               c
c                   a larger number of points are computed             c
c                         than used in the fd equation                 c
c**********************************************************************c
      m=points+10
      dim=m+1
      ioff=1
      do 20 i=1,2
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
   20 continue     
      title='           solve simple schroedinger equation '//
     1      'via finite difference'
      write(iout,*) title
      title='direct lu decomposition'
      if (varit) then
          title='variation-iteration'
      endif 
      write(iout,1) rmin,rmax,rdel,m,title
      do 30 ene=1,nener
c
c         energy is input in Rydbergs
c      
c         call posinp('$energy-'//itoc(ene),cpass)
c         call cardin(card)         
c         energy=fpkey(card,'energy',1.d0,' ')
         if (.not.urefe) then
             refe=energy(ene)
         endif    
         if (varit) then
             call vdrive(z(v),z(x),z(g),z(f),z(diag),z(sudiag),
     1                   z(spdiag),z(rhs),z(guess),z(exvc),z(exiter),
     2                   z(aold),z(anew),z(bold),z(bnew),z(temp),
     3                   ia(ipvt),rmin,rmax,rdel,energy(ene),refe,convg,
     4                   ovtol,value,points,m,iter,pottyp,drive,bcond,
     5                   gtype,mgrid,slndir,urefe,ops)
         else
             call ddrive(z(v),z(x),z(g),z(f),z(diag),z(sudiag),
     1                   z(spdiag),z(rhs),rmin,rmax,rdel,energy(ene),
     1                   value,m,points,pottyp,drive,bcond)
         endif
   30 continue
c
      call chainx(0)
c
    1 format(/,6x,'integration parameters:',/,20x,'r min =     ',e15.8,
     1         1x,'r max = ',e15.8,/,20x,'step size = ',e15.8,1x,
     2            'number steps = ',1x,i5,/,6x,'method of solution:',
     3             /,20x,a80,//)           
      stop
      end
