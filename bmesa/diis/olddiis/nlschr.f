*deck nlschr.f
c***begin prologue     nlschr
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           non-linear schroedinger equation, gross-pitaevski
c***author             schneider, barry (nsf)
c***source             timeprp
c***purpose            calculate eigenvalues of gross-pitaevski equation
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       nlschr
      subroutine nlschr(f,eig,rho,psi,err,pndvr,eigr,ham0,b,btmp,
     1                  sol,work1,work2,table,ipvt,natmin,natmax,hbar,
     2                  amass,omega,scatl,todiis,switch,cnverg,n,
     3                  maxit,nstep,dim,trunc,prnt,drctv)
      implicit integer (a-z)
      real*8 f, eig, rho, psi, err, pndvr, eigr, ham0, b, btmp, sol
      real*8 work1, work2, table, natmin, natmax, natoms
      real*8 hbar, amass, omega, scatl, todiis, cnverg, error, test
      real*8 gamma, norm0, pi, zero, one, delat, eigtr
      real*8 eig0, crit
      character*80 title
      logical drctv, prnt, flag, switch
      dimension f(n,n,maxit+1), eig(n), rho(n,n), psi(n)
      dimension err(n,n,maxit+1), pndvr(n,0:n-1), eigr(n), ham0(n,n)
      dimension b(maxit+1,maxit+1), btmp(maxit+1,maxit+1), sol(maxit+1)
      dimension work1(*), work2(*), table(dim,3), ipvt(n)
      common/io/inp, iout
      data pi/3.14159265358979323844d0/
      data zero, one/ 0.d0, 1.d0/
      if(switch) then
         todiis=1.d+10
      endif
      nobj=n*n
      nerr=nobj 
c      norm0=sqrt(2.d0*amass*omega/hbar)/(4.d0*pi)
c      norm0=sqrt(norm0)
      norm0=sqrt(2.d0*amass/hbar)/(4.d0*pi)
      norm0=sqrt(norm0)
      write(iout,1)
      write(iout,2) amass, omega, scatl
      write(iout,3) maxit, cnverg
      natoms=natmin
      delat=(natmax - natmin)/nstep
c
c     ham0 contains the kinetic energy plus harmonic trap potential
c
c     put the zero-atom ( trap ) wavefunction as the starting guess
c
      call init(ham0,f(1,1,1),eig,work2,psi,n)
      eig0=eig(1)
      count=0
      do 1000 num=1,nstep
         crit = todiis
         test = crit
         gamma=8.d0*pi*scatl*natoms
         write(iout,4) natoms
         write(iout,5) eig(1)
         eig0=eig(1)
         begin=1
         do 10 iter=1,maxit
c
c           form current approximation to fock and error matrix.
c
            call fock(f(1,1,begin),rho,eigr,err(1,1,begin),pndvr,ham0,
     1                psi,work1,norm0,gamma,test,n,begin,prnt)
            if(drctv) then
c
c                  diis procedure when appropriate
c            
               if(test.le.crit) then
                  call diis(work1,work2,f,err,b,btmp,psi,norm0,test,
     1                      cnverg,ipvt,begin,nobj,nerr,trunc,
     2                      maxit,prnt,flag)
                  if(flag) then
                     begin=1
                     call copy(work1,f(1,1,begin),nobj)
                  else
                     begin=begin+1
                     call copy(work1,f(1,1,begin),nobj)   
                  endif
               endif
               error=test   
            endif
c
c           diagonalize fock matrix for the current iteration
c           and find the new wavefunction
c
            call newfn(f(1,1,begin),eig,psi,work1,n)
            if(drctv) then
               error=test
            else
               error=sqrt( (eig0-eig(1))*(eig0-eig(1)) )
            endif
            write(iout,6) iter, eig(1), error
            if(error.le.cnverg) then
               count=count+1
               table(count,1)=natoms
               table(count,2)=eig(1)
               table(count,3)=error
               eigtr=eig(1)/omega
               write(iout,7) iter, natoms, eig(1), eigtr
               if(prnt) then
                  title='full spectrum'
                  call prntrm(title,eig,n,1,n,1,iout)
                  title='final wavefunction'
                  call prntrm(title,work1,n,1,n,1,iout)
               endif
               go to 1001
            endif
            eig0=eig(1)   
 10      continue
         write(iout,8) iter
 1001    natoms = natoms + delat
 1000 continue   
      write(iout,60)
      do 2000 i=1,count
         write(iout,70) table(i,1), table(i,2), table(i,3)
 2000 continue   
      return
    1 format(/,5x,'diagonalize time-independent non-linear schroedinger'
     1            ' equation')
    2 format(/,5x,'atomic mass        = ',e15.8,/,5x,
     1            'trap frequency     = ',e15.8,/,5x,
     2            'scattering length  = ',e15.8)
    3 format(/,5x,'maximum number of iterations = ',i4,/,5x,
     1            'convergence criterion        = ',e15.8)
    4 format(/,5x,'number of atoms = ',e15.8)
    5 format(/,5x,'starting approximation',2x,'energy = ',e15.8)
    6 format(/,5x,'iteration = ',i3,1x,'energy = ',e15.8,1x,
     1            'erms = ',e15.8)         
    7 format(/,5x,'converged after ',i4,' iterations',/,5x,
     1            'number of atoms = ',e15.8,1x,'energy = ',e15.8,/,5x,
     2            'energy in trap units = ',e15.8)
    8 format(/5x,'no convergence after ',i4,' iterations')    
 60   format(//,25x,'summary table',//,10x,'number of atoms',
     1           7x,'energy',12x,'convergence')
 70   format(7x,e15.8,5x,e15.8,5x,e15.8)
      end       
