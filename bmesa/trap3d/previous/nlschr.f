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
      subroutine nlschr(f,eig,rho,psi,err,vnl,p1,p2,p3,eigr1,eigr2,
     1                  eigr3,ham01,ham02,ham03,v,b,btmp,sol,work1,
     2                  work2,work3,table,ipvt,natmin,natmax,hbar,amass,
     3                  omega,scatl,todiis,cnverg,n,maxit,nstep,
     4                  number,ndim,nd,trunc,prnt,drctv,incore,onoff)
      implicit integer (a-z)
      real*8 f, eig, rho, psi, err, vnl, p1, p2, p3, eigr1, eigr2, eigr3
      real*8 ham01, ham02, ham03, v
      real*8 b, btmp, sol
      real*8 work1, work2, work3, table, natmin, natmax, natoms
      real*8 hbar, amass, omega, scatl, todiis, cnverg, error, test
      real*8 gamma, norm0, pi, zero, one, delat
      real*8 eig0, crit, eigtr, totrun, u0
      character*80 title
      character*3 itoc
      logical drctv, prnt, flag, incore, onoff
      dimension nd(3)
      dimension f(n,n,*), eig(n), rho(n,n)
      dimension err(n,n,*), vnl(n,*), psi(n,*)
      dimension p1(nd(1),0:nd(1)-1), p2(nd(2),0:nd(2)-1)
      dimension p3(nd(3),0:nd(3)-1)
      dimension eigr1(nd(1)), eigr2(nd(2)), eigr3(nd(3))
      dimension ham01(nd(1),nd(1)), ham02(nd(2),nd(2))
      dimension ham03(nd(3),nd(3)), v(n)
      dimension b(maxit+1,maxit+1), btmp(maxit+1,maxit+1), sol(maxit+1)
      dimension work1(*), work2(*), work3(*), table(number,3), ipvt(n)
      common/io/inp, iout
      data pi/3.14159265358979323844d0/
      data zero, one/ 0.d0, 1.d0/
      nobj=n*n
      nerr=nobj
      norm0=1.d0/(4.d0*pi)
      norm0=sqrt(norm0)
      totrun=omega*hbar
      u0=4.d0*pi*hbar*hbar*scatl/amass
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
      do 10 i=1,n
         ham01(i,i) = ham01(i,i) + v(i)
 10   continue   
      call init(ham01,f(1,1,1),eig,work2,psi(1,1),totrun,n)
      eig0=eig(1)
      count=0
      do 1000 num=1,nstep
         crit = todiis
         test = crit
         gamma=natoms*u0
         write(iout,4) natoms, gamma
         write(iout,5) eig(1)/totrun
         eig0=eig(1)
         begin=1
         do 20 iter=1,maxit
            indx=begin
            if(.not.incore) then
               indx=1
            endif
            call smul(psi(1,indx),psi(1,indx),norm0,n)
            call potnl(vnl(1,indx),eigr1,p1,psi(1,indx),gamma,n)
c            call dnmat(rho,psi(1,indx),n)
            if(onoff) then
               call fock(f,ham01,vnl(1,indx),n,begin,prnt)
c               call errmat(f,err,rho,test,n,begin,prnt)
               call errmat(f,err,psi(1,indx),test,n,begin,prnt)
            else
               call fock(f(1,1,indx),ham01,vnl(1,indx),n,begin,prnt)
c               call errmat(f(1,1,indx),err(1,1,indx),rho,test,n,
c     1                     begin,prnt)
               call errmat(f(1,1,indx),err(1,1,indx),psi(1,indx),test,n,
     1                     begin,prnt)
            endif
            error=test/(norm0*norm0*totrun)
            if(.not.incore) then
               call iosys('write real "object vector iteration = '
     1                    //itoc(begin)//'" to lamdat',n*n,f,
     2                      0,' ')                 
               call iosys('write real "error vector iteration = '
     1                    //itoc(begin)//'" to lamdat',n*n,
     2                      err,0,' ')                 
            endif
            if(drctv) then
c
c                  diis procedure when appropriate
c            
               if(test.le.crit) then
                  if(onoff) then
                     call ddiis(work1,work2,f,err,ham01,vnl,psi,rho,
     1                          work3,b,btmp,sol,totrun,test,cnverg,
     2                          ipvt,begin,n,trunc,maxit,prnt,flag)
                  else
                     call diis(work1,work2,f,err,b,btmp,sol,totrun,test,
     1                         cnverg,ipvt,begin,nobj,nerr,trunc,
     2                         maxit,prnt,flag,incore)
                  endif
                  if(incore) then
                     if(flag) then
                        begin=1
                        indx=begin
                        if(onoff) then
                           call copy(work1,f,nobj)
                        else
                           call copy(work1,f(1,1,begin),nobj)
                        endif
                     else
                        begin=begin+1
                        indx=begin
                        if(onoff) then
                           call copy(work1,f,nobj)
                        else
                           call copy(work1,f(1,1,begin),nobj)
                        endif
                     endif
                  else
                     if(flag) then
                        begin=1
                     else
                        begin=begin+1
                     endif
                     call copy(work1,f,nobj)
                  endif
               endif
               error=test   
            endif
c
c           diagonalize fock matrix for the current iteration
c           and find the new wavefunction
c
            if(onoff) then
               call newfn(f,eig,psi(1,indx),work1,n)
            else
               call newfn(f(1,1,indx),eig,psi(1,indx),work1,n)
            endif
            if(drctv) then
               error=test
            endif
            write(iout,6) iter, eig(1)/totrun, error
            if(error.le.cnverg) then
               eigtr=eig(1)/totrun
               count=count+1
               table(count,1)=natoms
               table(count,2)=eigtr
               table(count,3)=error
               write(iout,7) iter, natoms, eig(1), eigtr
               if(prnt) then
                  title='full spectrum'
                  call prntfm(title,eig,n,1,n,1,iout)
                  title='final wavefunction'
                  call prntfm(title,work1,n,1,n,1,iout)
               endif
               go to 1001
            endif
            eig0=eig(1)   
 20      continue
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
    4 format(/,5x,'number of atoms = ',e15.8,/,5x,
     1            'gamma           = ',e15.8)
    5 format(/,5x,'starting approximation',2x,'energy = ',e15.8)
    6 format(/,5x,'iteration = ',i3,1x,'energy = ',e15.8,1x,
     1            'erms = ',e15.8)         
    7 format(/,5x,'converged after ',i4,' iterations',/,5x,
     1            'number of atoms = ',e15.8,1x,'energy = ',e15.8,
     2       /,5x,'energy in trap units = ',e15.8)
    8 format(/5x,'no convergence after ',i4,' iterations')    
 60   format(//,25x,'summary table',//,10x,'number of atoms',
     1           7x,'energy',12x,'convergence')
 70   format(7x,e15.8,5x,e15.8,5x,e15.8)
      end       
