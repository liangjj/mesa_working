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
      subroutine nlschr(f,eig,psi,hpsi,vnl,p1,p2,p3,eigr1,eigr2,
     1                  eigr3,ham01,ham02,ham03,eig01,eig02,eig03,
     2                  u01,u02,u03,v,b,btmp,hh,pp,ph,hp,sol,hbuf,
     3                  ibuf,diag,trials,etrial,pvec,hpvec,
     4                  vec,resid,dvm,dvmt,t1,t2,t3,work1,work2,table,
     5                  ipvt,indx,natmin,natmax,hbar,amass,omega,scatl,
     6                  todiis,cnverg,thresh,nd,ndim,n,nroots,ntrials,
     7                  maxit,maxvec,lenbuf,nstep,number,doit,
     8                  incore,prtrl,prnt)
      implicit integer (a-z)
      real*8 f, eig, psi, hpsi, vnl, p1, p2, p3, eigr1, eigr2, eigr3
      real*8 ham01, ham02, ham03, eig01, eig02, eig03
      real*8 u01, u02, u03, v
      real*8 b, btmp, hh, pp, ph, hp, sol
      real*8 hbuf, diag, trials, etrial, pvec, hpvec, vec, resid
      real*8 dvm, dvmt, t1, t2, t3
      real*8 work1, work2, work3, table, natmin, natmax, natoms
      real*8 hbar, amass, omega, scatl, todiis, cnverg, error, test
      real*8 gamma, norm0, pi, zero, one, delat
      real*8 eig0, crit, eigtr, totrun, u0
      character*80 title
      character*3 itoc
      character*1 it
      logical drctv, prnt, flag, prtrl
      dimension nd(3)
      dimension f(n,n), eig(n), v(n)
      dimension vnl(n,*), psi(n,*), hpsi(n,*)
      dimension p1(nd(1),0:nd(1)-1), p2(nd(2),0:nd(2)-1)
      dimension p3(nd(3),0:nd(3)-1)
      dimension eigr1(nd(1)), eigr2(nd(2)), eigr3(nd(3))
      dimension ham01(nd(1),nd(1)), ham02(nd(2),nd(2))
      dimension ham03(nd(3),nd(3))
      dimension u01(nd(1),nd(1)), u02(nd(2),nd(2)), u03(nd(3),nd(3))
      dimension b(maxit+1,maxit+1), btmp(maxit+1,maxit+1), sol(maxit+1)
      dimension hh(maxit,maxit), pp(maxit,maxit), ph(maxit,maxit)
      dimension hp(maxit,maxit)
      dimension work1(*), work2(*), table(number,3), ipvt(n), indx(n,*)
      dimension hbuf(lenbuf), ibuf(2,lenbuf), diag(n), trials(n,ntrials)
      dimension etrial(ntrials), pvec(n,maxvec), hpvec(n,maxvec)
      dimension vec(n,maxvec), resid(n,maxvec)
      dimension dvm(maxvec,maxvec), dvmt(maxvec,maxvec)
      dimension t1(*), t2(*), t3(*)
      dimension table(number,3)
      dimension prnt(10)
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
c     make some trial vectors using the eigenvalues and eigenvectors
c     of the separable hamiltonian.  This is done without diagonalizing
c     a large hamiltonian but using the subsystem information and the
c     index vector
c
      call vtrial(u01,u02,u03,eig01,eig02,eig03,trials,etrial,
     1            indx(1,ndim+1),nd,ndim,n,ntrials,prtrl)
c
c     if the eigenvalue problem is solved iteratively we need to set up
c     the buffer information.  if the problem will be solved directly,
c     we store the hamiltonian information in hbuf as well, with the
c     understanding that it is big enough to hold all matrix elements.
c
      if(doit) then     
         it=itoc(ndim)
         title='buffered hamiltonian for '//it 
     1                                   //'d dimensional hamiltonian'
         call hamil(ham01,ham02,ham03,v,hbuf,ibuf,indx,diag,lenbuf,nd,
     2              ndim,n,.false.,ntot,incore,title)
      else
         call lschr(ham01,ham02,ham03,v,hbuf,dum,indx,n,nd,ndim,.false.,
     1              .false.)
      endif
c     put the zero-atom ( trap ) wavefunction as the starting guess
c
      call copy(trials(1,1),psi(1,1),n)
      eig0=etrial(1)
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
c           
c           wavefunction
            call smul(psi(1,begin),psi(1,begin),norm0,n)
            call potnl(vnl(1,begin),eigr1,p1,psi(1,begin),gamma,n)
            call fock(f,ham01,vnl(1,begin),n,begin,prnt)
            call fpsi(f,psi(1,begin),hpsi(1,begin),test,n,begin,prnt)
            error=test/(norm0*norm0*totrun)
            call errfac(psi,hpsi,hh,pp,ph,hp,n,begin,maxit)
            if(drctv) then
c
c                  diis procedure when appropriate
c            
               if(test.le.crit) then
                  call diis(work1,f,ham01,vnl,psi,hpsi,
     1                      b,btmp,hh,pp,ph,hp,sol,test,ipvt,begin,
     2                      n,maxit,prnt)
                  error=test/totrun
               endif
            endif
            begin=begin+1
c
c           diagonalize fock matrix for the current iteration
c           and find the new wavefunction
c
            call newfn(f,eig,psi(1,begin),work1,n)
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
