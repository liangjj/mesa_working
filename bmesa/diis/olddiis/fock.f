*deck fock.f
c***begin prologue     fock
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           non-linear schroedinger equation, gross-pitaevski
c***author             schneider, barry (nsf)
c***source             tstdiis
c***purpose            form current approximation to fock matrix
c***                   and the diis error vector
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       fock
      subroutine fock(f,rho,eig,evec,pndvr,ham0,psiin,work,
     1                norm0,gamma,error,n,iter,prnt)
      implicit integer (a-z)
      real*8 f, rho, eig, evec, pndvr, ham0, psiin, work
      real*8 norm0, gamma, error
      character*3 itoc
      character*80 title
      logical prnt
      dimension f(n,n), rho(n,n), eig(n), evec(n,n), pndvr(n,0:n-1)
      dimension ham0(n,n), psiin(n), work(*)
      common/io/inp, iout
      write(iout,1) iter
c   
c     scale wavefunction by norm0 because of change of variables
c
      call smul(psiin,psiin,norm0,n)
c
c     transform wavefunction
c
      call ebc(work,pndvr,psiin,n,n,1)
c     
c     zero density matrix
c
c
c     calculate density and non-linear term of fock matrix.      
c
      do 10 i=1,n
         do 20 j=1,i
            rho(i,j) = psiin(i)*psiin(j)
            rho(j,i) = rho(i,j)
 20      continue
 10   continue            
      do 30 i=1,n
         psiin(i) = work(i)*work(i)/(eig(i)*eig(i))
 30   continue
      call smul(psiin,psiin,gamma,n) 
c    
c     copy trap hamiltonian to fock matrix
c
      call copy(ham0,f,n*n)
c
c     add in non-linear term
c
      do 40 i=1,n
         f(i,i) = f(i,i) + psiin(i)    
 40   continue   
      if(prnt) then
         title='fock matrix iteration = '//itoc(iter)
         call prntrm(title,f,n,n,n,n,iout)
      endif
c
c     calculate error matrix = f*rho - rho*f = [f,rho]
c
      call ebc(evec,f,rho,n,n,n)
      call ambc(evec,rho,f,n,n,n)
      error=0.d0
      do 50 i=1,n
         do 60 j=1,i
            error = max(error,abs(evec(i,j)))
 60      continue
 50   continue   
      error=error/(norm0*norm0)
      write(iout,2) iter, error        
      if(prnt) then
         title='error matrix iteration = '//itoc(iter)
         call prntrm(title,evec,n,n,n,n,iout)
      endif                  
      return 
   1  format(/,5x,'fock and error matrix formation iteration =', i3)
   2  format(/,5x,'iteration = ',i3,2x,'error = ',e15.8)
      end       
