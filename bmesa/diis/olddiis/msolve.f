*deck lsolve.f
c***begin prologue     lsolve
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       lsolve
      subroutine lsolve(hamin,ham,rhs,energy,cnverg,ipvt,n,maxit)
      implicit integer (a-z)
      real*8 hamin, ham, rhs, energy, cnverg
      dimension hamin(n,n), ham(n,n), rhs(n,n), ipvt(n)
      common/io/inp, iout
      write(iout,1) energy
      call rzero(ham,n*n)
      do 10 i=1,n
         do 20 j=1,i
            ham(i,j) = - hamin(i,j) 
            ham(j,i) = ham(i,j)
 20      continue
         ham(i,i) = ham(i,i) + energy   
 10   continue
      call sgefa(ham,n,n,ipvt,info)
      if(info.ne.0) then
         call lnkerr('error in linear solve:singular matrix')
      endif
      do 30 i=1,nrhs
         call sgesl(ham,n,n,ipvt,rhs(1,i),0)   
 30   continue   
      nbeg=1
      nend=1
      call copy(ham(1,1),psi,n)
      do 10 iter=1,maxit
c      
c        orthonormalize the new trials to the old vectors
c       
         call gschmt(p,thresh,n,nbeg,nend,nout,.true.)
         if(nout.ne.0) then
            nend=nbeg+nout-1
            write(iout,3) iter, nend            
c
c           form current approximation to non-linear 
c                    schroedinger equation
c
            call smul(psi,psi,norm0,n)
            nrm=sdot(n,psi,1,psi,1)
            call ebc(vnl,p,psi,n,n,1)
            do 20 i=1,n
               vnl(i)=vnl(i)*vnl(i)/(eigr(i)*eigr(i))
 20         continue
            call copy(ham0,ham,n*n)
            do 30 i=1,n
               ham(i,i)=ham(i,i) + gamma*vnl(i)
 30         continue                    
            delta=sqrt( (eig(1)-eig0)*(eig(1)-eig0) )
            write(iout,4) iter, eig(1), nrm, delta
            eig0=eig(1)
  
            call tred2(n,n,ham,eig,vnl,ham)
            call tql2(n,n,eig,vnl,ham,ierr)

            if(delta.le.cnvrg) go to 50
 20      continue
         write(iout,5) ipass, iter
         return
 50      write(iout,6) ipass, iter, eig(1) 
 10   continue
c     compare numerical derivative at origin with analytic value.
      deriv=sdot(n,der0,1,ham,1)
      begin=sqrt(eig(1)/(natoms*gamma+1.d0))
      write(iout,7) begin, deriv
      return
    1 format(/,5x,'energy = ',e15.8)
    2 format(/,5x,'maximum number of iterations = ',i4,/,5x,
     1            'convergence criterion        = ',e15.8)    
    3 format(/,1x,'cycle = ',i3,5x,'size of vector space = ',i3) 
    4 format(12x,i3,8x,e15.8,1x,e15.8,1x,e15.8)
    5 format(/5x,'npass = ',i3,2x,'no convergence after ',i4,
     1           ' iterations')
    6 format(/,5x,'npass = ',i3,2x,'converged after ',i4,
     1            ' iterations',/,5x,'energy = ',e15.8)    
    7 format(/,5x,'exact derivative at origin       = ',e15.8,/,5x,
     1            'approximate derivative at origin = ',e15.8)   
      end       
