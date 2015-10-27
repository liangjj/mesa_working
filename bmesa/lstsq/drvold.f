*deck drvlst.f
c***begin prologue     drvlst
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           least squares
c***author             schneider, barry (nsf)
c***source             
c***purpose            least squares solution of linear equations.
c***                   
c***references         
c
c***routines called    
c***end prologue       drvlst
      subroutine drvlst(ham,rhs,pvec,hpvec,vec,tvec,b,btmp,sol,soltmp,
     1                  resid,cnvrg,eswtch,ipvt,n,maxit)
      implicit integer (a-z)
      real*8 ham, rhs, pvec, hpvec, vec, tvec, b, btmp, sol, soltmp
      real*8 resid, tmp, error, cnvrg, eswtch, sdot
      dimension ham(n,n), rhs(n), pvec(n,maxit), hpvec(n,maxit) 
      dimension vec(n), tvec(n), b(maxit,maxit)
      dimension btmp(maxit,maxit), sol(maxit), soltmp(maxit)
      dimension ipvt(maxit), resid(n)
      common/io/inp, iout
      write(iout,1) maxit, cnvrg
      write(iout,2)      
c     set up first vector as rhs
      call copy(rhs,pvec(1,1),n)
      do 10 iter=1,maxit

c
c        operate with hamiltonian on this vector
c
         call honv(ham,pvec(1,iter),hpvec(1,iter),n,1)
c
c        update the curent coefficient matrix and right hand side.
c            
         do 20 jter=1,iter
            b(jter,iter) = sdot(n,hpvec(1,jter),1,
     1                              hpvec(1,iter),1)
            b(iter,jter) = b(jter,iter)
 20      continue
c
         do 30 jter=1,iter
            do 40 kter=1,jter
               btmp(jter,kter) = b(jter,kter)
               btmp(kter,jter) = b(kter,jter)
 40         continue
 30      continue   
c
         sol(iter)=sdot(n,hpvec(1,iter),1,rhs,1)
         call copy(sol,soltmp,iter)
c
c        solve linear equations
c
         call sgefa(btmp,maxit,iter,ipvt,info)
         if(info.ne.0) then
            call lnkerr('matrix singular. quit')
         endif
         call sgesl(btmp,maxit,iter,ipvt,soltmp,0)
c                
c
c        calculate solution
c
         call ebcxx(vec,pvec,soltmp,n,iter,1,n,n,maxit)
c
c        calculate hamiltonian on solution
c
         call ebcxx(tvec,hpvec,soltmp,n,iter,1,n,n,maxit)
c
c        calculate residual
c
         error=0.d0
         do 50 i=1,n
            resid(i) = tvec(i) - rhs(i)
            error = error + resid(i)*resid(i)
 50      continue
         error=sqrt(error)
c
c        check convergence
c
         write(iout,3) iter, error
         if(error.le.cnvrg) then
            write(iout,4)
            return
         endif  
c
c        do a gauss-seidel iteration
c
         call seidel(ham,rhs,vec,pvec(1,iter+1),n,1)
c           
 10   continue
      return
 1    format(/,1x,'least squares solution of linear equations',/,10x,
     1            'maximum number of iterations  = ',i4,/,10x,
     2            'convergence criterion for rms = ',e15.8)      
 2    format(/,12x,'   iteration   ','   rms error   ')
 3    format(17x,i3,6x,e15.8)                
 4    format(/,1x,'solution has converged to tolerance')
      end       
