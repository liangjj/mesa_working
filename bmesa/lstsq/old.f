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
      call copy(rhs,tvec,n)
      discnt=0
      do 10 iter=1,maxit
c
c        do a gauss-seidel iteration
c
         call seidel(ham,rhs,tvec,vec,n,1)
         call vsub(tvec,tvec,vec,n)
         write(iout,*) 'I am Here'
         write(iout,*) tvec
         write(iout,*) vec
         error = sqrt( sdot(n,tvec,1,tvec,1)/sdot(n,vec,1,vec,1) )
c
c        check convergence
c
         write(iout,3) iter, error
         if(error.le.cnvrg) then
            write(iout,4)
            return
         elseif(error.le.eswtch) then
            discnt=discnt+1
            call copy(vec,pvec(1,discnt),n)
c        operate with hamiltonian on this vector
c
            call honv(ham,pvec(1,discnt),hpvec(1,discnt),n,1)
c
c        update the curent coefficient matrix and right hand side.
c            
            do 20 jter=1,discnt
               b(jter,discnt) = sdot(n,hpvec(1,jter),1,
     1                                 hpvec(1,discnt),1)
               b(discnt,jter) = b(jter,discnt)
 20         continue
c
            do 30 jter=1,discnt
               do 40 kter=1,jter
                  btmp(jter,kter) = b(jter,kter)
                  btmp(kter,jter) = b(kter,jter)
 40            continue
 30         continue   
c
            sol(discnt)=sdot(n,hpvec(1,discnt),1,rhs,1)
            call copy(sol,soltmp,discnt)
c
c        solve linear equations
c
            call sgefa(btmp,maxit,discnt,ipvt,info)
            if(info.ne.0) then
               call lnkerr('matrix singular. quit')
            endif
            call sgesl(btmp,maxit,discnt,ipvt,soltmp,0)
c                
c
c        calculate solution
c
            call ebcxx(vec,pvec,soltmp,n,discnt,1,n,n,maxit)
c
c        calculate hamiltonian on solution
c
            call ebcxx(tvec,hpvec,soltmp,n,discnt,1,n,n,maxit)
c
c        calculate residual
c
            error=0.d0
            do 50 i=1,n
               resid(i) = tvec(i) - rhs(i)
               error = error + resid(i)*resid(i)
 50         continue
            error=sqrt(error)
c
c        check convergence
c
            write(iout,3) iter, error
            if(error.le.cnvrg) then
               write(iout,4)
               return
            endif  
         else
            call copy(vec,tvec,n)
         endif
c
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
