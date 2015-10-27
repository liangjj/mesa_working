*deck drvdis.f
c***begin prologue     drvdis
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           diis, iterative solve
c***author             schneider, barry (nsf)
c***source             
c***purpose            driver for diis.
c***                   
c***references         
c
c***routines called    
c***end prologue       drvdis
      subroutine drvdis(ham,rhs,pvec,dlpvec,vec,dvec,b,btmp,sol,
     1                  cnvrg,ipvt,n,maxit,itrslv)
      implicit integer (a-z)
      real*8 ham, rhs, pvec, dlpvec, vec, dvec, b, btmp, sol
      real*8 tmp, error, cnvrg, sdot
      character*(*) itrslv
      dimension ham(n,n), rhs(n), pvec(n,maxit+1), dlpvec(n,maxit+1) 
      dimension vec(n), dvec(n), b(maxit+1,maxit+1)
      dimension btmp(maxit+1,maxit+1)
      dimension sol(maxit+1), ipvt(maxit+1)
      common/io/inp, iout
      write(iout,1) maxit, cnvrg
      write(iout,2)      
c     set up right hand side as unit vector here
      call rzero(rhs,n)
      rhs(1)=1.d0
c     set up first vector as rhs
      iter=1
      call copy(rhs,pvec(1,iter),n)
      do 10 iter=1,maxit
c      
c               do a gauss-seidel iteration for next vector
c
         call copy(pvec(1,iter),pvec(1,iter+1),n)
         do 20 i=1,n
            tmp=rhs(i) +ham(i,i)*pvec(i,iter+1)
            do 30 j=1,n
               tmp = tmp - ham(i,j)*pvec(j,iter+1)
 30         continue
            pvec(i,iter+1) = tmp/ham(i,i)
 20      continue
c            
c               form error vector
c
         do 40 i=1,n
            dlpvec(i,iter) = pvec(i,iter+1) - pvec(i,iter)
 40      continue
         if(itrslv.eq.'diis') then
c
c 
            call rdiis(pvec,dlpvec,vec,dvec,b,btmp,sol,error,ipvt,iter,
     1                 n,maxit)
         else
            error = sqrt(sdot(n,dlpvec(1,iter),1,dlpvec(1,iter),1))
         endif   
        write(iout,3) iter, error
        if(error.le.cnvrg) then
           return
        endif   
 10   continue
      return
 1    format(/,1x,'beginning diis solution of linear equations',/,10x,
     1            'maximum number of iterations  = ',i4,/,10x,
     2            'convergence criterion for rms = ',e15.8)      
 2    format(/,12x,'   iteration   ','   rms error   ')
 3    format(17x,i3,6x,e15.8)                
      end       
