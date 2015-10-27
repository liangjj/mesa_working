*deck linsys.f
c***begin prologue     linsys
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
c***end prologue       linsys
      subroutine linsys(ham,rhs,energy,matrix,xold,err,xcur,errcur,b,
     1                  btmp,sol,ipvt,todiis,switch,cnverg,n,m,maxit,
     2                  trunc,prnt,drctv,type,dirslv)
      implicit integer (a-z)
      complex*16 ham, rhs
      complex*16 matrix, xold, err, b, btmp, sol
      complex*16 xcur, errcur
      real*8 todiis, cnverg, error, test, crit
      real*8 energy
      character*80 title
      logical drctv, prnt, flag, switch, type, dirslv
      dimension ham(n,n), rhs(n,m), matrix(n,n), xold(n,maxit+1)
      dimension err(n,maxit+1), b(maxit+1,maxit+1)
      dimension btmp(maxit+1,maxit+1)
      dimension sol(maxit+1), ipvt(n), xcur(n), errcur(n)
      common/io/inp, iout
      if(dirslv) then
         call prep(ham,matrix,energy,n,prnt)
         call cgefa(matrix,n,n,ipvt,info)
         do 10 i=1,m
            call cgesl(matrix,n,n,ipvt,rhs(1,m),0)
 10      continue   
         title='solution to complex linear system'
         call prntcm(title,rhs,n,m,n,m,iout)       
      else
         if(switch) then
            todiis=1.d+10
         endif
         nobj=n
         nerr=nobj 
         write(iout,1)
         write(iout,2) maxit, cnverg, energy
c
c     ham contains the hamiltonian
c
c
c     prepare the system for the iteration
c
         call prep(ham,matrix,energy,n,prnt)
c      
         do 1000 num=1,m
c
            write(iout,3) num
c                     
c
c     put the rhs as the starting guess
c
            call cc2opy(rhs(1,num),xold(1,1),n)
            if(prnt) then
               title='initial guess'
               call prntcm(title,xold(1,1),n,1,n,1,iout)
            endif
            crit=todiis
            test=crit
            begin=1
            do 20 iter=1,maxit
c
c        compute the error vector and the maximum error
c
               call errvec(matrix,rhs(1,num),xold(1,begin),err(1,begin),
     1                     test,iter,n,prnt)
                  if(drctv) then
c
c           diis procedure when appropriate
c            
                  if(test.le.crit) then
                     call diis(xcur,errcur,xold,err,b,btmp,sol,test,
     1                         cnverg,ipvt,begin,nobj,nerr,trunc,
     2                         maxit,prnt,flag)
                     if(flag) then
                        begin=1
                        call cc2opy(xcur,xold(1,begin),nobj)
                     else
                        begin=begin+1
                        call cc2opy(xcur,xold(1,begin),nobj)   
                     endif
                  endif
               endif
               error=test
               write(iout,4) iter, begin, error
               if(error.le.cnverg) then
                  go to 30
               endif
c
c           refine the solution using gauss-seidel
c
               call seidel(matrix,rhs(1,num),xold(1,begin),xcur,iter,
     1                     n,type,prnt)
               call cc2opy(xcur,xold(1,begin),n)
 20         continue
            write(iout,6) iter
            return
 30         write(iout,5) iter
            title='converged solution'
            call prntcm(title,xcur,n,1,n,1,iout)
 1000    continue  
      endif
      return
    1 format(/,5x,'solve linear system using iteration and diis')
    2 format(/,5x,'maximum number of iterations = ',i4,/,5x,
     1            'convergence criterion        = ',e15.8,/,5x,
     2            'energy                       = ',e15.8)
    3 format(/,15x,'starting procedure for solution = ',i3) 
    4 format(/,5x,'actual iteration = ',i3,1x,'diis counter = ',i3,1x,
     1            'erms = ',e15.8)         
    5 format(/,5x,'converged after ',i4,' iterations')
    6 format(/5x,'no convergence after ',i4,' iterations')    
      end       
