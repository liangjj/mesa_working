*deck seidel.f
c***begin prologue     seidel
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
c***end prologue       seidel
      subroutine seidel(matrix,rhs,xold,xnew,iter,n,prnt)
      implicit integer (a-z)
      real*8 matrix, rhs, xold, xnew, tmp
      logical prnt
      character*80 title
      character*3 itoc
      dimension matrix(n,n), rhs(n)
      dimension xold(n), xnew(n)
      common/io/inp, iout
c     perform gauss-seidel iteration
      call copy(xold,xnew,n)
      do 10 i=1,n
         tmp = rhs(i) + matrix(i,i)*xnew(i)
         do 20 j=1,n
            tmp = tmp - matrix(i,j)*xnew(j)
 20      continue
         xnew(i) = tmp/matrix(i,i)
 10   continue   
      if(prnt) then
         title='solution matrix for iteration = '//itoc(iter)
         call prntrm(title,xnew,n,1,n,1,iout)
      endif
      return
      end       
