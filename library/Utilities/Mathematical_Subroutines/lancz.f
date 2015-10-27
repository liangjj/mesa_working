*deck lancz.f
c***begin prologue     lancz
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           eigenvalues, eigenvectors
c***author             schneider, barry (nsf)
c***source
c***purpose            lanczos diagonalization.
c***                   
c***description        b v  = (mat - a ) v   - b    v  
c***                    j j           j   j-1   j-1  j-2       
c***                                                                       
c***                                                          
c***                   v   = matrix of dimension iter*n
c***                         on entry first vector must be filled with
c***                         the starting vector to prime the Lanczos
c***                         recursion.
c***                   mat = the matrix you are trying to tridiagonalize
c***                   a,b = the recursion coefficients
c***                   scr = scratch vector of length n
c**                    vec = eigenvectors of tridiagonal matrix iter*iter
c***                   n   = vector length
c***                   iter = number of iterations
c***                   note as with all simple Lanczos methods, one can
c***                   get linear dependence problems and multiple copies
c***                   of eigenvalues if iter is large.
c***references         
c
c***routines called    
c***end prologue       lancz
      subroutine lancz(v,mat,a,b,vec,scr,n,iter)
      implicit integer (a-z)
      real*8 v, mat, a, b, vec, scr, sdot, anorm
      character*80 title
      dimension v(n,0:iter), mat(n,n), vec(iter,iter), scr(n)
      dimension a(1:iter), b(1:iter)
      common/io/inp, iout 
c
c
c         normalize the input vector
c
      anorm=sqrt(sdot( n,v(1,0),1,v(1,0),1) )
      anorm=1.d0/anorm
      call sscal(n,anorm,v(1,0),1)
c         we now have a normalized first function
      if (iter.gt.1) then
c
c         form mat times the first function
c
          call honv(scr,mat,v(1,0),n)
c
c         calculate a(1)
c
          a(1)=sdot(n,v(1,0),1,scr,1)
c
c         form mat times the first function - a(1) times the first function
c         and store it in the next polynomial
c
          do 10 i=1,n
             v(i,1)=scr(i)-a(1)*v(i,0)
 10       continue
c
c         calculate b(1)
c
          b(1)=sqrt( sdot(n,v(1,1),1,v(1,1),1) )
c
c         normalize the second polynomial
c
          call sscal(n,1.d0/b(1),v(1,1),1)  
      endif
      if (iter.gt.2) then
          do 20 i=2,iter
c
c            multiply the last calculated polynomial by mat
c
             call honv(scr,mat,v(1,i-1),n)
c
c            orthogonalize it to the two previous polynomials
c            calculating a(i) as we go
c
             a(i)=sdot(n,v(1,i-1),1,scr,1)
             do 30 j=1,n
                v(j,i)=scr(j) - a(i)*v(j,i-1)
 30          continue
             do 40 j=1,n
                v(j,i) = v(j,i) - b(i-1)*v(j,i-2)
 40          continue
c
c            calculate b(i)
c
             b(i) = sqrt( sdot(n,v(1,i),1,v(1,i),1) )  
c
c            normalize the polynomial and we are done
c
             call sscal(n,1.d0/b(i),v(1,i),1)
 20       continue
      endif
c
c     diagonalize the tridiagonal matrix for the eigenvalues
c
      call rzero(vec,iter*iter)
      do 50 i=1,iter
         vec(i,i)=1.d0
 50   continue
      call copy(b,scr(2),iter-1)   
      call imtql2(iter,iter,a(1),scr,vec,ierr)
      title='eigenvalues of lanczos polynomial'
      call prntrm(title,a(1),iter,1,iter,1,iout)      
      return
      end       
