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
c***references         
c
c***routines called    
c***end prologue       lancz
      subroutine lancz(v,mat,a,b,scr,n,iter,prn)
      implicit integer (a-z)
      real*8 v, mat, a, b, scr, sdot, anorm
      logical prn
      character*80 title
      dimension v(n,0:iter-1), mat(n,n), scr(n)
      dimension a(0:iter), b(0:iter)
      common/io/inp, iout
      call rzero(v(1,0),n)
      v(n,0)=1.d0       
      anorm=sqrt(sdot( n,v(1,0),1,v(1,0),1) )
      a(0)=1.d0/anorm
      call sscal(n,a(0),v(1,0),1)
c     we now have a normalized first function
      if (iter.gt.1) then
c         form mat times the first function
          call honv(v(1,0),scr,mat,n)
c         calculate a(1)
          a(1)=sdot(n,v(1,0),1,scr,1)
c         form mat times the first function - a(1) times the first function
c         and store it in the next polynomial
          do 10 i=1,n
             v(i,1)=scr(i)-a(1)*v(i,0)
 10       continue
c         calculate b(1)
          b(1)=sqrt( sdot(n,v(1,1),1,v(1,1),1) )
c         normalize the second polynomial
          call sscal(n,1.d0/b(1),v(1,1),1)  
      endif
      if (iter.gt.2) then
          do 20 i=2,iter-1
             im1=i-1
             im2=i-2                
c            multiply the last calculated polynomial by mat
             call honv(v(1,im1),scr,mat,n)
c            orthogonalize it to the two previous polynomials
c            calculating a(i) as we go
             a(i)=sdot(n,v(1,im1),1,scr,1)
             do 30 j=1,n
                v(j,i)=scr(j) - a(i)*v(j,im1)
 30          continue
             do 40 j=1,n
                v(j,i) = v(j,i) - b(im1)*v(j,im2)
 40          continue
c            calculate b(i)
             b(i) = sqrt( sdot(n,v(1,i),1,v(1,i),1) )  
c            normalize the polynomial and we are done
             call sscal(n,1.d0/b(i),v(1,i),1)
 20       continue
c         calculate the last a
          call honv(v(1,iter-1),scr,mat,n)          
          a(iter)=sdot(n,v(1,iter-1),1,scr,1)
      endif
      if (prn) then
          title='a coefficients'
          call prntrm(title,a(1),iter,1,iter,1,iout)
          title='b coefficients'
          call prntrm(title,b(1),iter-1,1,iter-1,1,iout)
          do 100 i=1,iter
             do 200 j=1,i
               anorm=sdot(n,v(1,i-1),1,v(1,j-1),1)
               write(iout,*) 'overlap =',anorm
 200         continue
 100      continue   
      endif
c     diagonalize the tridiagonal matrix for the eigenvalues
      call imtql1(iter,a(1),b(0),ierr)
      write(iout,*) 'ierr =',ierr
      return
      end       
