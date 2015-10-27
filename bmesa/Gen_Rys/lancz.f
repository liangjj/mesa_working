*deck lancz.f
c***begin prologue     lancz
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           eigenvalues, eigenvectors
c***author             schneider, barry (nsf)
c***source
c***purpose            generate orthogonal polynomials using recursion.
c***                   using the generated alpha and beta coefficients find
c***                   the points and weights of the generalized gauss
c***                   quadrature.
c***                   
c***description        b v  = (arg - a ) v   - b    v  
c***                    j j           j   j-1   j-1  j-2       
c***                                                                       
c***                   a = <v | arg |v   >
c***                    i    j-1      j-1 
c***                  
c***                   b = <v | arg - a |v   > - <v | b    |v   >  
c***                    i    i         i  i-1      i   i-1   i-2
c***references         
c
c***routines called    
c***end prologue       lancz
      subroutine lancz(v,arg,a,b,wt,refwt,mu,scr,n,iter)
      implicit integer (a-z)
      real*8 v, arg, a, b, wt, refwt, mu, scr, anorm, scaprd
      dimension v(n,0:iter), arg(n), scr(n), wt(n), refwt(n)
      dimension a(1:iter), b(1:iter)
      common/io/inp, iout 
c
c     fill with unit vector
c
      call vfill(v(1,0),1.d0,n )
      mu=scaprd(v(1,0),v(1,0),wt,refwt,n)
      write(iout,1) mu
c
c     normalize
c
      anorm=sqrt(1.d0/mu) 
      call sscal(n,anorm,v(1,0),1)
c
c     we now have the first function
c
      if (iter.gt.1) then
c
c         form argument times the first function
c
          call vmul(scr,arg,v(1,0),n)
c
c         calculate a(1)
c
          a(1)=scaprd(v(1,0),scr,wt,refwt,n)
c
c         form mat times the first function - a(1) times the first function
c         and store it in the next polynomial
c
          call vwxs(v(1,1),scr,v(1,0),a(1),-1,n)
c
c         calculate b(1)
c
          b(1)=sqrt( scaprd(v(1,1),v(1,1),wt,refwt,n) )
c
c         normalize the second polynomial
c
          call sscal(n,1.d0/b(1),v(1,1),1)  
      endif
      if (iter.gt.2) then
          do 10 i=2,iter
c
c            multiply the last calculated polynomial by mat
c
             call vmul(scr,arg,v(1,i-1),n)
c
c            orthogonalize it to the two previous polynomials
c            calculating a(i) as we go
c
             a(i)=scaprd(v(1,i-1),scr,wt,refwt,n)
             call vwxs(v(1,i),scr,v(1,i-1),a(i),-1,n)
             call vwxs(v(1,i),v(1,i),v(1,i-2),b(i-1),-1,n)
c
c            calculate b(i)
c
             b(i) = sqrt (scaprd(v(1,i),v(1,i),wt,refwt,n ) )
c 
c            normalize the polynomial and we are done
c
             call sscal(n,1.d0/b(i),v(1,i),1)
 10       continue
c
c         calculate the last a
c
      endif
      return
 1    format(/,1x,'weight integral = ',e15.8)
      end       
