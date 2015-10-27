*deck lpoly.f
c***begin prologue     lpoly
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           polynomials
c***author             schneider, barry (nsf)
c***source             
c***purpose            compute coefficients for orthogonal polynomials 
c***                   beginning with p  defined as:
c***                                   0
c***                   p(x) = (x-left)**alpha * (right-x)**beta
c***                    0
c***                        alpha and beta not necessarily integers                   
c***description        b p (x) = (x - a ) p   (x) - b   p   (x)
c***                    j j            j   j-1       j-1 j-2       
c***                                                          
c***references         
c
c***routines called    
c***end prologue       lpoly
      subroutine lpoly(p,x,weight,a,b,left,right,scr,n,npts,
     1                 alpha,beta)
      implicit integer (a-z)
      real*8 x, weight, p, a, b, scr, left, right, snorm, anorm
      real*8 alpha, beta
      dimension x(npts), weight(npts), p(npts,0:n-1), scr(npts)
      dimension a(0:n), b(0:n)
      common/io/inp, iout 
c     initialize.
      do 10 k=1,npts
         p(k,0) =  ( ( x(k)-left )**alpha ) * ( ( right -x(k) )**beta )
 10   continue              
      anorm=sqrt(snorm(p(1,0),p(1,0),weight,npts,.false.) )
      a(0)=1.d0/anorm
      call sscal(npts,a(0),p(1,0),1)
c     we now have a normalized first function
      if (n.gt.1) then
c         form x times the first function
          call vmul(scr,p(1,0),x,npts)
c         calculate a(1)
          a(1)=snorm(p(1,0),scr,weight,npts,.false.)
c         form x times the first function - a(1) times the first function
c         and store it in the next polynomial
          do 20 i=1,npts
             p(i,1)=scr(i)-a(1)*p(i,0)
 20       continue
c         calculate b(1)
          b(1)=sqrt( snorm(p(1,1),p(1,1),weight,npts,.false.) )
c         normalize the second polynomial
          call sscal(npts,1.d0/b(1),p(1,1),1)  
      endif
      if (n.gt.2) then
          do 30 i=2,n-1
             im1=i-1
             im2=i-2                
c            multiply the last calculated polynomial by x
             call vmul(scr,p(1,im1),x,npts)
c            orthogonalize it to the two previous polynomials
c            calculating a(i) as we go
             a(i)=snorm(p(1,im1),scr,weight,npts,.false.)
             do 40 j=1,npts
                p(j,i)=scr(j) - a(i)*p(j,im1)
 40          continue
             do 50 j=1,npts
                p(j,i) = p(j,i) - b(im1)*p(j,im2)
 50          continue
c            calculate b(i)
             b(i) = sqrt( snorm(p(1,i),p(1,i),weight,npts,.false.) )  
c            normalize the polynomial and we are done
             call sscal(npts,1.d0/b(i),p(1,i),1)
 30       continue
c         calculate the last a
          call vmul(scr,p(1,n-1),x,npts)          
          a(n)=snorm(p(1,n-1),scr,weight,npts,.false.)
      endif
      return
      end       
