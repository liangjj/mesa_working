*deck cpoly.f
c***begin prologue     cpoly
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           polynomials
c***author             schneider, barry (nsf)
c***source             
c***purpose            compute coefficients for orthogonal polynomials.
c***                   
c***description        b p (x) = (x - a ) p   (x) - b   p   (x)
c***                    j j            j   j-1       j-1 j-2       
c***                                                          
c***references         
c
c***routines called    
c***end prologue       cpoly
      subroutine cpoly(p,x,weight,a,b,left,right,scr,n,npts,
     1                 pleft,pright)
      implicit integer (a-z)
      real*8 x, weight, p, a, b, scr, left, right, snorm, anorm
      dimension x(npts), weight(npts), p(npts,0:n-1), scr(npts)
      dimension a(0:n), b(0:n)
      common/io/inp, iout 
c     first polynomial does not have to be one.
      call poly0(p(1,0),x,left,right,pleft,pright,0,npts)      
      anorm=sqrt(snorm(p(1,0),p(1,0),weight,npts,.false.) )
      a(0)=1.d0/anorm
      write(iout,*) a(0)
      call sscal(npts,a(0),p(1,0),1)
c     we now have a normalized first function
      if (n.gt.1) then
c         form x times the first function
          call vmul(scr,p(1,0),x,npts)
c         calculate a(1)
          a(1)=snorm(p(1,0),scr,weight,npts,.false.)
c         form x times the first function - a(1) times the first function
c         and store it in the next polynomial
          do 10 i=1,npts
             p(i,1)=scr(i)-a(1)*p(i,0)
 10       continue
c         calculate b(1)
          b(1)=sqrt( snorm(p(1,1),p(1,1),weight,npts,.false.) )
c         normalize the second polynomial
          call sscal(npts,1.d0/b(1),p(1,1),1)  
      endif
      if (n.gt.2) then
          do 20 i=2,n-1
             im1=i-1
             im2=i-2                
c            multiply the last calculated polynomial by x
             call vmul(scr,p(1,im1),x,npts)
c            orthogonalize it to the two previous polynomials
c            calculating a(i) as we go
             a(i)=snorm(p(1,im1),scr,weight,npts,.false.)
             do 30 j=1,npts
                p(j,i)=scr(j) - a(i)*p(j,im1)
 30          continue
             do 40 j=1,npts
                p(j,i) = p(j,i) - b(im1)*p(j,im2)
 40          continue
c            calculate b(i)
             b(i) = sqrt( snorm(p(1,i),p(1,i),weight,npts,.false.) )  
c            normalize the polynomial and we are done
             call sscal(npts,1.d0/b(i),p(1,i),1)
 20       continue
c         calculate the last a
          call vmul(scr,p(1,n-1),x,npts)          
          a(n)=snorm(p(1,n-1),scr,weight,npts,.false.)
      endif
      return
      end       
