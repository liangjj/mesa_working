*deck lgpoly.f
c***begin prologue     lgpoly
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           polynomials
c***author             schneider, barry (nsf)
c***source             
c***purpose            compute orthogonal polynomials and their first
c***                   and second derivatives using the recursion relation;
c***                   
c***description        b p (x) = (x - a ) p   (x) - b   p   (x)
c***                    j j            j   j-1       j-1 j-2
c***    
c***                   beginning with p  defined as:
c***                                   0
c***                   p(x) = (x-left)**alpha * (right-x)**beta
c***                    0
c***                        alpha and beta not necessarily integers     
c***references         
c
c***routines called    
c***end prologue       lgpoly
      subroutine lgpoly(p,dp,ddp,x,a,b,left,right,alpha,beta,
     1                  n,npts,noder)
      implicit integer (a-z)
      real*8 p, dp, ddp, x, a, b, left, right, alpha, beta
      real*8 tmp1, tmp2, tmp3, tmp4
      logical noder
      dimension p(npts,0:n-1), dp(npts,0:n-1), ddp(npts,0:n-1)
      dimension x(npts), a(0:n), b(0:n)
      common/io/inp, iout 
c     initialize
      do 10 k=1,npts
         p(k,0) =  ( ( x(k)-left )**alpha ) * ( ( right -x(k) )**beta )
 10   continue              
      call sscal(npts,a(0),p(1,0),1)         
      do 20 i=1,npts
         p(i,1) = (  x(i)-a(1) )*p(i,0)/b(1)
 20      continue
      if(n.ge.2) then
         do 30 i=2,n-1
            do 40 j=1,npts
               p(j,i) = ( ( x(j) -a(i) )*p(j,i-1) 
     1                                - 
     2                                  b(i-1)*p(j,i-2) )/b(i)
 40         continue
 30      continue   
      endif
      if(.not.noder) then
         do 50 i=1,npts
            tmp1=x(i)-left
            tmp2=right -x(i)
            tmp3=tmp1**alpha 
            tmp4=tmp2**beta
            dp(i,0) =  alpha*tmp3*tmp4/tmp1 - beta*tmp3*tmp4/tmp2
 50      continue                        
         call sscal(npts,a(0),dp(1,0),1)
         do 60 i=1,npts
            dp(i,1) = ( ( x(i) - a(1) )*dp(i,0) + p(i,0) )/b(1)
 60      continue
         do 70 i=1,npts
            tmp1=x(i)-left
            tmp2=right -x(i)
            tmp3=tmp1**alpha 
            tmp4=tmp2**beta
            ddp(i,0) =  alpha*(alpha-1)*tmp3*tmp4/(tmp1*tmp1)
     1                               - 
     2                  2.d0*alpha*beta*tmp3*tmp4/(tmp1*tmp2)
     3                               +
     4                  beta*(beta-1)*tmp3*tmp4/(tmp2*tmp2)
 70      continue                   
         call sscal(npts,a(0),ddp(1,0),1)
         do 80 i=1,npts
            ddp(i,1)=( ( x(i) - a(1) )*ddp(i,0) + 2.d0*dp(i,0) )/b(1)                
 80      continue
         if(n.ge.2) then
            do 90 i=2,n-1
               do 100 j=1,npts
                  dp(j,i) = ( ( x(j) - a(i) )*dp(j,i-1) 
     1                               - 
     2                         b(i-1)*dp(j,i-2) + p(j,i-1) )/b(i)
                  ddp(j,i) = ( ( x(j) - a(i) )*ddp(j,i-1) 
     1                               - 
     2                          b(i-1)*ddp(j,i-2) 
     3                               + 
     4                          2.d0*dp(j,i-1) )/b(i)
 100           continue   
 90         continue   
         endif
      endif         
      return
      end       
