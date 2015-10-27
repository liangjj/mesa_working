*deck gpoly.f
c***begin prologue     gpoly
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
c***references         
c
c***routines called    
c***end prologue       gpoly
      subroutine gpoly(p,dp,ddp,x,a,b,left,right,nleft,nright,
     1                 n,npts,noder)
      implicit integer (a-z)
      real*8 p, dp, ddp, x, a, b, left, right
      logical noder
      dimension p(npts,0:n-1), dp(npts,0:n-1), ddp(npts,0:n-1)
      dimension x(npts), a(0:n), b(0:n)
      common/io/inp, iout 
c     first polynomial does not have to be one.
c     initialize the first two members of the polynomial set.
      call poly0(p(1,0),x,left,right,nleft,nright,0,npts)
      call sscal(npts,a(0),p(1,0),1)         
      do 10 i=1,npts
         p(i,1) = (  x(i)-a(1) )*p(i,0)/b(1)
 10   continue
      if(n.ge.2) then
         do 20 i=2,n-1
            do 30 j=1,npts
               p(j,i) = ( ( x(j) -a(i) )*p(j,i-1) 
     1                                - 
     2                                  b(i-1)*p(j,i-2) )/b(i)
 30         continue
 20      continue   
      endif
      if(.not.noder) then
         call poly0(dp(1,0),x,left,right,nleft,nright,1,npts)
         call sscal(npts,a(0),dp(1,0),1)
         do 40 i=1,npts
            dp(i,1) = ( ( x(i) - a(1) )*dp(i,0) + p(i,0) )/b(1)
 40      continue
         call poly0(ddp(1,0),x,left,right,nleft,nright,2,npts)            
         call sscal(npts,a(0),ddp(1,0),1)
         do 50 i=1,npts
            ddp(i,1)=( ( x(i) - a(1) )*ddp(i,0) + 2.d0*dp(i,0) )/b(1)
 50      continue
         if(n.ge.2) then
            do 60 i=2,n-1
               do 70 j=1,npts
                  dp(j,i) = ( ( x(j) - a(i) )*dp(j,i-1) 
     1                               - 
     2                         b(i-1)*dp(j,i-2) + p(j,i-1) )/b(i)
                  ddp(j,i) = ( ( x(j) - a(i) )*ddp(j,i-1) 
     1                               - 
     2                          b(i-1)*ddp(j,i-2) 
     3                               + 
     4                          2.d0*dp(j,i-1) )/b(i)
 70            continue   
 60         continue   
         endif
      endif         
      return
      end       
