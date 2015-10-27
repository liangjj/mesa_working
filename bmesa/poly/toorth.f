*deck toorth.f
c***begin prologue     toorth
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           polynomials
c***author             schneider, barry (nsf)
c***source             
c***purpose            convert from polynomials to orthogonal functions
c***                   
c***                   
c***                                                          
c***references         
c
c***routines called    
c***end prologue       toorth
      subroutine toorth(p,dp,ddp,x,wt,dwt,ddwt,type,n,npts,noder)
      implicit integer (a-z)
      real*8 p, dp, ddp, x, wt, dwt, ddwt
      character*(*) type
      logical noder
      dimension p(npts,0:n-1), dp(npts,0:n-1), ddp(npts,0:n-1)
      dimension x(npts), wt(npts), dwt(npts), ddwt(npts)
      common/io/inp, iout
      if(noder) then 
         if(type.eq.'laguerre') then    
            do 10 i = 1, npts
               wt(i) = exp(-.5d0*x(i))
 10         continue 
         elseif(type.eq.'hermite') then
            do 20 i = 1, npts
               wt(i) = exp(-.5d0*x(i)*x(i))
 20         continue 
         elseif(type.eq.'one') then
             call vfill(wt,1.d0,npts)
         else          
             call lnkerr('error in weight function')   
         endif
         call vmmul(wt,p,p,npts,n)
      else
         if(type.eq.'laguerre') then    
            do 30 i = 1, npts
               wt(i) = exp(-.5d0*x(i))
 30         continue
            call vscale(dwt,wt,-.5d0,npts)
            call vscale(ddwt,wt,.25d0,npts)   
         elseif(type.eq.'hermite') then
            do 40 i = 1, npts
               wt(i) = exp(-.5d0*x(i)*x(i))
               dwt(i) = -x(i)*wt(i)
               ddwt(i) = -wt(i) - x(i)*dwt(i)
 40         continue
         endif
         if(type.ne.'one') then         
            do 50 i = 0, n-1
               do 60 j = 1, npts
c                the order of ddp before dp before p matters
                  ddp(j,i) = ddp(j,i)*wt(j) + 2.d0*dwt(j)*dp(j,i)
     1                                      +
     2                                    ddwt(j)*p(j,i)
                  dp(j,i) = dp(j,i)*wt(j) + p(j,i)*dwt(j)
                  p(j,i) = wt(j)*p(j,i)
   60          continue
   50       continue
         endif
      endif                                                                          
      return
      end       
