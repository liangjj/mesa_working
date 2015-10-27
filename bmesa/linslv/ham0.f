*deck ham0.f
c***begin prologue     ham0
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            kinetic energy in dvr representation.
c***                   
c***references         
c
c***routines called    
c***end prologue       ham0
      subroutine ham0(p,dp,ddp,wt,hamil,n,npts)
      implicit integer (a-z)
      real*8 p, dp, ddp, wt, hamil
      dimension p(npts,0:n-1), dp(npts,0:n-1), ddp(npts,0:n-1)
      dimension wt(npts), hamil(n,n) 
      common/io/inp, iout 
      call rzero(hamil,n*n)
      do 10 i=1,n
         do 20 j=1,i
            do 30 k=1,npts
               hamil(i,j) = hamil(i,j) - wt(k)*p(k,i-1)*ddp(k,j-1)
 30         continue
            hamil(i,j) = hamil(i,j) + p(npts,i-1)*dp(npts,j-1) 
     1                      - p(1,i-1)*dp(1,j-1)
            hamil(i,j)=.5d0*hamil(i,j)
            hamil(j,i) = hamil(i,j)
 20      continue   
 10   continue
      return
      end       
