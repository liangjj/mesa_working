*deck lgrply
      subroutine lgrply(p,dp,ddp,x,y,n,m,prnt)
c***begin prologue     lgrply
c***date written       940504   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***purpose            lagrange polynomials.
c***description
c***            
c               
c               
c***references
c
c***routines called
c
c***end prologue       lgrply
c
      implicit integer (a-z)
      real*8 p, dp, ddp, x, y
      real*8 sn, ssn, fac
      logical prnt
      character*80 title 
      dimension p(m,n), dp(m,n), ddp(m,n)
      dimension x(n), y(m)
      common /io/ inp, iout

c
c     generate polynomials and derivatives with respect to x
c
      do 10 i=1,m
         zerfac = 0
         do 20 j=1,n
            fac =  y(i) - x(j) 
            if(abs(fac).le.1.d-14) then
               zerfac = j
            endif   
 20      continue
         do 30 j=1,n
            p(i,j) = 1.d0
            do 40 k=1,j-1
               p(i,j) = p(i,j)*( y(i) - x(k) )
     1                        / ( x(j) - x(k) )
 40         continue    
            do 50 k=j+1,n
               p(i,j) = p(i,j)*( y(i) - x(k) )
     1                        / ( x(j) - x(k) )
 50         continue
            if(abs(p(i,j)).gt.1.d-14) then
               sn = 0.d0
               ssn = 0.d0
               do 60 k=1,j-1
                  fac = 1.d0/( y(i) - x(k) )
                  sn = sn + fac
                  ssn = ssn + fac*fac
 60            continue
               do 70 k=j+1,n
                  fac = 1.d0/( y(i) - x(k) )
                  sn = sn + fac
                  ssn = ssn + fac*fac
 70            continue                                 
               dp(i,j) = sn*p(i,j)               
               ddp(i,j) = sn*dp(i,j) - ssn*p(i,j)
            else
               sn = 1.d0
               ssn = 0.d0
               do 80 k=1,j-1
                  if(k.ne.zerfac) then
                     fac = 1.d0/( x(j) - x(k) )
                     sn = sn*fac*( y(i) - x(k) )
                     ssn = ssn + 1.d0/(y(i) - x(k))
                  endif
 80            continue 
               do 90 k=j+1,n
                  if(k.ne.zerfac) then
                     fac = 1.d0/( x(j) - x(k) )
                     sn = sn*fac*( y(i) - x(k) )
                     ssn = ssn + 1.d0/( y(i) - x(k) )             
                  endif
 90            continue
               dp(i,j) = sn/( x(j) - x(zerfac) )
               ddp(i,j) = 2.d0*ssn*dp(i,j)
            endif                    
 30      continue
c
 10   continue
       if(prnt) then
          title='polynomials'
          call prntfm(title,p,m,n,m,n,iout)
          title='derivative of polynomials'
          call prntfm(title,dp,m,n,m,n,iout)
          title='second derivative of polynomials'
          call prntfm(title,ddp,m,n,m,n,iout)
      endif          
      return
 1    format(/,5x,'quadrature type = ',a8)
      end
