*deck lgngr2
      subroutine lgngr2(p,dp,ddp,x,y,nx,ny,prnt,drctv)
c***begin prologue     lgrply
c***date written       940504   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***purpose            lagrange polynomials of a quadratic variable
c***                   at arbitrary points.
c***description
c***            
c               
c               
c***references
c
c***routines called
c
c***end prologue       lgngr2
c
      implicit integer (a-z)
      real*8 p, dp, ddp
      real*8 x, y
      real*8 sn, ssn, fac
      logical prnt
      character*80 title 
      character*(*) drctv
      dimension p(ny,nx), dp(ny,nx), ddp(ny,nx)
      dimension x(nx), y(ny)
      common /io/ inp, iout
c
c     generate polynomials and derivatives with respect to x
c
      do 10 i=1,ny
         zerfac = 0
         do 20 j=1,nx
            fac =  y(i)*y(i) - x(j)*x(j) 
            if(abs(fac).le.1.d-10) then
               zerfac = j
            endif   
 20      continue
         do 30 j=1,nx
            p(i,j) = 1.d0
            do 40 k=1,j-1
               p(i,j) = p(i,j)*( y(i)*y(i) - x(k)*x(k) )
     1                        / ( x(j)*x(j) - x(k)*x(k) )
 40         continue    
            do 50 k=j+1,nx
               p(i,j) = p(i,j)*( y(i)*y(i) - x(k)*x(k) )
     1                        / ( x(j)*x(j) - x(k)*x(k) )
 50         continue
            if(drctv.ne.'functions-only') then
               if(abs(p(i,j)).gt.1.d-10) then
                  sn = 0.d0
                  ssn = 0.d0
                  do 60 k=1,j-1
                     fac = 1.d0/( y(i)*y(i) - x(k)*x(k) )
                     sn = sn + fac
                     ssn = ssn + fac*fac
 60               continue
                  do 70 k=j+1,nx
                     fac = 1.d0/( y(i)*y(i) - x(k)*x(k) )
                     sn = sn + fac
                     ssn = ssn + fac*fac
 70               continue                                 
                  dp(i,j) = 2.d0*y(i)*sn*p(i,j)               
                  ddp(i,j) = 2.d0*y(i)*sn*dp(i,j) 
     1                                + 
     2                       2.d0 * ( sn - y(i)*ssn ) *p(i,j)
               else
                  first=j
                  second=zerfac
                  if(first.gt.second) then
                     first=zerfac
                     second=j
                  endif
                  sn = 1.d0
                  ssn = 0.d0
                  do 80 k=1,first-1
                     fac = 1.d0/( x(j)*x(j) - x(k)*x(k) )
                     sn = sn*fac*( y(i)*y(i) - x(k)*x(k) )
                     ssn = ssn + 1.d0/(y(i)*y(i) - x(k)*x(k))
 80               continue 
                  do 90 k=first+1,second-1
                     fac = 1.d0/( x(j)*x(j) - x(k)*x(k) )
                     sn = sn*fac*( y(i)*y(i) - x(k)*x(k) )
                     ssn = ssn + 1.d0/(y(i)*y(i) - x(k)*x(k))
 90               continue
                  do 100 k=second+1,nx
                     fac = 1.d0/( x(j)*x(j) - x(k)*x(k) )
                     sn = sn*fac*( y(i)*y(i) - x(k)*x(k) )
                     ssn = ssn + 1.d0/(y(i)*y(i) - x(k)*x(k))
 100              continue
                  dp(i,j) = 2.d0*y(i)*sn/( x(j)*x(j) 
     1                                         - 
     2                                     x(zerfac)*x(zerfac) )
                  ddp(i,j) = 2.d0*sn/( x(j)*x(j) 
     1                                     - 
     2                                 x(zerfac)*x(zerfac) )
                  ddp(i,j) = ddp(i,j) + 2.d0*y(i)*ssn*dp(i,j)
               endif                    
            endif
 30      continue
c
 10   continue
      if(prnt) then
          title='quadratic polynomials'
          call prntfm(title,p,ny,nx,ny,nx,iout)
          if(drctv.ne.'functions-only') then
             title='derivative of quadratic polynomials'
             call prntfm(title,dp,ny,nx,ny,nx,iout)
             title='second derivative of quadratic polynomials'
             call prntfm(title,ddp,ny,nx,ny,nx,iout)
          endif           
      endif
      return
      end















