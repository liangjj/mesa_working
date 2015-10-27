*deck gfdiff.f
c***begin prologue     gfdiff
c***date written       950704   (yymmdd)
c***revision date               (yymmdd)
c***keywords           generalized finite difference, band
c***author             schneider, barry(nsf)
c***source             
c***purpose            calculate bands in generalized finite difference 
c***                   approximation
c***
c***
c***references
c
c***routines called
c***end prologue      gfdiff
      subroutine gfdiff(band,r,f,n,order,bw,type)
      implicit integer (a-z)
      dimension band(n,-bw:bw), r(n), f(n)
      real*8 band, r, stp, f, alpha, beta, gamma, delta, omega
      character*(*) type
      common /io/ inp, iout
      if(type.eq.'finite-difference') then
         if (order.eq.3) then
             do 10 i=2,n-1
                alpha=r(i)-r(i-1)
                beta=r(i+1)-r(i)
                band(i,-1) =   - 1.d0/(alpha*(alpha+beta))
                band(i,1)  =   - 1.d0/(beta*(alpha+beta))
                band(i,0)  =     1.d0/(alpha*beta) +f(i)
 10          continue    
         else
              call lnkerr('bad call to finite difference routine')
         endif
      elseif(type.eq.'numerov') then
         do 20 i=2,n-1
            alpha=r(i)-r(i-1)
            beta=r(i+1)-r(i)
            gamma = ( alpha*alpha + alpha*beta -beta*beta )
     1                                        /12.d0
            delta  = ( alpha*alpha + 3.d0*alpha*beta + beta*beta )
     1                                       /12.d0
            omega = alpha*( -alpha*alpha + alpha*beta +beta*beta )
     1                                        /12.d0
            band(i,-1) = beta * ( 1.d0 - gamma*f(i-1) )
            band(i,0) = - ( alpha + beta ) * ( 1.d0 + delta*f(i) )
            band(i,1) = alpha * ( 1.d0 - omega*f(i+1) )
   20     continue
      else
            call lnkerr('fdiff called with bad type')
      endif
      return
      end

