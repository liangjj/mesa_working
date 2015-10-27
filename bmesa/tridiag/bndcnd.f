*deck bndcnd.f
c***begin prologue     bndcnd
c***date written       950704   (yymmdd)
c***revision date               (yymmdd)
c***keywords           finite difference, band, boundary conditions
c***author             schneider, barry(nsf)
c***source             
c***purpose            apply boundary conditions in finite difference
c***                   approximation.
c***
c***description        the formula used is
c***                   y'(n) = 0 = ( 3.*y(n) - 4.*y(n-1) +y(n-2) )/(2.*stp)
c***                   which enables us to eliminate y(n) 
c***
c***references
c
c***routines called
c***end prologue      bndcnd
      subroutine bndcnd(band,n,order,bw,bc)
      implicit integer (a-z)
      dimension band(n,-bw:bw)
      real*8 band
      character*(*) bc
      common /io/ inp, iout
c
c     regularity at the origin is taken as the left hand boundary
c     condition.  this requires the function to be zero.
c     the right hand boundary condition is either the function or its
c     derivative is zero.
c
      if (order.eq.3) then
          if(bc.eq.'zero-function') then
             return
          elseif(bc.eq.'zero-derivative') then
             band(n-1,-1) = band(n-1,-1) - band(n-1,1)/3.d0
             band(n-1,0)  = band(n-1,0)  + band(n-1,1)*4.d0/3.d0
c              band(n-1,0) = band(n-1,0) + band(n-1,1)
          else
             call lnkerr('error in boundary condition routine')
          endif 
      elseif(order.eq.5) then
          call lnkerr('order five not yet implimented')
      else
            call lnkerr('bad call to finite difference routine')
      endif                        
      return
      end

