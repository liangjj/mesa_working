*deck fdiff.f
c***begin prologue     fdiff
c***date written       950704   (yymmdd)
c***revision date               (yymmdd)
c***keywords           finite difference, band
c***author             schneider, barry(nsf)
c***source             
c***purpose            calculate bands in finite difference approximation
c***
c***
c***references
c
c***routines called
c***end prologue      fdiff
      subroutine fdiff(band,f,stp,n,order,bw,type)
      implicit integer (a-z)
      dimension band(n,-bw:bw), temp(10), f(n)
      real*8 band, stp, temp, f
      character*(*) type
      common /io/ inp, iout
      if(type.eq.'standard-finite-difference') then
         if (order.eq.3) then
              temp(1)=1.d0/(stp*stp) 
              temp(2)=-2.d0*temp(1)
              temp(1)=-.5d0*temp(1)
              temp(2)=-.5d0*temp(2)
              call vfill(band(1,0),temp(2),n)
              call vfill(band(2,-1),temp(1),n-1)
              call vfill(band(1,1),temp(1),n-1)
              call vadd(band(1,0),band(1,0),f,n)
         elseif(order.eq.5) then
              temp(1)=1.d0/(12.d0*stp*stp)                  
              temp(2)=16.d0*temp(1)
              temp(3)=30.d0*temp(1)
              temp(1)=-.5d0*temp(1)
              temp(2)=-.5d0*temp(2)
              temp(3)=-.5d0*temp(3)
              call vfill(band(1,0),-temp(3),n)
              call vfill(band(1,1),temp(2),n-1)
              call vfill(band(1,2),-temp(1),n-2)
              call vfill(band(2,-1),temp(2),n-1)
              call vfill(band(3,-2),-temp(1),n-2)
              call vadd(band(1,0),band(1,0),f,n)
         else
              call lnkerr('bad call to finite difference routine')
         endif
      elseif(type.eq.'standard-numerov') then
            do 20 i=1,n-1
               band(i,0) = -( 24.d0 + 20.d0*stp*stp*f(i) )
 20         continue   
            do 30 i=2,n-1
               band(i,-1) = ( 12.d0 - 2.d0*stp*stp*f(i-1) )
 30         continue   
            do 40 i=2,n-1
               band(i,1) = (  12.d0 - 2.d0*stp*stp*f(i+1) )
 40         continue   
      else
            call lnkerr('fdiff called with bad type')
      endif
      return
      end

