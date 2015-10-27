*deck init
c***begin prologue     init
c***date written       xxxxxx   (yymmdd)
c***revision date      890422   (yymmdd)
c***keywords           m6004, link 6004, bessel, spline
c***author             schneider, barry (lanl)
c***source             m6004
c***purpose            initialize backward recursion for bessel functions of
c***                   integer values
c***description        
c***references       
c
c***routines called
c***end prologue       init
      subroutine init(x,min,mout)
      implicit integer (a-z)
      common /io/ inp, iout
      real*8 x
      if(x.gt.1.d0) then
         mout=min
      else
         mout=-log(r1mach(1)/x)
         if(mout.gt.min) then
            write(iout,1) min, mout
            call lnkerr('quit')
         endif
	 write(iout,1) min, mout
      endif
      return
 1    format(/,1x,'min = ',i4,1x,'mout = ',i4)
      end



