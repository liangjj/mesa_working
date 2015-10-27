*deck bndryn.f
c***begin prologue     bndryn
c***date written       921206   (yymmdd)
c***revision date               (yymmdd)
c***keywords           numerov, boundary condition
c***author             schneider, barry(nsf)
c***source             colloc1d
c***purpose            modify numerov matrix to account for
c***                   boundary conditions. 
c***description        the formula used is
c***                   y'(n) = [ { y(n) - y(n-1) }/stp +stp*y''/2. ]
c***                   we use y'(n) = logder*y(n) and the differential
c***                   equation to get:
c***                   [ 1. - stp*logder - stp*stp*f(n)/2.]*y(n) =
c***                    y(n-1) - stp*stp*g(n)*.5
c***                   which enables us to elininate y(n) 
c***
c***references         
c***routines called  
c***end prologue
      subroutine bndryn(band,rhs,f,g,logder,stp,work,n,bw,npt)
      implicit integer (a-z)
      real*8  band, rhs, f, g, logder, stp, work
      dimension band(npt,-bw:bw), rhs(npt), f(npt), g(npt), work(*)
      common /io/ inp, iout
      if (n.eq.3) then
         work(1) = band(npt-1,-1)
         work(2) = band(npt-1,0)
         work(3) = band(npt-1,1)
         work(4) = rhs(npt-1)
         work(5) = 1.d0 - logder*stp -.5d0*stp*stp*f(npt)
         work(6) = .5d0*stp*stp*g(npt)
         rhs(npt-1) = rhs(npt-1) +work(3)*work(6)/work(5)
         band(npt-1,0) = band(npt-1,0) + work(3)/work(5)
      elseif (n.eq.5) then
          call lnkerr('fifth order numerov not yet implimented')
      else
          call lnkerr('this order numerov not available')
      endif
      return
      end



