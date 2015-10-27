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
c***                   which enables us to eliminate y(n) 
c***
c***references         
c***routines called  
c***end prologue
      subroutine bndryn(bandy,rhs,f,g,logder,stp,work,n,bw,ntot,nq)
      implicit integer (a-z)
      real*8  bandy, rhs, f, g, logder, stp, work
      dimension bandy(nq,-bw:bw), rhs(ntot), f(ntot), g(ntot), work(*)
      common /io/ inp, iout
      if (n.eq.3) then
         work(1) = bandy(ntot-1,-1)
         work(2) = bandy(ntot-1,0)
         work(3) = bandy(ntot-1,1)
         work(4) = rhs(ntot-1)
         work(5) = 1.d0 - logder*stp -.5d0*stp*stp*f(ntot)
         work(6) = .5d0*stp*stp*g(ntot)
         rhs(ntot-1) = rhs(ntot-1) +work(3)*work(6)/work(5)
         bandy(ntot-1,0) = bandy(ntot-1,0) + work(3)/work(5)
      elseif (n.eq.5) then
          call lnkerr('fifth order numerov not yet implimented')
      else
          call lnkerr('this order numerov not available')
      endif
      return
      end
