*deck matqp.f
c***begin prologue     matqp
c***date written       950704   (yymmdd)
c***revision date               (yymmdd)
c***keywords           lattice, grid, matrix
c***author             schneider, barry(nsf)
c***source             
c***purpose            calculate lattice representation of 
c***                          y''  + f(x) y = g.
c***
c***
c***references
c
c***routines called
c***end prologue      matqp
      subroutine matqp(ham,r,f,g,rhs,n,np,nq,ntot)
      implicit integer (a-z)
      dimension ham(ntot,ntot), r(ntot), f(ntot), g(ntot)
      dimension rhs(ntot)
      real*8 ham, r, f, g, rhs
      common /io/ inp, iout
      if (n.eq.3) then
c     fill in the sub-diagonal elements of the three point numerov formula     
          ham(np+1,np) = 12.d+00 + stp*stp*f(np)
          rhs(np+1)    =           stp*stp*g(np)
      elseif (n.eq.5) then
c     fill in the sub-diagonal elements of the five point numerov formula    
          ham(np+1,np)   = 1920. + 688.*stp*stp*f(np)
          ham(np+1,np-1) = 465.  + 23. *stp*stp*f(np-1)
          ham(np+2,np)   = 465.  + 23. *stp*stp*f(np)
          rhs(np+1)      = 688.d0*stp*stp*g(np) + 23.d0*stp*stp*g(np-1)
          rhs(np+2)      = 23.d0*stp*stp*g(np)
      else          
          call lnkerr('error in numerov formula being used')
      endif          
      return
      end
