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
      subroutine matqp(hamqp,r,f,g,rhs,n,np,nq,ntot,prnt)
      implicit integer (a-z)
      dimension hamqp(ntot,np), r(0:ntot), f(0:ntot), g(0:ntot)
      dimension rhs(1:ntot)
      real*8 hamqp, r, f, g, rhs, alpha, beta, b
      logical prnt
      character*80 title
      common /io/ inp, iout
c     the r, f and g arays begin at 0 rather than one because the first
c     point in q-space was defined in the main routine to be the last point
c     in p-space and we need to access these arrays to define the matrix
c     elements      
      if (n.eq.3) then
c     fill in the sub-diagonal elements of the three point numerov formula
          alpha = r(1)-r(0)  
          beta  = r(2)-r(1)
          b = beta * ( alpha*alpha + alpha*beta - beta*beta ) /12.d0
          hamqp(1,np) =  beta + b*f(0)     
          rhs(1)    =  b*g(0)
          if(prnt) then
             title='matrix coupling q to p-space'
             call prntrm(title,hamqp,nq,np,ntot,np,iout)
          endif               
      elseif (n.eq.5) then
c     fill in the sub-diagonal elements of the five point numerov formula    
c          ham(np+1,np)   = 1920. + 688.*stp*stp*f(np)
c          ham(np+1,np-1) = 465.  + 23. *stp*stp*f(np-1)
c          ham(np+2,np)   = 465.  + 23. *stp*stp*f(np)
c          rhs(np+1)      = 688.d0*stp*stp*g(np) + 23.d0*stp*stp*g(np-1)
c          rhs(np+2)      = 23.d0*stp*stp*g(np)
      else          
          call lnkerr('error in numerov formula being used')
      endif          
      return
      end
