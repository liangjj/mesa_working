*deck matqq.f
c***begin prologue     matqq
c***date written       950704   (yymmdd)
c***revision date               (yymmdd)
c***keywords           lattice, grid, matrix
c***author             schneider, barry(nsf)
c***source             
c***purpose            calculate numerov part of lattice representation of 
c***                          y''  + f(x) y = g.
c***
c***
c***references
c
c***routines called
c***end prologue      matqq
      subroutine matqq(hamqq,r,f,g,rhs,bandy,bandg,work,energy,rlast,l,
     1                 n,np,nq,ntot,order,bw,fixed,prnt)
      implicit integer (a-z)
      dimension hamqq(ntot,nq), r(0:nq), f(0:nq), g(0:nq)
      dimension rhs(1:nq), bandy(1:nq,-bw:bw), bandg(1:nq,-bw:bw)
      dimension work(nq)
      real*8 hamqq, r, f, g, rhs, bandy, bandg, work
      real*8 j, y, rlast, jp, jpp, yp, ypp, energy, rr, k, logder
      logical fixed, prnt
      character*80 title
      common /io/ inp, iout
      if (n.eq.3) then
          call num3pt(bandy,bandg,r(1),f(1),nq,fixed)
          rhs(1) = rhs(1) + bandg(1,0)*g(1) + 
     1                              bandg(1,1)*g(2)
          do 40 i=2,nq-1
             rhs(i) = ( bandg(i,-1)*g(i-1) + bandg(i,0)*g(i) + 
     1                                bandg(i,1)*g(i+1) )
 40       continue
c         we assume that the inhomogeneity vanishes outside the last
c         point.                      
          rhs(nq) = ( bandg(nq,-1)*g(nq-1) + bandg(nq,0)*g(nq) )
c     apply the boundary conditions by expressing the value of the
c     wavefunction at the last point in terms of known quantities
c     and the wavefunction at previous points.
          k=sqrt(energy)
          rr=rlast*k
          call rbes('ricatti-bessel',l,rr,j,jp,jpp,y,yp,ypp)
          logder=k*(yp/y)
          call bndryn(bandy,rhs,f(1),g(1),logder,stp,work,order,bw,
     1                nq,nq)
          work(5)=y
          do 10 i=1,nq
             hamqq(i,i)=bandy(i,0)
 10       continue
          do 20 i=1,nq-1
             hamqq(i,i+1)=bandy(i,1)
 20       continue
          do 30 i=2,nq
             hamqq(i,i-1)=bandy(i,-1)
 30       continue
          if (prnt) then
              title='q-space rhs'
              call prntrm(title,rhs,nq,1,nq,1,iout)
              title='q-space hamiltonian'
              call prntrm(title,hamqq,nq,nq,ntot,nq,iout)
          endif
      elseif(n.eq.5) then
          call num5pt(bandy,bandg,f,stp,nq)
          do 50 i=1,nq
             hamqq(i,i)=bandy(i,0)
 50       continue
          do 60 i=1,nq-1
             hamqq(i,i+1)=bandy(i,1)
 60       continue
          do 70 i=1,nq-2
             hamqq(i,i+2)=bandy(i,2)
 70       continue
          do 80 i=2,nq
             hamqq(i,i-1)=bandy(i,-1)
 80       continue
          do 90 i=3,nq
             hamqq(i,i-2)=bandy(i,-2)
 90       continue
          rhs(1) = rhs(1) + 
     1             stp*stp*( 2358.d0*g(1) +688.d0*g(2) +23.d0*g(3) )
          rhs(2) = rhs(2) + 
     1             stp*stp*( 2358.d0*g(2) +688.d0*g(1) +688.d0*g(3) +
     2                       23.d0*g(4) )
          do 100 i=3,nq-2
             rhs(i) = 23.d0*stp*stp*( g(i-2) + g(i+2) ) +
     1                688.*stp*stp*( g(i+1) + g(i-1) )  +
     2                2358.*stp*stp*g(i)  
 100     continue            
         rhs(nq-1) = stp*stp*( 23.d0*g(nq-3) + 688.d0*g(nq-2) +
     1                         2358.d0*g(nq-1) ) 
         rhs(nq) = stp*stp*( 23.d0*g(nq-2) + 688.d0*g(nq-1) +  
     1                       2358.d0*g(nq) ) 
      else
          call lnkerr('quit')
      endif    
      return
      end
