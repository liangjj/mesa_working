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
      subroutine matqq(hqq,f,g,rhs,band,stp,n,m,nq,bw)
      implicit integer (a-z)
      dimension hqq(m,nq), f(nq), g(nq), rhs(nq), band(nq,-bw:bw)
      real*8 hqq, f, g, rhs, band, stp
      common /io/ inp, iout
      if (n.eq.3) then
          call num3pt(band,f,stp,nq)
          do 10 i=1,nq
             hqq(i,i)=band(i,0)
 10       continue
          do 20 i=1,nq-1
             hqq(i,i+1)=band(i,1)
 20       continue
          do 30 i=2,nq
             hqq(i,i-1)=band(i,-1)
 30       continue
c         the g(0) part of rhs(1) has already been added into rhs(1)
c         in the subroutine matqp
          rhs(1)=rhs(1)+stp*stp*( 10.d0*g(1) + g(2) )
          do 40 i=2,nq-1
             rhs(i) = stp*stp*( g(i-1) + 10.d0*g(i) +g(i+1) )
 40       continue
c         we assume that the inhomogeneity vanishes outside the last
c         point.                      
          rhs(nq)= stp*stp*( g(nq-1) + 10.d0*g(nq) )
      elseif(n.eq.5) then
          call num5pt(band,f,stp,nq)
          do 50 i=1,nq
             hqq(i,i)=band(i,0)
 50       continue
          do 60 i=1,nq-1
             hqq(i,i+1)=band(i,1)
 60       continue
          do 70 i=1,nq-2
             hqq(i,i+2)=band(i,2)
 70       continue
          do 80 i=2,nq
             hqq(i,i-1)=band(i,-1)
 80       continue
          do 90 i=3,nq
             hqq(i,i-2)=band(i,-2)
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
