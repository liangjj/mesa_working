*deck matpq.f
c***begin prologue     matpq
c***date written       950704   (yymmdd)
c***revision date               (yymmdd)
c***keywords           lattice, grid, matrix
c***author             schneider, barry(nsf)
c***source             
c***purpose            to replace the last collocation point equation
c***                   with one enforcing continuity of the derivative
c***                   at the last point.  this is needed when 
c***                   collocation and numerov are combined in order to
c***                   get a well defined set of equations.
c***
c***description        the formula y'(n) = ( y(n+1) - y(n) )/stp
c***                                                -.5*stp*y''(n)
c***                   is set equal to the lattice form
c***                               y'(n) = sum d(n,m) * y(m)
c***                                    and
c***                               y''  + f(x) y = g. 
c***                   is used to replace the second derivative.
c***                   the last equation is dropped and the matrix 
c***                   appropriately modified.
c***
c***
c***references
c
c***routines called
c***end prologue      matpq
      subroutine matpq(hampp,hampq,dmat,r,fp,gp,rhsp,np,nq,n,prnt)
      implicit integer (a-z)
      dimension hampp(n,np), hampq(n,nq), dmat(np,np), r(0:n) 
      dimension fp(*), gp(*), rhsp(*)
      real*8 hampp, hampq, dmat, r, fp, gp, rhsp, stp
      logical prnt
      character*80 title
      common /io/ inp, iout
      do 10 i=1,np
         hampp(np,i)=dmat(np,i)
 10   continue
      stp=r(1)-r(0)
      hampp(np,np) = hampp(np,np) + ( 1.d0/stp - .5d0*stp*fp(np) )
      rhsp(np) = -.5d0*stp*gp(np)
      hampq(np,1) = -1.d0/stp
      if (prnt) then
          title='p-space hamiltonian'
          call prntrm(title,hampp,np,np,n,n,iout)
          title='pq-coupling block'
          call prntrm(title,hampq,np,nq,n,n,iout)
      endif   
      return
      end
