*deck matpp.f
c***begin prologue     matpp
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
c***end prologue      matpp
      subroutine matpp(hpp,dmat,ddmat,f,g,rhs,work,energy,r,rlast,
     1                 ipvt,l,np,prnth,prnts,dir)
      implicit integer (a-z)
      dimension hpp(np,np), dmat(np,np), ddmat(np,np), f(np)
      dimension g(np), rhs(np), ipvt(np), r(np), work(np)
      real*8 hpp, f, g, dmat, ddmat, rhs, j, y, rlast, work
      real*8 jp, jpp, yp, ypp, energy, rr, k, r, logder
      character*80 title, dir
      logical prnth, prnts
      common /io/ inp, iout
c     fill the hamiltonian with the collocation sub-matrix and the
c     collocation right hand side.
      do 10 i=1,np
         do 20 j=1,np
            hpp(i,j)=ddmat(i,j)
 20      continue
         hpp(i,i)=hpp(i,i)+f(i)   
 10   continue         
      call copy(g,rhs,np)
c     if the calcualtion is entirely collocation base we proceed to
c     solve the equations and extract the relevant information.
c     if there is a numerove or finite difference part, just return to
c     the main part of the code and continue.     
      if (dir.eq.'collocation') then    
c     apply the boundary conditions by dropping the last equation
c     and incorporating the constraints into the other equations to
c     get an inhomogeneous system to solve.
          k=sqrt(energy)
          rr=rlast*k
          call rbes('ricatti-bessel',l,rr,j,jp,jpp,y,yp,ypp)
          logder=k*(yp/y)
          call bndryc(hpp,dmat,logder,work,np)
          if (prnth) then
              title='matrix'
              call prntrm(title,hpp,np,np,np,np,iout)
              title='rhs'
              call prntrm(title,rhs,np,1,np,1,iout)
          endif                          
          call sgefa(hpp,np,np-1,ipvt,info)
          if (info.ne.0) then
              call lnkerr('error in hamiltonian matrix')
          else          
              call sgesl(hpp,np,np-1,ipvt,rhs,0)
          endif
          call phase(rhs,r,work,energy,y,np,prnts)
      endif          
      return
      end
