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
     1                 ipvt,l,np,ntot,prnth,prnts,dir,nobc)
      implicit integer (a-z)
      dimension hpp(ntot,np), dmat(np,np), ddmat(np,np), f(np)
      dimension g(np), rhs(np), ipvt(np), r(np), work(np)
      real*8 hpp, f, g, dmat, ddmat, rhs, j, y, rlast, work
      real*8 jp, jpp, yp, ypp, energy, rr, k, r, logder
      character*80 title, dir
      logical prnth, prnts, nobc
      common /io/ inp, iout
c     fill the hamiltonian with the collocation sub-matrix and the
c     collocation right hand side.
      do 10 i=1,np
         do 20 j=1,np
            hpp(i,j)=ddmat(i,j)
 20      continue
         hpp(i,i)=hpp(i,i)+f(i)   
 10   continue
      if (prnth) then
          title='matrix before boundary conditions applied'
          call prntrm(title,hpp,np,np,ntot,ntot,iout)
      endif                   
      call copy(g,rhs,np)
c     if the calculation is entirely collocation based we proceed to
c     solve the equations and extract the relevant information.
c     if there is a numerov or finite difference part, just return to
c     the main part of the code and continue.     
      if (dir.eq.'collocation') then    
c     apply the boundary conditions by dropping the last equation
c     and incorporating the constraints into the other equations to
c     get a modified inhomogeneous system to solve.  this modified system
c     uses the known logarithmic derivative to remove the equation for the
c     function at the last point with the logarithmic derivative and the
c     values of the functions at the other points.  see subroutine bndryc
c     for details.
          if(.not.nobc) then
             k=sqrt(energy)
             rr=rlast*k
             call rbes('ricatti-bessel',l,rr,j,jp,jpp,y,yp,ypp)
             logder=k*(yp/y)
             call bndryc(hpp,dmat,logder,work,np)
             neq=np-1
          else
             k=sqrt(energy)
             rr=rlast*k
             call rbes('ricatti-bessel',l,rr,j,jp,jpp,y,yp,ypp)   
             neq=np
          endif     
          if (prnth) then
              title='matrix'
              call prntrm(title,hpp,np,np,np,np,iout)
              title='rhs'
              call prntrm(title,rhs,np,1,np,1,iout)
          endif                          
          call sgefa(hpp,np,neq,ipvt,info)
          if (info.ne.0) then
              call lnkerr('error in hamiltonian matrix')
          else          
              call sgesl(hpp,np,neq,ipvt,rhs,0)
          endif
          call phase(rhs,r,work,energy,y,np,prnts,nobc)
      endif          
      return
      end
