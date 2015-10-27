*deck bndryc.f
c***begin prologue     bndryc
c***date written       921206   (yymmdd)
c***revision date               (yymmdd)
c***keywords           collocation, boundary condition
c***author             schneider, barry(nsf)
c***source             colloc1d
c***purpose            modify collocation matrix to account for
c***                   boundary conditions. 
c***description        the equation here is the lattice version of
c***                              y''  + f(x) y = g
c***                   subject to log-derivative boundary conditions.  
c***
c***references         
c***routines called  
c***end prologue
      subroutine bndryc(ham,dmat,logder,tvec,npt)
      implicit integer (a-z)
      real*8  ham, dmat, logder, tvec, fac
      dimension ham(npt,npt), dmat(npt,npt), tvec(npt)
      common /io/ inp, iout
c     we eliminate the value of the wavefunction at the last point using
c     the boundary condition that the log-derivative of the wavefunction
c     is known at the boundary.  this is expressed using the lattice
c     representation of the first derivative which expresses that quantity
c     in terms of the values of the function on the lattice points.
      fac=1.d0/( logder-dmat(npt,npt) )
      do 10 i=1,npt-1
         tvec(i)=fac*dmat(npt,i)
   10 continue
      do 20 i=1,npt-1
         do 30 j=1,npt-1
            ham(i,j)=ham(i,j)+ham(i,npt)*tvec(j)
   30    continue
   20 continue                           
      return
      end



