*deck filham.f
c***begin prologue     filham
c***date written       921206   (yymmdd)
c***revision date               (yymmdd)
c***keywords           collocation, hamiltonian
c***author             schneider, barry(nsf)
c***source             colloc1d
c***purpose            fill total hamiltonian with sub-block
c***                    
c***
c***references         
c***routines called  
c***end prologue
      subroutine filham(ham,hampp,n,np)
      implicit integer (a-z)
      real*8  ham, hampp
      dimension ham(n,n), hampp(np,np)
      common /io/ inp, iout
      do 10 i=1,np
         do 20 j=1,np
            ham(i,j)=hampp(i,j)
   20    continue
   10 continue                                
      return
      end



