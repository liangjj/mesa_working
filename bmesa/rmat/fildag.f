*deck fildag
c***begin prologue     fildag
c***date written       940213   (yymmdd)
c***revision date               (yymmdd)
c***keywords           hamiltonian, matrix
c***author             schneider, barry (nsf)
c***source             
c***purpose            fill diagonal blocks of the hamiltonian matrix.
c***description        four types of channel blocks can appear
c***                   corresponding to the closed-closed, closed-open,
c***                   open-closed and open-open portions of the hamiltonian
c***                   matrix. the hamiltonian matrix is passed for all four
c***                   parts of this block in the calling routine.  the t
c***                   matrix for the block is passed in its entirety.  it is
c***                   assumed that the t matrix is blocked in the same fashion
c***                   with all closed orbitals first in the list.
c***                   
c
c***routines called
c***end prologue       fildag
      subroutine fildag(h,t,ec,ni,n)
      implicit integer (a-z)
      common /io/ inp, iout
      real*8 h, t, ec
      dimension h(n,n), t(ni,ni)
      do 10 i=1,ni
         h(i,i)=h(i,i)+ec
         do 20 j=1,i
            h(i,j)=h(i,j)+t(i,j)
            h(j,i)=h(i,j)
   20    continue                                                     
   10 continue
      return
      end
