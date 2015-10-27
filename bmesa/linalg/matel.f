c $Header$
*deck matel.f
c***begin prologue     matel
c***date written       930201   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           matel, link 6201, wavefunction
c***author             schneider, barry (nsf)
c***source             m6203
c***purpose            matrix elements of ( del**2 + k**2 )
c*** 
c***description        matrix elements of the wave operator for the
c***                   free and bound orbitals for a given physical
c***                   channel are computed.
c***                   
c***                   
c***                    
c***references         
c***routines called
c***end prologue       matel
      subroutine matel (psln,rhsfn,bndfn,plm,phim,bndfn)
      implicit integer(a-z)
      real*8 psln, rhsfn, plm, phim, wph, wthet, wr
      common /io/ inp, iout
      dimension psln(nr,nl), rhsfn(nr,nl), plm(nth,nl), phim(nph)

      return
      end


