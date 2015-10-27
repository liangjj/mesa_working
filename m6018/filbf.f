*deck @(#)filbf.f	1.1 9/8/91
c***begin prologue     filbf
c***date written       930417   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           filbf, link 6018
c***author             schneider, barry (nsf)  
c***source             m6018
c***purpose            fill bound-free hamiltonian and overlap matrices
c***                   
c***references        
c
c***routines called    iosys, util and mdutil
c***end prologue       filbf
      subroutine filbf (ovbf,hpvb,opb,omb,hpb,hmb,nl,nb)
      implicit integer(a-z)
      complex*16 ovbf, hpvb, opb, hpb
      real*8 omb, hmb
      dimension ovbf(nl,nb), hpvb(nl,nb), omb(nl,nb), hmb(nl,nb)
      dimension opb(nl,nb), hpb(nl,nb)
      common /io/ inp,iout
      do 10 bfn=1,nb
         do 20 nolm1=1,nl
            opb(nolm1,bfn)=ovbf(nolm1,bfn)
            omb(nolm1,bfn)=imag(opb(nolm1,bfn))
            hpb(nolm1,bfn)=hpvb(nolm1,bfn)
            hmb(nolm1,bfn)=imag(hpb(nolm1,bfn))
   20    continue
   10 continue
      return
      end
