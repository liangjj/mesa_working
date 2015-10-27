*deck @(#)ovfil.f	1.1 9/8/91
c***begin prologue     ovfil
c***date written       930417   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           ovfil, link 6018
c***author             schneider, barry (nsf)  
c***source             m6018
c***purpose            fill bound-free hamiltonian and overlap matrices
c***                   
c***references        
c
c***routines called    iosys, util and mdutil
c***end prologue       ovfil
      subroutine ovfil (ovbf,hpvb,ovpb,hpb,orblst,nl,nb,nmo,maxlm,
     1                  dimmo)
      implicit integer(a-z)
      complex*16 ovbf, hpvb, ovpb, hpb
      dimension ovbf(maxlm,nmo), hpvb(maxlm,nmo)
      dimension ovpb(nl,nb), hpb(nl,nb), orblst(dimmo)
      common /io/ inp,iout
      do 10 bfn=1,nb
         orb=orblst(bfn)
         do 20 nolm1=1,nl
            ovpb(nolm1,bfn)=ovbf(nolm1,orb)
            hpb(nolm1,bfn)=hpvb(nolm1,orb)
   20    continue
   10 continue
      return
      end
