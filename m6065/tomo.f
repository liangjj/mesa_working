*deck @(#)tomo.f	1.1 9/8/91
c***begin prologue     tomo
c***date written       890529   (yymmdd)
c***revision date               (yymmdd)
c***keywords           transformation, kohn
c***author             schneider, barry (lanl)
c***source             m6005
c***purpose            bound-free ao to mo transformation
c***
c***references         none
c
c***routines called    ecbc
c***end prologue       tomo
 
      subroutine tomo(trans,ovlmb,ovlmbm,nolam,nkept,nmo)
      implicit integer (a-z)
      real *8 trans
      complex *16 ovlmb, ovlmbm
      dimension ovlmb(nolam,nkept), ovlmbm(nolam,nmo)
      dimension trans(nkept,nmo)
      common/io/inp,iout
      call ecbc(ovlmbm,ovlmb,trans,nolam,nkept,nmo)       
      return
      end




