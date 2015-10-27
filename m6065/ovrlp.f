*deck @(#)ovrlp.f
c***begin prologue     ovrlp
c***date written       xxxxxx   (yymmdd)
c***revision date      920409   (yymmdd)
c***keywords           kohn integrals
c***author             schneider, barry (nsf)
c***source             m6060
c***purpose            overlaps of bound functions with vlamdas
c***references
c
c***routines called    iosys, util and mdutil
c***end prologue       ovrlp
      subroutine ovrlp(ovlmb,basis,vlamda,nao,npt,nolam)
      implicit integer (a-z)
      real *8 basis
      complex *16 vlamda, ovlmb
      dimension vlamda(npt,nolam), ovlmb(nolam,nao), basis(npt,nao)
      common /io/ inp, iout
c---------------------------------------------------------------------c
c             do overlap of all bound functions with vlamdas          c
c---------------------------------------------------------------------c
      call apcbtc(ovlmb,vlamda,basis,nolam,npt,nao)
      return
      end





