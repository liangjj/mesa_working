*deck @(#)mcabcx.f	1.1  11/30/90
      subroutine mcabcx(cc,cd,xjk,nbfc,nbfd,naoc,naod,
     $     norbc,norbd,temp,tv,nfav,isym, ncoc,ncod,buf)
c
c***begin prologue     mcabcx
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)mcabcx.f	1.1   11/30/90
c
c***purpose
c
c***description
cc
c----------------------------------------------------------------
c
c     this program performs a limited transformation
c     of the an exchange integral block
c       one of the indices is transformed over the active orbitals.
c       the other index is transformed over all of the orbitals
c
c
c       cc    molecular orbitals
c       cd    molecular orbitals
c       xjk    exchange integral block
c       temp  temporary storage for the transformed fock operator
c       tv    temporary storage .. vector of length  nbf
c       naoc  the number of active orbitals
c       naod  the number of active orbitals
c       norbc the number of orbitals
c       norbd the number of orbitals
c       nbfc  the number of basis functions
c       nbfd  the number of basis functions
c       nfav  output unit
c
c----------------------------------------------------------------
c
c***references
c
c***routines called    (none)
c
c***end prologue       mcabcx
c
      implicit real*8 (a-h,o-z)
c
      real*8 buf(*)
      dimension cc(nbfc,2),cd(nbfd,2),xjk(nbfc,2),temp(2),tv(2)
c
      common / number / zero,pt5,one,two,four,eight
      common / nhexpk / nhexab
c
      call ebc(temp,xjk,cd(1,ncod+1),nbfc,nbfd,naod)
      call ebtc(buf,cc,temp,norbc,nbfc,naod)
c
      if(isym.ne.0) call lnkerr('  isym ne 0 in mcabcx ')
c
      return
      end
